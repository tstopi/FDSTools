import unittest

import numpy as np

import fds_mesher as fm
from fds_mesher import parse
from fds_mesher.tests import cases


def blocks_xb(res):
    """(m,6) block bounds from the generated snippet."""
    return np.array([parse.get_floats(nl.params["XB"])
                     for nl in parse.parse_namelists(res.snippet)
                     if nl.name == "MESH"])


class SolidRock(unittest.TestCase):
    def test_1_straight_rock_tunnel(self):
        res = fm.run(cases.rock_tunnel(), 0.2, block=(10, 5, 5))
        self.assertEqual(res.errors, [])
        xb = blocks_xb(res)
        # 8 x 2 x 2 blocks cover the tunnel only (of 8 x 4 x 4 in the box)
        self.assertEqual(len(xb), 32)
        self.assertGreaterEqual(xb[:, 2].min(), 1 - 1e-9)
        self.assertLessEqual(xb[:, 3].max(), 3 + 1e-9)
        self.assertEqual(res.warnings, [])
        # vents on the exterior of the union
        self.assertAlmostEqual(xb[:, 0].min(), 0)
        self.assertAlmostEqual(xb[:, 1].max(), 16)

    def test_5_closed_cavity(self):
        txt = cases.rock_tunnel() + cases.hole((8, 10, 3.2, 3.8, 1, 3))
        res = fm.run(txt, 0.2, block=(10, 5, 5))
        self.assertEqual(res.errors, [])
        self.assertEqual(res.reach.n_unreachable, 10 * 3 * 10)
        self.assertEqual(len(blocks_xb(res)), 32)   # cavity not meshed


class ThinWall(unittest.TestCase):
    def test_2_thin_wall_tunnel(self):
        res = fm.run(cases.shell_tunnel(), 0.2, block=(10, 5, 5))
        self.assertEqual(res.errors, [])
        self.assertEqual(res.warnings, [])
        self.assertGreater(res.reach.n_unreachable, 0)
        xb = blocks_xb(res)
        # interior is y,z in 3..5 -> blocks [3,4],[4,5] only
        self.assertEqual(set(xb[:, 2]), {3.0, 4.0})
        self.assertEqual(set(xb[:, 4]), {3.0, 4.0})
        self.assertEqual(len(xb), 8 * 2 * 2)

    def test_3_wall_thinner_than_dx(self):
        res = fm.run(cases.shell_tunnel(t=0.01), 0.2, block=(10, 5, 5))
        self.assertEqual(res.errors, [])
        self.assertFalse(any("leak" in w for w in res.warnings))
        r = res.reach.reachable
        self.assertLessEqual(r.sum(), 80 * 12 * 12)
        self.assertGreater(r.sum(), 0)
        ys, zs = np.nonzero(r.any(axis=0))[0], np.nonzero(r.any(axis=1))[0]
        self.assertTrue(ys.min() >= 14 and ys.max() <= 25)

    def test_4_leak_warns(self):
        res = fm.run(cases.shell_tunnel(gap=True), 0.2, block=(10, 5, 5))
        self.assertTrue(any("leak" in w for w in res.warnings), res.warnings)


class Parser(unittest.TestCase):
    def test_10_parser(self):
        txt = ("! &OBST XB=0,1,0,1,0,1 / commented out\n"
               "&HEAD CHID='a/b', TITLE='x / y ! not a comment' /\n"
               "&MESH IJK=10,10,10,\n      XB=0,2, 0,2,\n 0,2 / ! trailing\n"
               "&VENT MB='XMIN', SURF_ID='OPEN' /\n"
               "&VENT MB='ZMAX', SURF_ID='open' /\n"
               "Free text with an apostrophe: don't\n"
               "&OBST XB=0,1,0,1,0,1, SURF_ID='A=B' /\n")
        nls = parse.parse_namelists(txt)
        self.assertEqual([n.name for n in nls],
                         ["HEAD", "MESH", "VENT", "VENT", "OBST"])
        self.assertEqual(parse.get_str(nls[0].params["CHID"]), "a/b")
        self.assertIn("/", parse.get_str(nls[0].params["TITLE"]))
        self.assertEqual((nls[1].line, nls[1].line_end), (3, 5))
        self.assertEqual(parse.get_str(nls[4].params["SURF_ID"]), "A=B")
        m = parse.parse_text(txt)
        self.assertEqual([v.mb for v in m.vents], ["XMIN", "ZMAX"])
        self.assertEqual(len(m.obsts), 1)   # the commented OBST is ignored
        self.assertEqual(m.meshes[0].xb.tolist(), [0, 2, 0, 2, 0, 2])

    def test_11_round_trip(self):
        res = fm.run(cases.rock_tunnel(), 0.2, block=(10, 5, 5))
        again = parse.parse_text(res.output)
        self.assertEqual(len(again.meshes), len(res.layout.block_starts()))
        self.assertEqual(sum(1 for ln in res.output.splitlines()
                             if ln.startswith("!MESH-AUTO")), 1)
        got = sorted(tuple(m.xb) for m in again.meshes)
        starts = res.layout.block_starts()
        g = res.grid
        want = sorted(tuple(g.coord(s[a] + o, a) for a in range(3)
                            for o in (0, res.layout.b[a])) for s in starts)
        want = [(w[0], w[1], w[2], w[3], w[4], w[5]) for w in want]
        self.assertEqual(len(got), len(want))
        for a, b in zip(got, sorted(want)):
            np.testing.assert_allclose(a, b)
