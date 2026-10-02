import unittest

import numpy as np

import fds_mesher as fm
from fds_mesher import parse
from fds_mesher.tests import cases
from fds_mesher.tests.test_phase1 import blocks_xb


class Phase2(unittest.TestCase):
    def test_6_l_tunnel_with_shaft(self):
        res = fm.run(cases.l_shaft_tunnel(), 0.2, block=(5, 5, 5))
        self.assertEqual(res.errors, [])
        r = res.reach.reachable
        self.assertEqual(r.sum(), res.reach.n_reachable)
        self.assertGreater(r.sum(), 0)
        xb = blocks_xb(res)
        self.assertLess(len(xb), 0.5 * 10 * 10 * 8)
        self.assertAlmostEqual(xb[:, 5].max(), 8.0)           # shaft top
        self.assertTrue(((xb[:, 0] < 1) & (xb[:, 2] < 3)).any())   # portal leg
        self.assertTrue(((xb[:, 1] > 7) & (xb[:, 3] > 8)).any())   # leg 2 / shaft
        # every carved cell is reachable, nothing else
        solid_free = np.zeros(r.shape, bool)
        solid_free[:40, 5:15, 5:15] = True
        solid_free[30:40, 5:45, 5:15] = True
        solid_free[30:40, 35:45, 5:] = True
        np.testing.assert_array_equal(r, solid_free)

    def test_7_geom_tunnel_matches_shell(self):
        res = fm.run(cases.geom_tunnel(), 0.2, block=(10, 5, 5))
        self.assertEqual(res.errors, [])
        self.assertEqual(res.warnings, [])
        xb = blocks_xb(res)
        self.assertEqual(set(xb[:, 2]), {3.0, 4.0})
        self.assertEqual(set(xb[:, 4]), {3.0, 4.0})
        self.assertEqual(len(xb), 8 * 2 * 2)
        area = np.pi * 0.9 ** 2 / 0.04          # cells per cross-section
        per_slice = res.reach.n_reachable / 80
        # surface voxels eat about one cell of radius
        self.assertTrue(0.6 * area < per_slice < area, (per_slice, area))
        self.assertGreater(res.reach.n_unreachable, 0)
        # stride-3 FACES give the same voxels
        txt3 = cases.head((0, 16, 0, 8, 0, 8), (80, 40, 40)) + cases.geom_text(
            *cases.tube_geom(), stride=3) + cases.open_vent((0, 0, 3.3, 4.7, 3.3, 4.7)) \
            + cases.open_vent((16, 16, 3.3, 4.7, 3.3, 4.7))
        res3 = fm.run(txt3, 0.2, block=(10, 5, 5))
        np.testing.assert_array_equal(res3.solid, res.solid)

    def test_8_mult(self):
        txt = ("&MULT ID='row', DX=2, I_UPPER=4 /\n"
               "&MULT ID='grid', DX=1, DY=3, I_UPPER=2, J_UPPER=1 /\n"
               "&MULT ID='seq', DZ=0.5, N_LOWER=1, N_UPPER=3 /\n"
               "&OBST XB=0,1,0,1,0,1, MULT_ID='row' /\n"
               "&HOLE XB=0,1,0,1,0,1, MULT_ID='grid' /\n"
               "&VENT XB=0,0,0,1,0,1, SURF_ID='OPEN', MULT_ID='seq' /\n")
        m = parse.parse_text(txt)
        self.assertEqual([b.xb[0] for b in m.obsts], [0, 2, 4, 6, 8])
        self.assertEqual(len(m.holes), 6)
        self.assertEqual(sorted({(h.xb[0], h.xb[2]) for h in m.holes}),
                         [(0, 0), (0, 3), (1, 0), (1, 3), (2, 0), (2, 3)])
        self.assertEqual([v.xb[4] for v in m.vents], [0.5, 1.0, 1.5])
        with self.assertRaises(NotImplementedError):
            parse.parse_text("&MULT ID='m', DXB=1,1,0,0,0,0 /\n"
                             "&OBST XB=0,1,0,1,0,1, MULT_ID='m' /\n")
        # expanded OBSTs are voxelised
        res = fm.run("&MESH IJK=50,5,5, XB=0,10,0,1,0,1 /\n" + txt.split("\n")[0]
                     + "\n&OBST XB=0,1,0,1,0,1, MULT_ID='row' /\n"
                     "&VENT MB='YMIN', SURF_ID='OPEN' /\n", 0.2, block=(5, 5, 5))
        self.assertEqual(int(res.solid.sum()), 5 * 5 * 5 * 5)

    def test_seed_surf(self):
        txt = cases.rock_tunnel().replace("SURF_ID='OPEN'", "SURF_ID='SUPPLY'", 1)
        res = fm.run(txt, 0.2, block=(10, 5, 5), min_block=5, seed_surf=["supply"])
        self.assertEqual(len(res.reach.vents), 2)
        with self.assertRaises(ValueError):
            fm.run(txt.replace("SURF_ID='OPEN'", "SURF_ID='X'"), 0.2, block=(10, 5, 5))

    def test_unsupported_geom(self):
        with self.assertRaises(NotImplementedError) as c:
            parse.parse_text("&GEOM ID='ball', SPHERE_ORIGIN=0,0,0, SPHERE_RADIUS=1 /")
        self.assertIn("ball", str(c.exception))
        with self.assertRaises(NotImplementedError):
            parse.parse_text("&GEOM ID='f', BINARY_FILE='a.bingeom' /")
