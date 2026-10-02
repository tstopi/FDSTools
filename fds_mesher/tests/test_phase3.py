import itertools
import unittest

import numpy as np

import fds_mesher as fm
from fds_mesher import layout, parse
from fds_mesher.tests import cases
from fds_mesher.tests.test_phase1 import blocks_xb


class Phase3(unittest.TestCase):
    def test_9_overhang_rule(self):
        kw = dict(block=(10, 10, 10))
        res = fm.run(cases.odd_rock_tunnel(), 0.2, **kw)
        self.assertEqual(res.errors, [])
        xb = blocks_xb(res)
        # flush at the portal faces
        self.assertAlmostEqual(xb[:, 0].min(), 0)
        self.assertAlmostEqual(xb[:, 1].max(), 16)
        # y (37 cells) cannot be tiled exactly: lattice overhangs the domain
        self.assertTrue(xb[:, 2].min() < -1e-9 or xb[:, 3].max() > 7.4 + 1e-9)
        fill = [nl for nl in parse.parse_namelists(res.snippet) if nl.name == "OBST"]
        self.assertGreater(len(fill), 0)
        # fillers lie inside the meshed region and outside the voxel domain
        for nl in fill:
            f = np.array(parse.get_floats(nl.params["XB"]))
            outside = (f[2] < -1e-9 or f[3] > 7.4 + 1e-9 or f[4] < -1e-9
                       or f[5] > 7 + 1e-9 or f[0] < -1e-9 or f[1] > 16 + 1e-9)
            self.assertTrue(outside)
        self.assertEqual(res.warnings, [])
        res2 = fm.run(cases.odd_rock_tunnel(), 0.2, fill_overhang=False, **kw)
        self.assertEqual(res2.errors, [])
        self.assertNotIn("&OBST", res2.snippet)
        # total filled volume equals overhang volume of the retained blocks
        vol = sum(np.prod(np.diff(np.array(parse.get_floats(nl.params["XB"])
                                           ).reshape(3, 2), axis=1))
                  for nl in fill)
        mesh_vol = len(xb) * np.prod(np.array(res.layout.b)) * 0.2 ** 3
        grid_box = np.array([0, 16, 0, 7.4, 0, 7]).reshape(3, 2)
        ov = 0.0
        for r in xb:
            r = r.reshape(3, 2)
            inter = np.prod(np.maximum(0, np.minimum(r[:, 1], grid_box[:, 1])
                                       - np.maximum(r[:, 0], grid_box[:, 0])))
            ov += np.prod(r[:, 1] - r[:, 0]) - inter
        self.assertAlmostEqual(vol, ov, places=4)

    def test_free_faces_force_flush(self):
        r = np.zeros((40, 37, 35), bool)
        r[:, 10:20, 10:20] = True
        free = layout.free_faces(r, [])
        self.assertEqual(free, [[False, False], [True, True], [True, True]])

    def test_search_matches_brute_force(self):
        rng = np.random.default_rng(1)
        r = np.zeros((36, 30, 25), bool)
        r[:, 7:13, 4:9] = True
        r[20:30, 7:28, 4:9] = True
        free = [[False, False], [True, True], [True, True]]
        cands = layout.candidates(r.shape, 300, 4, 8)
        best = None
        for b in cands:
            offs = [layout.axis_offsets(r.shape[a], b[a], max(1, b[a] // 8), *free[a])
                    for a in range(3)]
            for oe in itertools.product(*offs):
                edges = tuple(e for _, e in oe)
                mask = layout.retained_blocks(r, edges)
                n = int(mask.sum())
                key = (layout._cells(mask, edges), n, max(b) / min(b))
                best = key if best is None or key < best else best
        lay, _ = layout.search(r, free, cells_per_mesh=300, min_block=4)
        n = lay.n_blocks
        self.assertEqual((lay.total_cells, n,
                          max(lay.b) / min(lay.b)), best)
        self.assertTrue(lay.blocks.any())

    def test_search_without_block(self):
        res = fm.run(cases.rock_tunnel(), 0.2, cells_per_mesh=500, min_block=5)
        self.assertEqual(res.errors, [])
        lay = res.layout
        self.assertTrue(all(layout.is_smooth(n) for n in lay.b))

    def test_mpi(self):
        res = fm.run(cases.rock_tunnel(), 0.2, block=(10, 5, 5), mpi=4)
        procs = [int(nl.params["MPI_PROCESS"]) for nl in
                 parse.parse_namelists(res.snippet) if nl.name == "MESH"]
        self.assertEqual(sorted(set(procs)), [0, 1, 2, 3])
        self.assertEqual(procs, sorted(procs))
        self.assertEqual(procs.count(0), 8)

    def test_validation_catches_bad_input(self):
        res = fm.run(cases.rock_tunnel(), 0.2, block=(10, 5, 5))
        bad = res.snippet.replace("IJK=10,5,5", "IJK=10,5,7", 1)
        errs = layout.validate(bad, res.grid, res.reach.reachable, res.reach.vents)
        self.assertTrue(errs)
        lines = res.snippet.splitlines()
        short = "\n".join(lines[:5] + lines[6:])      # drop one block
        errs = layout.validate(short, res.grid, res.reach.reachable, res.reach.vents)
        self.assertTrue(any("not in exactly one block" in e for e in errs))


class EndBlocks(unittest.TestCase):
    def test_axis_edges(self):
        e = layout.axis_edges
        np.testing.assert_array_equal(e(40, 10, 0, False, False), [0, 10, 20, 30, 40])
        # remainder 3 < 5: merged into the previous block
        np.testing.assert_array_equal(e(83, 10, 0, False, False),
                                      [0, 10, 20, 30, 40, 50, 60, 70, 83])
        # remainder 7 >= 5: short last block
        self.assertEqual(list(e(87, 10, 0, False, False))[-2:], [80, 87])
        # free high face: overhang, no cut
        self.assertEqual(e(83, 10, 0, False, True)[-1], 90)
        self.assertIsNone(e(83, 10, 3, False, False))
        self.assertEqual(list(e(4, 10, 0, False, False)), [0, 4])

    def test_portals_on_indivisible_length(self):
        # 16.6 m at dx 0.2 = 83 cells (prime): no block size tiles x flush
        for length, odd in ((16.6, 13), (17.4, 7)):
            txt = cases.rock_tunnel(L=length)
            res = fm.run(txt, 0.2, block=(10, 5, 5))
            self.assertEqual(res.errors, [])
            xb = blocks_xb(res)
            self.assertAlmostEqual(xb[:, 0].min(), 0)
            self.assertAlmostEqual(xb[:, 1].max(), length)
            self.assertEqual(res.layout.odd_sizes()[0], [odd])
            if odd == 13:
                self.assertTrue(any("not factorable" in w for w in res.warnings))
            self.assertIn("end blocks", fm.writer.format_report(res, 0.2))
        res = fm.run(cases.rock_tunnel(L=16.6), 0.2, cells_per_mesh=500,
                     min_block=5)
        self.assertEqual(res.errors, [])
