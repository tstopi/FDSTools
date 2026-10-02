"""Solid voxel mask from OBST boxes and HOLEs (GEOM added in phase 2)."""

import numpy as np


def voxelize(model, grid):
    """Return the boolean solid mask (True = solid) at cell size grid.dx.

    OBSTs are set first, then HOLEs clear them (order independent, as in FDS).
    """
    solid = np.zeros(grid.shape, bool)
    for b in model.obsts:
        _fill(solid, grid, b.xb, True)
    for b in model.holes:
        _fill(solid, grid, b.xb, False)
    return solid


def _fill(solid, grid, xb, value):
    r = grid.snap_box(xb)
    if r is not None:
        solid[r[0][0]:r[0][1], r[1][0]:r[1][1], r[2][0]:r[2][1]] = value
