"""Solid voxel mask from OBST boxes, HOLEs and GEOM triangle surfaces."""

import numpy as np

# Triangles with an edge longer than this many sample steps are split first,
# which bounds the per-triangle sample lattice.
_MAX_N = 48
_CHUNK_POINTS = 4_000_000


def voxelize(model, grid):
    """Return the boolean solid mask (True = solid) at cell size grid.dx.

    OBSTs are set first, then HOLEs clear them (order independent, as in FDS).
    GEOM triangles are then added; HOLEs do not carve GEOM.
    """
    solid = np.zeros(grid.shape, bool)
    for b in model.obsts:
        _fill(solid, grid, b.xb, True)
    for b in model.holes:
        _fill(solid, grid, b.xb, False)
    for g in model.geoms:
        raster_triangles(solid, grid, g.verts[g.faces])
    return solid


def _fill(solid, grid, xb, value):
    r = grid.snap_box(xb)
    if r is not None:
        solid[r[0][0]:r[0][1], r[1][0]:r[1][1], r[2][0]:r[2][1]] = value


def _bary(n):
    """Barycentric weights (m,3) of a triangle subdivided n times per edge."""
    i, j = np.meshgrid(np.arange(n + 1), np.arange(n + 1), indexing="ij")
    keep = (i + j) <= n
    i, j = i[keep], j[keep]
    return np.stack([(n - i - j), i, j], axis=1) / n


def _split_long(tri, max_len):
    """Bisect the longest edge of triangles until all edges <= max_len."""
    while True:
        e = np.stack([np.linalg.norm(tri[:, (k + 1) % 3] - tri[:, k], axis=1)
                      for k in range(3)], axis=1)
        long_ = e.max(axis=1) > max_len
        if not long_.any():
            return tri
        t, k = tri[long_], e[long_].argmax(axis=1)
        r = np.arange(len(t))
        a, b, c = t[r, k], t[r, (k + 1) % 3], t[r, (k + 2) % 3]
        m = 0.5 * (a + b)
        tri = np.concatenate([tri[~long_],
                              np.stack([a, m, c], axis=1),
                              np.stack([m, b, c], axis=1)])


def raster_triangles(solid, grid, tri):
    """Mark every voxel passing a sample of the triangles (nt,3,3) as solid.

    Triangles are sampled on a barycentric lattice with spacing <= dx/2, so
    consecutive samples fall in the same or a face/edge/corner-adjacent voxel
    and the marked surface is 26-connected (watertight for a 6-connected
    flood fill). Vectorised over batches of triangles of equal lattice size.
    """
    dx, org = grid.dx, grid.origin
    shape = np.array(grid.shape)
    h = dx / 2
    # drop triangles whose bounding box misses the grid
    lo_b, hi_b = tri.min(axis=1), tri.max(axis=1)
    inside = ((hi_b >= org) & (lo_b <= org + shape * dx)).all(axis=1)
    tri = _split_long(tri[inside], _MAX_N * h)
    if not len(tri):
        return
    e = np.stack([np.linalg.norm(tri[:, (k + 1) % 3] - tri[:, k], axis=1)
                  for k in range(3)], axis=1).max(axis=1)
    n = np.maximum(1, np.ceil(e / h - 1e-9)).astype(np.int64)
    order = np.argsort(n, kind="stable")
    tri, n = tri[order], n[order]
    flat = solid.reshape(-1)
    ny, nz = int(shape[1]), int(shape[2])
    bounds = np.flatnonzero(np.diff(n)) + 1
    for lo, hi in zip(np.r_[0, bounds], np.r_[bounds, len(n)]):
        w = _bary(int(n[lo]))
        step = max(1, _CHUNK_POINTS // len(w))
        for s in range(lo, hi, step):
            t = tri[s:min(s + step, hi)]
            p = np.einsum("tcd,mc->tmd", t, w)
            ijk = np.floor((p.reshape(-1, 3) - org) / dx).astype(np.int64)
            ok = ((ijk >= 0) & (ijk < shape)).all(axis=1)
            ijk = ijk[ok]
            flat[(ijk[:, 0] * ny + ijk[:, 1]) * nz + ijk[:, 2]] = True
