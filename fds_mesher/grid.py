"""Voxel grid definition: origin, cubic cell size dx, shape, snapping."""

from dataclasses import dataclass

import numpy as np

_EPS = 1e-9


@dataclass
class Grid:
    origin: np.ndarray
    dx: float
    shape: tuple
    bounds_source: str = ""

    @property
    def n_voxels(self):
        return int(np.prod(self.shape, dtype=np.int64))

    def coord(self, i, axis):
        """World coordinate of grid line i on an axis (rounded to 1e-9)."""
        return round(float(self.origin[axis]) + i * self.dx, 9)

    def upper(self):
        return self.origin + self.dx * np.array(self.shape)

    def line_index(self, c, axis):
        """Nearest grid line index to world coordinate c."""
        return int(np.floor((c - self.origin[axis]) / self.dx + 0.5 + _EPS))

    def snap_interval(self, a, b, axis):
        """FDS-like snap of [a, b] to a half-open voxel range (i0, i1).

        Both ends round to the nearest grid line. If they coincide (object
        thinner than a cell), one voxel is kept, on the lower-index side of
        the line, so thin walls still block flow (FDS itself would drop or
        thicken them). The range is clipped to the grid; may be empty.
        """
        i0, i1 = self.line_index(a, axis), self.line_index(b, axis)
        if i1 <= i0:
            i0, i1 = i0 - 1, i0
        return max(i0, 0), min(i1, self.shape[axis])

    def snap_box(self, xb):
        """Voxel ranges ((i0,i1),(j0,j1),(k0,k1)) of an XB, or None if empty."""
        r = tuple(self.snap_interval(xb[2 * a], xb[2 * a + 1], a)
                  for a in range(3))
        return None if any(lo >= hi for lo, hi in r) else r


def domain_bounds(model, bounds=None):
    """Return (xb array of 6, source label)."""
    if bounds is not None:
        return np.array(bounds, float), "--bounds"
    if model.meshes:
        xb = np.array([m.xb for m in model.meshes])
        return _union(xb), "&MESH"
    boxes = [b.xb for b in model.obsts + model.vents if b.xb is not None]
    boxes += [_verts_bbox(g.verts) for g in model.geoms]
    if not boxes:
        raise ValueError("cannot determine domain bounds: no MESH, OBST, "
                         "GEOM or VENT found (use --bounds)")
    return _union(np.array(boxes)), "OBST/GEOM/VENT bbox"


def _verts_bbox(v):
    return np.array([v[:, 0].min(), v[:, 0].max(), v[:, 1].min(),
                     v[:, 1].max(), v[:, 2].min(), v[:, 2].max()])


def _union(xb):
    return np.array([xb[:, 0].min(), xb[:, 1].max(), xb[:, 2].min(),
                     xb[:, 3].max(), xb[:, 4].min(), xb[:, 5].max()])


def build_grid(model, dx, bounds=None, origin=None, max_voxels=2e9):
    """Build the Grid; bounds are snapped outward to the dx lattice.

    Returns (grid, warnings). The grid origin is the snapped domain minimum,
    which lies on the lattice of the user origin (or domain minimum).
    """
    if dx <= 0:
        raise ValueError("dx must be positive")
    xb, source = domain_bounds(model, bounds)
    warnings = []
    org = np.array(origin, float) if origin is not None else xb[0::2].copy()
    lo = np.floor((xb[0::2] - org) / dx + 1e-6).astype(np.int64)
    hi = np.ceil((xb[1::2] - org) / dx - 1e-6).astype(np.int64)
    hi = np.maximum(hi, lo + 1)
    moved = np.abs(np.concatenate([org + lo * dx - xb[0::2],
                                   org + hi * dx - xb[1::2]]))
    if moved.max() > 1e-6:
        warnings.append("domain bounds were snapped outward to the dx grid "
                        f"(max move {moved.max():.4g} m)")
    shape = tuple(int(n) for n in hi - lo)
    n = float(np.prod(shape, dtype=np.float64))
    if n > max_voxels:
        raise ValueError(f"grid has {n:.3g} voxels (~{n / 1e9:.1f} GB as bool), "
                         f"more than --max-voxels {max_voxels:.3g}")
    return Grid(org + lo * dx, dx, shape, source), warnings
