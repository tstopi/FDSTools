"""Vent seeds, flood fill of reachable air, and leak report.

Air is the complement of the solid mask. Seeds are the air voxels next to each
OPEN vent (plus any --seed-surf vents); air connected to a seed through face
neighbours (6-connectivity) is "reachable". A rasterised thin surface that is
26-connected cannot be crossed by a 6-connected flood, so shells and solids
need no inside/outside test.
"""

from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage

_MB = {"XMIN": (0, 0), "XMAX": (0, 1), "YMIN": (1, 0), "YMAX": (1, 1),
       "ZMIN": (2, 0), "ZMAX": (2, 1)}


@dataclass
class VentFace:
    """A vent snapped to a grid plane. `ranges` are voxel ranges per axis;
    on `axis` it is the half-open range of the air-side seed layer(s)."""
    axis: int
    plane: int          # grid line index along axis
    ranges: tuple       # ((i0,i1),(j0,j1),(k0,k1)) of the seed layers
    side: str           # 'lo' / 'hi' (domain boundary) or 'interior'
    line: int
    surf: str


@dataclass
class ReachResult:
    reachable: np.ndarray
    vents: list
    n_solid: int
    n_reachable: int
    n_unreachable: int
    warnings: list = field(default_factory=list)


def resolve_vents(model, grid, seed_surfs=()):
    """Snap seeding vents to grid planes. Returns (list[VentFace], warnings)."""
    wanted = {"OPEN"} | {s.upper() for s in seed_surfs}
    faces, warns = [], []
    shape = grid.shape
    for v in model.vents:
        if (v.surf or "").upper() not in wanted:
            continue
        if v.mb is not None:
            if v.mb not in _MB:
                raise ValueError(f"&VENT line {v.line}: bad MB '{v.mb}'")
            axis, hi = _MB[v.mb]
            plane = shape[axis] if hi else 0
            fp = [(0, shape[a]) for a in range(3)]
        else:
            flat = [a for a in range(3) if abs(v.xb[2 * a + 1] - v.xb[2 * a]) < 1e-6]
            if len(flat) != 1:
                raise ValueError(f"&VENT line {v.line}: XB must be planar "
                                 "(exactly one equal pair)")
            axis = flat[0]
            plane = grid.line_index(v.xb[2 * axis], axis)
            if plane < 0 or plane > shape[axis]:
                warns.append(f"vent at line {v.line} lies outside the grid; ignored")
                continue
            fp = [grid.snap_interval(v.xb[2 * a], v.xb[2 * a + 1], a)
                  if a != axis else (0, shape[a]) for a in range(3)]
            if any(lo >= hi_ for lo, hi_ in fp):
                warns.append(f"vent at line {v.line} has an empty footprint; ignored")
                continue
        if plane == 0:
            side, layer = "lo", (0, 1)
        elif plane == shape[axis]:
            side, layer = "hi", (shape[axis] - 1, shape[axis])
        else:
            side, layer = "interior", (plane - 1, plane + 1)
            warns.append(f"vent at line {v.line} is not on the domain boundary; "
                         "seeding air on both sides")
        fp[axis] = layer
        faces.append(VentFace(axis, plane, tuple(fp), side, v.line, v.surf))
    return faces, warns


def flood(solid, vents, leak_fraction=0.9):
    """Flood-fill air from the vent seeds. Returns a ReachResult."""
    if not vents:
        raise ValueError("no OPEN vents to seed from; add OPEN vents or "
                         "--seed-surf")
    warns = []
    labels, n = ndimage.label(~solid)
    keep = np.zeros(n + 1, bool)
    for v in vents:
        sl = tuple(slice(lo, hi) for lo, hi in v.ranges)
        seeds = labels[sl]
        seeds = seeds[seeds > 0]
        if seeds.size == 0:
            warns.append(f"vent at line {v.line}: footprint is fully solid")
        keep[np.unique(seeds)] = True
    reachable = keep[labels]
    del labels
    n_solid = int(np.count_nonzero(solid))
    n_reach = int(np.count_nonzero(reachable))
    total = solid.size
    if n_reach > leak_fraction * total:
        warns.append(f"possible leak: reachable air fills {100 * n_reach / total:.0f} %"
                     " of the box (gap in a thin wall, or whole-face MB OPEN vent?)")
    return ReachResult(reachable, vents, n_solid, n_reach,
                       total - n_solid - n_reach, warns)
