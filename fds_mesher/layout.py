"""Mesh block layout on the voxel grid, and the output self-validation.

Correctness argument. Reachable air is closed under face adjacency by
construction: every air neighbour of a reachable voxel is reachable. So at any
face between a retained block and a removed block, the voxel on the retained
side is either solid or unreachable air that is separated from reachable air
by solid, and the removed side holds no reachable air. Removing blocks that
contain no reachable voxel therefore never turns reachable air into a domain
boundary. The only other boundary that matters is where OPEN vents sit; they
must stay on the outside of the retained union, which is why the lattice may
overhang the domain only on faces that carry neither vents nor reachable air.
"""

from dataclasses import dataclass

import numpy as np

from .parse import get_floats, parse_namelists


@dataclass
class Layout:
    """Identical blocks b=(bi,bj,bk) on a lattice starting at index `start`
    (<= 0 per axis; negative means overhang). `blocks` is the retained mask."""
    b: tuple
    start: tuple
    blocks: np.ndarray

    @property
    def n_blocks(self):
        return int(np.count_nonzero(self.blocks))

    @property
    def cells_per_block(self):
        return int(np.prod(self.b))

    def block_starts(self):
        """(m,3) integer index of each retained block's low corner."""
        idx = np.argwhere(self.blocks)
        return np.array(self.start) + idx * np.array(self.b)


def is_smooth(n):
    for p in (2, 3, 5):
        while n % p == 0:
            n //= p
    return n == 1


def free_faces(reachable, vents):
    """free[axis][side]: True if no reachable voxel touches that domain face
    and no vent lies on it (the lattice may overhang there)."""
    free = [[True, True] for _ in range(3)]
    for a in range(3):
        for side, idx in ((0, 0), (1, reachable.shape[a] - 1)):
            if reachable.take(idx, axis=a).any():
                free[a][side] = False
    for v in vents:
        if v.side == "lo":
            free[v.axis][0] = False
        elif v.side == "hi":
            free[v.axis][1] = False
    return free


def _reduce(a, axis, b, o):
    """any() over blocks of size b along axis, lattice start at -o."""
    n = a.shape[axis]

    def sl(s, e):
        idx = [slice(None)] * 3
        idx[axis] = slice(s, e)
        return a[tuple(idx)]

    parts, pos = [], 0
    if o > 0:
        pos = min(b - o, n)
        parts.append(sl(0, pos).any(axis=axis, keepdims=True))
    m = (n - pos) // b
    if m > 0:
        shp = list(a.shape)
        shp[axis:axis + 1] = [m, b]
        parts.append(sl(pos, pos + m * b).reshape(shp).any(axis=axis + 1))
        pos += m * b
    if pos < n:
        parts.append(sl(pos, n).any(axis=axis, keepdims=True))
    return parts[0] if len(parts) == 1 else np.concatenate(parts, axis=axis)


def retained_blocks(reachable, b, offset):
    """Retained-block mask for block size b and lattice offsets (cells)."""
    r = reachable
    for a in range(3):
        r = _reduce(r, a, b[a], offset[a])
    return r


def axis_ok(n, b, o, low_free, high_free):
    """Overhang rule for one axis: lattice start -o, block size b."""
    if o > 0 and not low_free:
        return False
    high = -(-(n + o) // b) * b - n - o
    return high == 0 or high_free


def fixed_layout(reachable, free, block, offset=(0, 0, 0)):
    """Layout for a user-given block size (phase 1: zero offset)."""
    for a in range(3):
        if not is_smooth(block[a]):
            raise ValueError(f"block size {block[a]} is not 2^a 3^b 5^c")
        if not axis_ok(reachable.shape[a], block[a], offset[a],
                       *free[a]):
            raise ValueError(
                f"block {tuple(block)} does not tile axis {'xyz'[a]} "
                f"({reachable.shape[a]} cells) flush to a face that carries "
                "a vent or reachable air")
    return Layout(tuple(block), tuple(-o for o in offset),
                  retained_blocks(reachable, block, offset))


# ============================================================
# SELF-VALIDATION
# ============================================================

def validate(text, grid, reachable, vents):
    """Check the generated text against the voxel result. Returns error list."""
    errs = []
    meshes = [nl for nl in parse_namelists(text) if nl.name == "MESH"]
    if not meshes:
        return ["no &MESH generated"]
    org, dx = grid.origin, grid.dx
    lo = np.empty((len(meshes), 3), np.int64)
    hi = np.empty_like(lo)
    ijk = np.empty_like(lo)
    for m, nl in enumerate(meshes):
        xb = np.array(get_floats(nl.params["XB"])).reshape(3, 2)
        ijk[m] = [int(v) for v in get_floats(nl.params["IJK"])]
        f = (xb - org[:, None]) / dx
        if np.abs(f - np.round(f)).max() > 1e-4:
            errs.append(f"{nl.params.get('ID')}: bounds are not multiples of dx "
                        "from the origin")
        lo[m], hi[m] = np.round(f[:, 0]), np.round(f[:, 1])
    # (4) IJK factorisation and consistency with XB
    for m in range(len(meshes)):
        if not all(is_smooth(int(n)) for n in ijk[m]):
            errs.append(f"mesh {m + 1}: IJK {tuple(ijk[m])} not 2^a 3^b 5^c")
        if (hi[m] - lo[m] != ijk[m]).any():
            errs.append(f"mesh {m + 1}: IJK does not match XB at dx")
    # (2) identical blocks on one lattice, so no overlap unless duplicated
    if (ijk != ijk[0]).any():
        errs.append("blocks differ in size")
    elif ((lo - lo[0]) % ijk[0]).any():
        errs.append("blocks are not on a common lattice")
    elif len(np.unique(lo, axis=0)) != len(lo):
        errs.append("duplicate (overlapping) blocks")
    if errs:
        return errs
    # (1) every reachable voxel in exactly one block
    cov = np.zeros(grid.shape, np.uint8)
    for m in range(len(meshes)):
        l, h = np.maximum(lo[m], 0), np.minimum(hi[m], grid.shape)
        if (l < h).all():
            cov[l[0]:h[0], l[1]:h[1], l[2]:h[2]] += 1
    if (cov[reachable] != 1).any():
        errs.append(f"{int(np.count_nonzero(cov[reachable] != 1))} reachable "
                    "voxels are not in exactly one block")
    if cov.max() > 1:
        errs.append("blocks overlap")
    # (3) boundary vents on the exterior of the retained union
    for v in vents:
        if v.side == "interior":
            continue
        r = [list(x) for x in v.ranges]
        r[v.axis] = [-1, 0] if v.side == "lo" else \
            [grid.shape[v.axis], grid.shape[v.axis] + 1]
        r = np.array(r)
        inter = (np.maximum(lo, r[:, 0]) < np.minimum(hi, r[:, 1])).all(axis=1)
        if inter.any():
            errs.append(f"vent at line {v.line} is not on the exterior of the "
                        "retained meshes (a block extends beyond it)")
    return errs
