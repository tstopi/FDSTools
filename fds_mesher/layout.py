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

On a face that must stay flush, the lattice is flush at the low end (offset 0)
and the last block on the high end is cut to end at the face. A remainder under
half a block is merged into the previous block instead, so every axis has at
most one block that differs from the nominal size.
"""

import time
from dataclasses import dataclass

import numpy as np

from .parse import get_floats, parse_namelists


@dataclass
class Layout:
    """Nominal block size b=(bi,bj,bk); `edges[a]` are the block boundaries
    (cell indices, may lie outside [0, n] where the lattice overhangs) along
    axis a. `blocks` is the retained mask."""
    b: tuple
    edges: tuple
    blocks: np.ndarray

    @property
    def start(self):
        return tuple(int(e[0]) for e in self.edges)

    @property
    def n_blocks(self):
        return int(np.count_nonzero(self.blocks))

    @property
    def cells_per_block(self):
        return int(np.prod(self.b))

    @property
    def total_cells(self):
        return _cells(self.blocks, self.edges)

    def odd_sizes(self):
        """Per axis: block sizes (cells) that differ from the nominal size."""
        return [sorted({int(d) for d in np.diff(e)} - {self.b[a]})
                for a, e in enumerate(self.edges)]

    def block_boxes(self):
        """(lo, hi): (m,3) integer index corners of each retained block."""
        idx = np.argwhere(self.blocks)
        lo = np.stack([self.edges[a][idx[:, a]] for a in range(3)], axis=1)
        hi = np.stack([self.edges[a][idx[:, a] + 1] for a in range(3)], axis=1)
        return lo, hi


def _cells(mask, edges):
    s = [np.diff(e).astype(np.int64) for e in edges]
    return int(np.einsum("ijk,i,j,k->", mask.astype(np.int64), *s))


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


def axis_edges(n, b, o, low_free, high_free):
    """Block boundaries on one axis for lattice start -o, block size b.

    None if o > 0 on a face that must stay flush. On a flush high face the
    last block is cut at n, or merged into the previous one when the
    remainder is under half a block.
    """
    if o > 0 and not low_free:
        return None
    m = -(-(n + o) // b)
    e = -o + b * np.arange(m + 1)
    if e[-1] > n and not high_free:
        if m > 1 and n - e[-2] < b / 2:
            e = e[:-1]
        e[-1] = n
    return e


def _reduce(a, axis, edges):
    """any() over the blocks given by `edges` along axis."""
    starts = np.clip(edges[:-1], 0, a.shape[axis])
    return np.logical_or.reduceat(a, starts, axis=axis)


def retained_blocks(reachable, edges):
    """Retained-block mask for per-axis block edges."""
    r = reachable
    for a in range(3):
        r = _reduce(r, a, edges[a])
    return r


def smooth_numbers(lo, hi):
    """Integers in [lo, hi] of the form 2^a 3^b 5^c."""
    return [n for n in range(max(lo, 1), hi + 1) if is_smooth(n)]


def candidates(shape, cells_per_mesh, min_block, max_aspect, block=None):
    """Block-size triples (bi,bj,bk), most plausible window first.

    Triples have cell count in [0.5, 1.5] x cells_per_mesh and aspect
    <= max_aspect; the window and aspect limit are relaxed if nothing fits.
    """
    if block is not None:
        for a in range(3):
            if not is_smooth(block[a]):
                raise ValueError(f"block size {block[a]} is not 2^a 3^b 5^c")
        return [tuple(int(x) for x in block)]
    axes = []
    for n in shape:
        v = smooth_numbers(min(min_block, n), n)
        axes.append(v or [max(smooth_numbers(1, n))])
    g = np.stack([m.ravel() for m in np.meshgrid(*axes, indexing="ij")], axis=1)
    cells = g.prod(axis=1)
    aspect = g.max(axis=1) / g.min(axis=1)
    for f in (1, 2, 4, 8, 16):
        ok = (cells >= 0.5 * cells_per_mesh / f) & \
             (cells <= 1.5 * cells_per_mesh * f) & (aspect <= max_aspect * f)
        if ok.any():
            break
    else:
        return [tuple(int(max(a)) for a in axes)]
    g, cells = g[ok], cells[ok]
    keep = np.argsort(np.abs(np.log(cells / cells_per_mesh)))[:1500]
    return [tuple(int(x) for x in r) for r in g[keep]]


def axis_offsets(n, b, step, low_free, high_free):
    """Valid lattice overhangs (cells) -> edges on one axis, as (o, edges)
    pairs. step 0 means flush only."""
    cand = {0} if step <= 0 else set(range(0, b, step)) | {(-n) % b}
    out = [(o, axis_edges(n, b, o, low_free, high_free)) for o in sorted(cand)]
    return [(o, e) for o, e in out if e is not None]


def search(reachable, free, *, cells_per_mesh=300_000, max_aspect=8.0,
           min_block=8, block=None, offset_step=None, time_budget=110.0):
    """Choose block size and lattice offset minimising retained cells.

    Candidates are evaluated hierarchically (any() reduction along x, then y,
    then z) with caching of the partial reductions; a partial result gives a
    lower bound on the retained block count (each non-empty (x,y) column
    needs its own block), which prunes against the incumbent. Score is
    (cells, blocks, aspect). Returns (Layout, info dict).
    """
    shape = reachable.shape
    by = {}
    for bi, bj, bk in candidates(shape, cells_per_mesh, min_block,
                                 max_aspect, block):
        by.setdefault(bi, {}).setdefault(bj, []).append(bk)
    t0 = time.time()
    best, n_eval, truncated = None, 0, False

    def steps(b):
        return max(1, b // 8) if offset_step is None else offset_step

    for bi in sorted(by):
        for oi, ei in axis_offsets(shape[0], bi, steps(bi), *free[0]):
            if time.time() - t0 > time_budget:
                truncated = True
                break
            r1 = _reduce(reachable, 0, ei)
            si = np.diff(ei)
            for bj, ks in by[bi].items():
                for oj, ej in axis_offsets(shape[1], bj, steps(bj), *free[1]):
                    r2 = _reduce(r1, 1, ej)
                    cols = r2.any(axis=2)                    # non-empty columns
                    area = int(np.einsum("ij,i,j->", cols.astype(np.int64),
                                         si, np.diff(ej)))
                    for bk in ks:
                        offs_k = axis_offsets(shape[2], bk, steps(bk), *free[2])
                        lower = area * min(int(np.diff(e).min()) for _, e in offs_k)
                        if best is not None and lower > best[0][0]:
                            continue
                        for ok, ek in offs_k:
                            mask = _reduce(r2, 2, ek)
                            edges = (ei, ej, ek)
                            n_eval += 1
                            asp = max(bi, bj, bk) / min(bi, bj, bk)
                            key = (_cells(mask, edges),
                                   int(np.count_nonzero(mask)), asp)
                            if best is None or key < best[0]:
                                best = (key, (bi, bj, bk), edges, mask)
        if truncated:
            break
    if best is None:
        raise ValueError("no valid block layout found; try another --dx, "
                         "--block or --min-block")
    _, b, edges, mask = best
    lay = Layout(b, edges, mask)
    return lay, {"evaluated": n_eval, "seconds": time.time() - t0,
                 "truncated": truncated}


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
    # (4) IJK consistent with XB (2-3-5 factors are checked in non_smooth)
    for m in range(len(meshes)):
        if (hi[m] - lo[m] != ijk[m]).any():
            errs.append(f"mesh {m + 1}: IJK does not match XB at dx")
    # (2) no two blocks overlap (including outside the voxel domain)
    for m in range(len(meshes) - 1):
        over = (np.maximum(lo[m], lo[m + 1:]) <
                np.minimum(hi[m], hi[m + 1:])).all(axis=1)
        if over.any():
            errs.append(f"mesh {m + 1} overlaps mesh {m + 2 + int(np.argmax(over))}")
            break
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


def non_smooth(layout):
    """Distinct retained-block IJK triples that are not all 2^a 3^b 5^c."""
    lo, hi = layout.block_boxes()
    return sorted({tuple(int(v) for v in r) for r in hi - lo
                   if not all(is_smooth(int(v)) for v in r)})
