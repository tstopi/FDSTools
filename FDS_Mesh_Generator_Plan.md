# FDS Automatic Mesh Generator — Implementation Plan

Turns an FDS input file into a set of `&MESH` lines that tile only the **air reachable
from OPEN vents**, using identical, aligned blocks. Target: tunnels, stations,
caverns and ventilation networks.

(Part of the FDSTools monorepo. Plain Python CLI — no Blender.)

---

## 0. Decisions (fixed — do not revisit)

| Topic | Decision |
|---|---|
| Cell size `dx` | **User input** (`--dx`). Cubic cells. The tool never chooses dx. |
| Flood-fill seeds | **OPEN vents only** (`SURF_ID='OPEN'`). Optional `--seed-surf` adds other SURF IDs (e.g. supply/exhaust). The domain box is treated as solid. **No margin. No automatic portal detection.** |
| Geometry style | Both "solid rock" models and "thin-wall shell" models must work. |
| Mesh layout | **Identical blocks** (same IJK, same size) on one lattice. **No merging.** |
| Dependencies | `numpy`, `scipy` only. No trimesh, no VTK. |
| What is optimised | Block size in cells (`bi,bj,bk`) and lattice origin shift. Objective = retained cells. |

### Why these differ from the earlier draft
- Seeding from domain faces with a margin floods all exterior air (all faces are air).
- Optimising dx always picks the coarsest; dx is a physics choice (D*/dx).
- "MPI-friendly" in FDS means: per-mesh I, J, K each factor as 2^a·3^b·5^c (FFT
  pressure solver) and equal cells per mesh (load balance) — not whole-domain counts.
- Merging contradicts equal-sized meshes.
- Coarse voxels (2·dx, mesh/4) leak through thin walls and close doorways.

---

## 1. Core idea

1. Build a boolean voxel grid at **exactly the FDS cell size `dx`**, aligned to a global
   grid origin. 1 voxel = 1 future FDS cell.
2. **Rasterise** obstacles as solid: OBST boxes (snapped like FDS), minus HOLEs, and the
   **surface** voxels of GEOM triangles. No inside/outside tests are needed — the
   interior of a closed solid is simply never reached by the flood fill. This handles
   watertight solids and open thin shells identically.
3. Seed air voxels adjacent to each OPEN vent; 6-connected flood fill
   (`scipy.ndimage.label` on the air mask, keep labels that contain a seed).
4. Choose block size + lattice offset that minimises the cells in blocks containing
   ≥1 reachable voxel. Emit those blocks as `&MESH`.

Correctness argument (put in module docstring): a reachable voxel's air neighbours are
reachable, so every face between a retained and a removed block is solid on the
reachable side → removing blocks never creates a new boundary through reachable air.

---

## 2. Package layout

```
fds_mesher/
├── __init__.py
├── __main__.py        # CLI (argparse) → main()
├── parse.py           # namelist reader + OBST/HOLE/VENT/GEOM/MESH/MULT extraction
├── grid.py            # Grid dataclass (origin, dx, shape), world<->index, snapping
├── voxelize.py        # solid mask from OBST/HOLE/GEOM
├── reach.py           # seeds from vents, flood fill, leak report
├── layout.py          # block-size candidates, offset search, block retention
├── writer.py          # &MESH text, output file, report
└── tests/
    ├── __init__.py
    ├── cases.py       # synthetic FDS input builders
    └── test_*.py      # unittest (run: python -m unittest discover fds_mesher/tests)
```

Keep it to these modules. Use `unittest` (stdlib), runnable without Blender.

---

## 3. Module specs

### 3.1 `parse.py`
- Generic namelist reader: strip `!` comment lines; split on `&NAME ... /` while
  respecting quoted strings (a `/` inside `'...'` must not end a namelist — the regex
  in `geom_to_ast_devc.py` gets this wrong). Return `list[(name, dict[key->raw str])]`.
- Value parsing helpers: floats list, string, int; case-insensitive keys; `T`/`.TRUE.`.
- Extract:
  - `MESH`: `XB`, `IJK`, `MULT_ID` → list of boxes (used for domain bounds only).
  - `MULT`: `ID`, `DX,DY,DZ`, `DXB`, `I_LOWER/I_UPPER` (J, K), `N_LOWER/N_UPPER`.
    Implement expansion for `DX/DY/DZ` + `I/J/K_LOWER/UPPER` and `N_LOWER/N_UPPER`;
    raise a clear `NotImplementedError` for `DXB`/`DX0`-style variants not implemented.
  - `OBST`: `XB`, `MULT_ID`, `SURF_ID`. Skip `REMOVABLE`/`CTRL_ID`/`DEVC_ID`? **No** —
    treat as present (conservative: solid at t=0); log a warning count.
  - `HOLE`: `XB`, `MULT_ID`. Holes with `CTRL_ID`/`DEVC_ID` are treated as **open**
    (air) — conservative for coverage; warn.
  - `VENT`: `XB` or `MB` (`XMIN`…`ZMAX`), `SURF_ID`, `MULT_ID`. `MB` resolves against
    the domain bounds.
  - `GEOM`: `ID`, `VERTS`, `FACES` (stride 3 or 4; reuse the auto-stride logic from
    `geom_to_ast_devc.py`, ported, not imported). `BINARY_FILE`, `ZVALS`, `SPHERE_*`,
    `CYLINDER_*`, `XB`-box GEOMs → `NotImplementedError` with the GEOM ID in the message.
- `&MOVE`, `&INIT`, `&ZONE` ignored.

### 3.2 `grid.py`
- `Grid(origin: (3,), dx: float, shape: (nx,ny,nz))`.
- Domain bounds, in priority: (1) `--bounds x0,x1,y0,y1,z0,z1`; (2) union of existing
  `&MESH` XB; (3) bbox of OBST ∪ GEOM ∪ VENT. Print which was used.
- Origin = `--origin` if given, else domain min. Bounds snapped **outward** to the dx
  grid from the origin. Warn if the snap moved a bound by > 1e-6.
- Snapping of an interval `[a,b]` on an axis: FDS-like — round both ends to the
  nearest grid line. If they coincide (thin object), keep it as **one voxel** on the
  side of the lower index (so thin walls still block flow). Document the deviation.

### 3.3 `voxelize.py`
- `solid = np.zeros(shape, bool)`.
- OBST: snapped index ranges → `solid[i0:i1, j0:j1, k0:k1] = True`.
- HOLE: same, set `False`. Apply after all OBST (FDS order-independent semantics).
  **Do not** let HOLEs carve GEOM.
- GEOM triangles: mark every voxel the triangle passes through. Implementation:
  sample each triangle with barycentric points at spacing ≤ `dx/2` (vectorised per
  triangle batch), convert to indices, set True. Clip to grid. Must handle 1e6
  triangles in reasonable time — batch triangles, avoid Python per-point loops.
- Memory: bool array; estimate `nx*ny*nz` bytes up front and abort with a message if
  > `--max-voxels` (default 2e9).

### 3.4 `reach.py`
- Seeds: for each OPEN vent (and `--seed-surf` vents), snap to a plane on the grid.
  A vent must lie on a face perpendicular to one axis (one XB pair equal); otherwise
  error. Seed = the air voxel layer **inside the domain** adjacent to that plane,
  within the vent's footprint. If the vent is on the domain boundary, inside = toward
  the domain interior; if interior (FDS allows only on obstruction faces / exterior),
  seed both sides' air voxels and warn.
- `air = ~solid`; `labels, n = scipy.ndimage.label(air)` (default 6-connectivity
  structure); `reachable = isin(labels, seed_labels)`.
- Errors/warnings:
  - No OPEN vents → error ("nothing to seed; add OPEN vents or --seed-surf").
  - A vent whose footprint is fully solid → warn with vent index/line.
  - **Leak heuristic**: if reachable volume > `--leak-fraction` (default 0.9) of the
    bounding-domain air, warn "possible leak: reachable air fills the box" — typical of
    a thin-wall model whose walls have a gap, or of whole-face `MB` OPEN vents.
  - Report unreachable air volume (closed cavities removed).

### 3.5 `layout.py`
- **Block-size candidates**: per axis, integers `n` with `n = 2^a 3^b 5^c`,
  `--min-block` ≤ n ≤ axis cell count (default min 8). Triples `(bi,bj,bk)` with
  `bi*bj*bk` in `[0.5, 1.5] × --cells-per-mesh` (default 300 000) and aspect
  `max/min ≤ --max-aspect` (default 8; tunnels need elongated blocks).
- **Offset search**: for each candidate, try lattice offsets `(oi,oj,ok)` in cells,
  `0 ≤ o < b` per axis, on a stride (`--offset-step`, default `max(1, b//8)`) to bound
  cost.
- **Retention**: pad `reachable` so the lattice tiles it; reshape to
  `(Ni,bi,Nj,bj,Nk,bk)` and `.any(axis=(1,3,5))` → retained block mask. Score =
  `retained_blocks * bi*bj*bk`; tie-break fewer blocks, then smaller aspect.
- **Overhang rule (critical)**: the lattice may extend beyond the domain bounds only
  on a domain face that **no reachable voxel touches** and **no OPEN vent lies on**.
  On faces with reachable air or vents, the lattice must be flush with the bound
  (otherwise the vent becomes interior or new air appears). Implement by restricting
  offsets/block sizes per axis-side; candidates that violate it are skipped.
  Overhang cells are outside the voxel grid; FDS will see them as air pockets enclosed
  by solid/exterior. Option `--fill-overhang` (default on) emits `&OBST` boxes filling
  overhang regions of retained blocks so they are solid.
- Pruning for speed: compute a lower bound `ceil(reachable_count / block_cells)`
  blocks; skip candidates whose best possible score already exceeds the incumbent.
- Return chosen `(b, offset, retained_mask)`.

### 3.6 `writer.py`
- `&MESH ID='M_0001', IJK=bi,bj,bk, XB=... /` with coordinates formatted with
  `repr`-safe precision (round to 1e-6 and strip trailing zeros) so adjacent blocks
  share identical bounds as text.
- Optional `--mpi N`: assign `MPI_PROCESS` round-robin along the dominant axis
  (blocks sorted by PCA-free key: the axis with most retained blocks), contiguous runs.
- Output modes:
  - default: write `<input>_mesh.fds` = original text with every `&MESH ... /` line
    commented out (`!MESH-AUTO ` prefix) and the new block (plus filler OBSTs) inserted
    where the first original `&MESH` was, or after `&HEAD` if none.
  - `--snippet`: write only the generated lines to stdout/`-o`.
- Report (stdout):
  ```
  dx: 0.20 m   domain: 2000.0 x 20.0 x 12.0 m   bounds from: &MESH
  voxels: 50.0 M   solid: 31.2 M   reachable air: 14.1 M   unreachable air: 4.7 M
  block: 240 x 40 x 32 cells (48.0 x 8.0 x 6.4 m)   offset: (0,0,0)
  meshes: 58   total cells: 17.8 M   cells/mesh: 307 200
  efficiency (reachable/meshed): 79 %   removed vs. full box: 64 %
  est. memory: ~17.8 GB (1 kB/cell)
  ```

### 3.7 `__main__.py` CLI
```
python -m fds_mesher INPUT.fds --dx 0.2
    [--cells-per-mesh 300000] [--max-aspect 8] [--min-block 8]
    [--bounds x0,x1,y0,y1,z0,z1] [--origin x,y,z]
    [--seed-surf SURF_ID ...] [--leak-fraction 0.9]
    [--offset-step N] [--max-voxels 2e9]
    [--fill-overhang/--no-fill-overhang] [--mpi N]
    [--snippet] [-o OUT]
    [--dump-voxels out.npz]   # reachable/solid masks for debugging
```
Exit codes: 0 ok, 1 input error, 2 validation failure.

---

## 4. Self-validation (run always, before writing)
1. Every reachable voxel lies in exactly one retained block.
2. Blocks don't overlap; all bounds are integer multiples of dx from the origin.
3. Every OPEN vent lies on the exterior boundary of the retained-block union
   (not on an internal block interface, not inside a block).
4. Each block's IJK factors into 2, 3, 5.
Any failure → exit 2 with details.

---

## 5. Tests (synthetic, small dx so they run in < 10 s total)

`cases.py` builds FDS text programmatically. Required cases:

1. **Straight solid-rock tunnel**: box domain fully OBST, tunnel carved with `&HOLE`,
   OPEN vents at both portals. Expect: retained blocks cover the tunnel only; vents on
   exterior; all §4 checks pass.
2. **Thin-wall tunnel**: tunnel walls as 1-cell OBST shells inside an empty box, OPEN
   vents only at the tunnel ends. Expect: air outside the shell is unreachable and
   not meshed.
3. **Thin wall thinner than dx** (e.g. 0.01 m at dx 0.2): must still block (snapping rule).
4. **Leak**: same as 2 with a 1-cell gap in the wall → leak warning fires.
5. **Closed cavity**: sealed room beside the tunnel → excluded, reported as unreachable.
6. **L-shaped tunnel + vertical shaft** with an OPEN vent at the shaft top → meshes
   follow both legs and the shaft.
7. **GEOM tunnel**: a triangulated closed tube (generate a polygonal cylinder shell
   with end caps removed) → equivalent result to case 2.
8. **MULT**: OBST array via `MULT_ID` expands correctly.
9. **Overhang rule**: domain size not divisible by any candidate block → lattice is
   flush at portal faces; filler OBSTs emitted where it overhangs.
10. **Parser**: `/` inside a quoted string, comments, multi-line namelists, `MB=` vents.
11. **Round-trip**: generated file re-parses; MESH count/positions match.

Optional (skip if `fds` not on PATH): run `fds` with `T_END=0` on case 1 output and
assert it starts without mesh errors.

---

## 6. Phases & acceptance

**Phase 1 — MVP** (parse OBST/HOLE/VENT/MESH, voxelise OBST, flood fill, fixed
user-given block size `--block bi,bj,bk` with zero offset, writer, §4 checks).
Accept: tests 1–5, 10, 11 pass.

**Phase 2** — GEOM rasterisation, MULT expansion, `--seed-surf`.
Accept: tests 6–8 pass; 1e6-triangle GEOM voxelises in < 60 s.

**Phase 3** — block-size/offset search, overhang rule + filler OBSTs, `--mpi`, report.
Accept: test 9 passes; search on a 50 M-voxel grid completes in < 2 min.

---

## 7. Out of scope
- Mesh merging, non-uniform/refined meshes, per-mesh dx.
- Choosing dx from HRR.
- Automatic portal detection from domain faces.
- GEOM `BINARY_FILE`, terrain (`ZVALS`), primitive GEOMs, `&MOVE`.
- Time-dependent geometry (devices/controls) beyond the conservative rules above.
- A GUI or Blender integration.

## 8. Known limitations to document in README section
- Voxelisation at dx means openings narrower than ~1 cell may close (as in FDS itself).
- Whole-face `MB=` OPEN vents seed all exterior air touching that face; use `XB` vents
  covering only the portal if the outside should be excluded.
- Removed blocks become FDS exterior (default INERT wall) — correct only because they
  contain no reachable air; the check in §4 enforces this.
