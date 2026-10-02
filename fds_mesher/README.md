# fds_mesher

Turns an FDS input file into `&MESH` lines that tile only the air reachable from
OPEN vents, using identical, aligned blocks. Target: tunnels, stations, caverns
and ventilation networks. Needs `numpy` and `scipy` only.

## Usage

```
python -m fds_mesher INPUT.fds --dx 0.2
    [--cells-per-mesh 300000] [--max-aspect 8] [--min-block 8]
    [--block bi,bj,bk]                 # fixed block size, skips the size search
    [--bounds x0,x1,y0,y1,z0,z1] [--origin x,y,z]
    [--seed-surf SURF_ID ...] [--leak-fraction 0.9]
    [--offset-step N] [--max-voxels 2e9]
    [--fill-overhang | --no-fill-overhang] [--mpi N]
    [--snippet] [-o OUT] [--dump-voxels out.npz]
```

`--dx` is the FDS cell size (cubic cells) and is never chosen by the tool. The
default output is `<input>_mesh.fds`: the original text with every `&MESH`
commented out (`!MESH-AUTO ` prefix) and the generated meshes inserted in its
place. `--snippet` writes only the generated lines (stdout, or `-o`). A report
is printed on every run. The block size search needs a seed vent: the `OPEN`
vents (plus any `--seed-surf` surfaces) are the only flood-fill seeds; the
domain box is treated as solid.

Exit codes: 0 ok, 1 input error, 2 self-validation failure. The self-validation
(every reachable voxel in exactly one block, no overlap, bounds on the dx grid,
vents on the exterior of the retained meshes, IJK factors 2/3/5) runs on every
invocation, before writing.

`overhang` in the report is how many cells the block lattice starts before the
domain minimum on each axis. The lattice may overhang only on domain faces with
no OPEN vent and no reachable air; overhanging mesh cells are filled with
`&OBST` (disable with `--no-fill-overhang`). `--mpi N` assigns `MPI_PROCESS`
in N contiguous runs along the axis with the most blocks.

Tests: `python -m unittest discover fds_mesher/tests -v`.

## Known limitations

- Voxelisation is at dx, so openings narrower than about one cell may close (as
  in FDS itself). Objects thinner than a cell become one voxel (lower-index
  side) so thin walls still block flow.
- Whole-face `MB=` OPEN vents seed all exterior air touching that face; use `XB`
  vents covering only the portal if the outside should be excluded. A vent that
  is not on the domain boundary seeds air on both sides and warns.
- Removed blocks become FDS exterior (default INERT wall). That is correct only
  because they contain no reachable air, which the validation enforces.
  Original `&OBST` etc. outside the retained meshes are left in the file; FDS
  ignores them. References to the old mesh IDs (`MESH_ID=`) are not rewritten.
- Removable/controlled OBST are treated as present, controlled HOLEs as open.
- `&MULT` supports DX/DY/DZ with I/J/K_LOWER/UPPER and N_LOWER/N_UPPER; DXB,
  DX0-style and `*_SKIP` raise `NotImplementedError`.
- GEOM: only `VERTS`/`FACES` triangle meshes (FACES stride 3 or 4). `BINARY_FILE`,
  `ZVALS`, primitive GEOMs and `&MOVE` are not supported. HOLEs do not carve GEOM.
- Leak detection is a heuristic (reachable air above `--leak-fraction` of the box).
- The whole grid is held as bool arrays (about 1 byte/voxel, plus ~5 bytes/voxel
  during the flood fill); `--max-voxels` aborts above 2e9.
