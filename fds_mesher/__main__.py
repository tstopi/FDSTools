"""Command line: python -m fds_mesher INPUT.fds --dx 0.2 --block 10,5,5"""

import argparse
import sys
from pathlib import Path

import numpy as np

from . import run, writer


def _csv(s, n, cast=float):
    v = [cast(x) for x in s.split(",")]
    if len(v) != n:
        raise argparse.ArgumentTypeError(f"expected {n} comma-separated values")
    return v


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="python -m fds_mesher",
        description="Generate &MESH lines tiling only the air reachable from "
                    "OPEN vents.")
    ap.add_argument("input")
    ap.add_argument("--dx", type=float, required=True, help="cell size (m)")
    ap.add_argument("--cells-per-mesh", type=int, default=300_000)
    ap.add_argument("--max-aspect", type=float, default=8.0)
    ap.add_argument("--min-block", type=int, default=8)
    ap.add_argument("--block", type=lambda s: _csv(s, 3, int),
                    help="fixed block size in cells bi,bj,bk (skips the "
                         "block-size search; offsets are still searched)")
    ap.add_argument("--offset-step", type=int,
                    help="lattice offset step in cells (default max(1,b//8); "
                         "0 = zero offset only)")
    ap.add_argument("--bounds", type=lambda s: _csv(s, 6))
    ap.add_argument("--origin", type=lambda s: _csv(s, 3))
    ap.add_argument("--seed-surf", nargs="*", default=[], metavar="SURF_ID")
    ap.add_argument("--leak-fraction", type=float, default=0.9)
    ap.add_argument("--max-voxels", type=float, default=2e9)
    ap.add_argument("--fill-overhang", action=argparse.BooleanOptionalAction,
                    default=True, help="emit solid OBSTs for mesh cells "
                    "outside the voxel domain (default on)")
    ap.add_argument("--mpi", type=int, help="assign MPI_PROCESS to N processes")
    ap.add_argument("--snippet", action="store_true")
    ap.add_argument("-o", "--output")
    ap.add_argument("--dump-voxels", metavar="OUT.npz")
    a = ap.parse_args(argv)

    try:
        text = Path(a.input).read_text()
        res = run(text, a.dx, bounds=a.bounds, origin=a.origin,
                  seed_surf=a.seed_surf, leak_fraction=a.leak_fraction,
                  cells_per_mesh=a.cells_per_mesh, max_aspect=a.max_aspect,
                  min_block=a.min_block, block=a.block,
                  offset_step=a.offset_step, max_voxels=a.max_voxels,
                  fill_overhang=a.fill_overhang, mpi=a.mpi)
    except (OSError, ValueError, NotImplementedError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1

    report = writer.format_report(res, a.dx)
    if a.dump_voxels:
        np.savez_compressed(a.dump_voxels, solid=res.solid,
                            reachable=res.reach.reachable,
                            origin=res.grid.origin, dx=a.dx)
    if res.errors:
        print(report, file=sys.stderr)
        return 2
    if a.snippet and not a.output:
        sys.stdout.write(res.snippet)
        print(report, file=sys.stderr)
        return 0
    out = Path(a.output) if a.output else \
        Path(a.input).with_name(Path(a.input).stem + "_mesh.fds")
    out.write_text(res.snippet if a.snippet else res.output)
    print(report)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
