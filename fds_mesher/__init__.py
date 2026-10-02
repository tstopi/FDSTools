"""Automatic FDS &MESH generator: tiles only the air reachable from OPEN vents."""

from dataclasses import dataclass, field

from . import grid as _grid, layout as _layout, parse, reach, voxelize, writer


@dataclass
class Result:
    model: object
    grid: object
    solid: object
    reach: object
    layout: object
    snippet: str
    output: str
    warnings: list = field(default_factory=list)
    errors: list = field(default_factory=list)


def run(text, dx, *, bounds=None, origin=None, seed_surf=(), leak_fraction=0.9,
        block=None, max_voxels=2e9, mpi=None):
    """Full pipeline on FDS input text. Raises ValueError on input problems."""
    model = parse.parse_text(text)
    warns = list(model.warnings)
    grid, w = _grid.build_grid(model, dx, bounds, origin, max_voxels)
    warns += w
    solid = voxelize.voxelize(model, grid)
    vents, w = reach.resolve_vents(model, grid, seed_surf)
    warns += w
    rr = reach.flood(solid, vents, leak_fraction)
    warns += rr.warnings
    if rr.n_reachable == 0:
        raise ValueError("no reachable air: every vent footprint is solid")
    if block is None:
        raise ValueError("block size required (--block)")
    free = _layout.free_faces(rr.reachable, vents)
    lay = _layout.fixed_layout(rr.reachable, free, block)
    snippet = writer.build_snippet(grid, lay, mpi)
    errors = _layout.validate(snippet, grid, rr.reachable, vents)
    output = writer.merge_into_text(model, snippet)
    return Result(model, grid, solid, rr, lay, snippet, output, warns, errors)
