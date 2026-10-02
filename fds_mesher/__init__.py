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
    info: dict = field(default_factory=dict)


def run(text, dx, *, bounds=None, origin=None, seed_surf=(), leak_fraction=0.9,
        cells_per_mesh=300_000, max_aspect=8.0, min_block=8, block=None,
        offset_step=None, max_voxels=2e9, fill_overhang=True, mpi=None,
        time_budget=110.0):
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
    free = _layout.free_faces(rr.reachable, vents)
    lay, info = _layout.search(
        rr.reachable, free, cells_per_mesh=cells_per_mesh,
        max_aspect=max_aspect, min_block=min_block, block=block,
        offset_step=offset_step, time_budget=time_budget)
    if info["truncated"]:
        warns.append(f"layout search hit the {time_budget:g} s budget; "
                     "result may not be optimal")
    snippet = writer.build_snippet(grid, lay, mpi, fill_overhang)
    errors = _layout.validate(snippet, grid, rr.reachable, vents)
    output = writer.merge_into_text(model, snippet)
    return Result(model, grid, solid, rr, lay, snippet, output, warns, errors, info)
