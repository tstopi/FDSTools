# -*- coding: utf-8 -*-
"""Phases 4-6 — SDF volume creation and morphology.

The signed distance field is the workhorse of the pipeline. Default to
band/exterior density (not interior fill) so open CAD shells don't silently
under-fill.
"""


def mesh_to_sdf(obj, props):
    """Voxelize *obj* into a signed distance field grid. TODO(phase-4)."""
    raise NotImplementedError("volume.mesh_to_sdf is a Phase 4 stub")


def dilate(grid, distance):
    """Offset the isosurface outward by *distance* (clearance). TODO(phase-5)."""
    if not distance:
        return grid
    raise NotImplementedError("volume.dilate is a Phase 5 stub")


def morphological_close(grid, radius):
    """Dilate then erode by *radius*: fill small holes, drop thin internals.

    Replaces the plan's separate suppression + hole-fill phases. TODO(phase-6).
    """
    if not radius:
        return grid
    raise NotImplementedError("volume.morphological_close is a Phase 6 stub")
