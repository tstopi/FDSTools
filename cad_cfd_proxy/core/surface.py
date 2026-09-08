# -*- coding: utf-8 -*-
"""Phases 8-10 — surface reconstruction, smoothing, decimation.

Volume → Mesh yields a watertight, manifold isosurface by construction; the
real risk downstream is *self-intersection* from smoothing/decimation, checked
in :mod:`.validate`.
"""


def volume_to_mesh(grid, props):
    """Adaptive isosurface extraction to a mesh object. TODO(phase-8)."""
    raise NotImplementedError("surface.volume_to_mesh is a Phase 8 stub")


def smooth(obj, props):
    """Smooth position + optional Laplacian smooth, in place. TODO(phase-9)."""
    raise NotImplementedError("surface.smooth is a Phase 9 stub")


def decimate_to_target(obj, target_faces, tolerance=0.05):
    """Decimate to within *tolerance* of *target_faces* by bisecting the
    collapse ratio (monotonic → converges in a few passes). TODO(phase-10)."""
    raise NotImplementedError("surface.decimate_to_target is a Phase 10 stub")
