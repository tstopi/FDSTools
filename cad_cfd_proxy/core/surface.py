# -*- coding: utf-8 -*-
"""Phases 8-10 — surface reconstruction, smoothing, decimation.

Phase 8 (implemented) runs a Volume→Mesh modifier on the SDF grid and bakes the
result to a real mesh. The output is a watertight, manifold isosurface by
construction; the downstream risk is *self-intersection* introduced by
smoothing/decimation, which Phase 12 will check.

Phases 9 (smoothing) and 10 (decimation) are graceful no-ops for now — they
return the mesh unchanged so the pipeline runs end-to-end — and get real
implementations when those phases land.
"""

import bpy

from . import collect, volume

PROXY_NAME = "CADCFD_proxy"


def volume_to_mesh(context, grid, props):
    """Extract an adaptive isosurface mesh from *grid*. Returns the proxy object."""
    work_coll = collect._ensure_work_collection(context)
    mesh = bpy.data.meshes.new(PROXY_NAME)
    obj = bpy.data.objects.new(PROXY_NAME, mesh)
    work_coll.objects.link(obj)

    mod = obj.modifiers.new("volume_to_mesh", "VOLUME_TO_MESH")
    mod.object = grid.volume_object
    volume._set_first(mod, ("grid_name",), "density")
    volume._set_first(mod, ("threshold",), props.surface_threshold)
    volume._set_first(mod, ("adaptivity",), props.adaptivity)
    # Use the grid's own resolution rather than re-sampling.
    volume._set_first(mod, ("resolution_mode",), "GRID")

    _bake_modifiers(context, obj)
    return obj


def smooth(obj, props):
    """Smooth position + optional Laplacian smooth, in place. TODO(phase-9)."""
    # Pending: no-op so the pipeline runs end-to-end.
    return obj


def decimate_to_target(obj, target_faces, tolerance=0.05):
    """Decimate to within *tolerance* of *target_faces* by bisecting the
    collapse ratio (monotonic → converges in a few passes). TODO(phase-10)."""
    # Pending: no-op so the pipeline runs end-to-end.
    return obj


# --- helpers ---------------------------------------------------------------

def _bake_modifiers(context, obj):
    """Apply *obj*'s modifier stack by evaluating the depsgraph, leaving *obj*
    with real, modifier-free geometry (and no cross-object dependencies)."""
    depsgraph = context.evaluated_depsgraph_get()
    eval_obj = obj.evaluated_get(depsgraph)
    baked = bpy.data.meshes.new_from_object(eval_obj, depsgraph=depsgraph)
    obj.modifiers.clear()
    old = obj.data
    obj.data = baked
    if old.users == 0:
        bpy.data.meshes.remove(old)
