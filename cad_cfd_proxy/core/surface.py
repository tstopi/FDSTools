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
    """Laplacian-style vertex smoothing, in place (Phase 9).

    Gentle by default. Smoothing can introduce self-intersections, which Phase
    12 validation reports. No-op at zero iterations.
    """
    import bmesh
    iterations = props.smooth_iterations
    if iterations <= 0:
        return obj
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    for _ in range(iterations):
        bmesh.ops.smooth_vert(
            bm, verts=bm.verts, factor=props.smooth_strength,
            use_axis_x=True, use_axis_y=True, use_axis_z=True)
    bm.to_mesh(obj.data)
    bm.free()
    obj.data.update()
    return obj


def decimate_to_target(context, obj, target_faces, tolerance=0.05, max_iter=12):
    """Decimate *obj* to within *tolerance* of *target_faces* (triangles).

    Collapse decimation triangulates and its output count is ~monotonic in the
    ratio, so we seed the ratio with the linear estimate then bisect. Does
    nothing if the target is at or above the current triangle count (we can't
    add detail). Bakes the modifier when done.
    """
    current = _tri_count(obj)
    if target_faces <= 0 or current <= target_faces * (1.0 + tolerance):
        return obj

    mod = obj.modifiers.new("decimate", "DECIMATE")
    mod.decimate_type = "COLLAPSE"

    lo, hi = 0.0, 1.0
    ratio = min(1.0, float(target_faces) / current)  # linear seed
    best_ratio = ratio
    best_err = None
    for _ in range(max_iter):
        mod.ratio = ratio
        faces = _evaluated_face_count(context, obj)
        err = abs(faces - target_faces)
        if best_err is None or err < best_err:
            best_err, best_ratio = err, ratio
        if err <= tolerance * target_faces:
            break
        if faces > target_faces:
            hi = ratio
        else:
            lo = ratio
        ratio = 0.5 * (lo + hi)

    # Bake the best ratio seen, not whatever the loop last computed-but-never-
    # measured (the trailing midpoint would otherwise be discarded unmeasured).
    mod.ratio = best_ratio
    _bake_modifiers(context, obj)
    return obj


# --- helpers ---------------------------------------------------------------

def _tri_count(obj):
    """Triangle count of a mesh with mixed polygon sizes (fan triangulation)."""
    return sum(len(p.vertices) - 2 for p in obj.data.polygons)


def _evaluated_face_count(context, obj):
    """Polygon count of *obj* after its modifier stack is evaluated."""
    depsgraph = context.evaluated_depsgraph_get()
    return len(obj.evaluated_get(depsgraph).data.polygons)


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
