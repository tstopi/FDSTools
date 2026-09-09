# -*- coding: utf-8 -*-
"""Phases 4-6 — SDF volume creation and morphology.

The signed distance field is the workhorse of the pipeline.

Phase 4 (implemented here) joins the cleaned working meshes into one mesh and
runs a **Mesh to Volume** modifier on a Volume object, producing an OpenVDB
grid. The modifier computes an interior SDF band and (optionally) fills the
interior — filling is what lets Volume→Mesh yield a solid envelope rather than
a double-walled shell, at the cost of needing reasonably closed input (Phase 6
morphological close makes that robust on open CAD shells).

Phases 5 (clearance) and 6 (morphological close) use the API-only route:
Blender's public API exposes a density grid, not an SDF Python can offset, so
each morphological step is a **volume→mesh → offset-along-normals → mesh→volume**
round-trip. Blender's re-voxelization cleans up the self-intersections a raw
normal offset creates at concave corners, so this approximates true SDF
dilation/erosion. Close = dilate then erode. Each op no-ops at zero distance so
a default run passes straight through.
"""

import bmesh
import bpy

from . import collect

VOLUME_NAME = "CADCFD_volume"
JOINED_NAME = "CADCFD_joined"


class VolumeResult(object):
    """Handle for the generated volume and the geometry it came from.

    Carried between phases; ``volume_object`` holds the live Mesh-to-Volume
    modifier so later phases can tweak it or read the evaluated grids.
    """

    def __init__(self, volume_object, source_mesh, voxel_size):
        self.volume_object = volume_object
        self.source_mesh = source_mesh
        self.voxel_size = voxel_size

    def evaluated_grids(self, context):
        """Return the evaluated OpenVDB grids (list) for inspection/validation."""
        depsgraph = context.evaluated_depsgraph_get()
        vol_eval = self.volume_object.evaluated_get(depsgraph)
        return list(vol_eval.data.grids)


def mesh_to_sdf(context, sources, props):
    """Voxelize *sources* into an SDF volume. Returns a :class:`VolumeResult`."""
    work_coll = collect._ensure_work_collection(context)
    mesh_obj = _join_sources(sources, work_coll)
    return _voxelize_object(context, mesh_obj, props)


def dilate(context, grid, distance, props):
    """Offset the surface outward by *distance* (clearance) via a round-trip."""
    if not distance:
        return grid
    return _offset_and_revoxelize(context, grid, distance, props)


def morphological_close(context, grid, radius, props):
    """Dilate then erode by *radius*: fill holes/gaps and bridge parts smaller
    than ~2·radius, and suppress thin features. Replaces the plan's separate
    internal-suppression + hole-fill phases."""
    if not radius:
        return grid
    grid = _offset_and_revoxelize(context, grid, radius, props)    # dilate
    grid = _offset_and_revoxelize(context, grid, -radius, props)   # erode
    return grid


# --- helpers ---------------------------------------------------------------

def _voxelize_object(context, mesh_obj, props):
    """Build a Volume object + Mesh-to-Volume modifier from *mesh_obj*."""
    work_coll = collect._ensure_work_collection(context)
    volume = bpy.data.volumes.new(VOLUME_NAME)
    vol_obj = bpy.data.objects.new(VOLUME_NAME, volume)
    work_coll.objects.link(vol_obj)

    mod = vol_obj.modifiers.new("mesh_to_volume", "MESH_TO_VOLUME")
    mod.object = mesh_obj
    # Property names have drifted across Blender versions; set defensively so
    # the compat layer's "probe, don't assume" rule holds here too.
    _set_first(mod, ("resolution_mode",), "VOXEL_SIZE")
    _set_first(mod, ("voxel_size",), props.voxel_size)
    _set_first(mod, ("density",), 1.0)
    _set_first(mod, ("interior_band_width", "exterior_band_width"),
               props.interior_band_width)
    _set_first(mod, ("fill_volume", "use_fill_volume"), props.fill_volume)

    return VolumeResult(vol_obj, mesh_obj, props.voxel_size)


def _offset_and_revoxelize(context, grid, distance, props):
    """One morphological step: mesh the grid, offset every vertex along its
    normal by *distance* (negative = inward), then re-voxelize the result.

    Re-voxelization is what makes this robust — the raw normal offset self-
    intersects at concave corners, and rebuilding the SDF discards those.
    """
    # Local import avoids a circular import (surface imports volume).
    from . import surface
    mesh_obj = surface.volume_to_mesh(context, grid, props)
    _offset_along_normals(mesh_obj, distance)
    return _voxelize_object(context, mesh_obj, props)


def _offset_along_normals(obj, distance):
    """Move every vertex of *obj* along its outward normal by *distance*.

    Normals are recalculated outward first so +distance always dilates and
    -distance always erodes, regardless of the winding volume_to_mesh produced.
    """
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    bmesh.ops.recalc_face_normals(bm, faces=bm.faces)
    bm.normal_update()
    for v in bm.verts:
        v.co += v.normal * distance
    bm.to_mesh(obj.data)
    bm.free()
    obj.data.update()

def _join_sources(sources, work_coll):
    """Merge the working mesh objects into one mesh object.

    Phase 2 already baked world transforms into the mesh data, so joining the
    datablocks directly is correct.
    """
    bm = bmesh.new()
    for obj in sources:
        if obj.type == "MESH":
            bm.from_mesh(obj.data)
    mesh = bpy.data.meshes.new(JOINED_NAME)
    bm.to_mesh(mesh)
    bm.free()
    obj = bpy.data.objects.new(JOINED_NAME, mesh)
    work_coll.objects.link(obj)
    return obj


def _set_first(owner, names, value):
    """Set the first attribute in *names* that exists on *owner*.

    Returns the name used, or ``None`` if none matched (silently — a missing
    optional property is not fatal, and the diagnostics panel covers hard reqs).
    """
    for name in names:
        if hasattr(owner, name):
            setattr(owner, name, value)
            return name
    return None
