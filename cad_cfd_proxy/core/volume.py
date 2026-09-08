# -*- coding: utf-8 -*-
"""Phases 4-6 — SDF volume creation and morphology.

The signed distance field is the workhorse of the pipeline.

Phase 4 (implemented here) joins the cleaned working meshes into one mesh and
runs a **Mesh to Volume** modifier on a Volume object, producing an OpenVDB
grid. The modifier computes an interior SDF band and (optionally) fills the
interior — filling is what lets Volume→Mesh yield a solid envelope rather than
a double-walled shell, at the cost of needing reasonably closed input (Phase 6
morphological close makes that robust on open CAD shells).

Phases 5 (dilate/clearance) and 6 (morphological close) operate on the same
grid and are stubbed until those phases land. True SDF offset/morphology needs
grid-level access (OpenVDB); the stubs no-op when their distance is zero so a
default run passes straight through.
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


# --- helpers ---------------------------------------------------------------

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
