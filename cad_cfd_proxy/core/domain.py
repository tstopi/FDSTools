# -*- coding: utf-8 -*-
"""Phase 7 — internal fluid volume.

Produces the **flow domain**: a bounding box (or a user-supplied domain object)
with the solid subtracted — the fluid region *around* the object, as a
wind-tunnel / external-aerodynamics domain wants it and as snappyHexMesh meshes
it. The obstacle is first simplified via the SDF pipeline (watertight/manifold),
which also makes the boolean robust; then a boolean DIFFERENCE carves it out of
the domain box.

Returns a mesh object directly (not a grid): re-voxelizing the carved region
with fill enabled would refill the obstacle-shaped void and destroy it, so the
boolean result is the fluid proxy.

Scope note: this is the flow-around-object domain. Extracting the enclosed
interior of a hollow part (duct/manifold internals) is a different operation
(no fill + flood-fill from a seed) and is out of scope for v1.
"""

import bmesh
import bpy
import mathutils

from . import collect, volume

DOMAIN_NAME = "CADCFD_domain"


def extract_fluid_volume(context, sources, solid_grid, props):
    """Return the fluid-domain proxy mesh object = domain − simplified solid."""
    from . import surface  # local import avoids a circular import

    work_coll = collect._ensure_work_collection(context)
    # Simplify the obstacle first: watertight/manifold makes the boolean robust.
    obstacle = surface.volume_to_mesh(context, solid_grid, props)
    domain_obj = _domain_object(context, sources, props, work_coll)

    mod = domain_obj.modifiers.new("fluid_bool", "BOOLEAN")
    mod.operation = "DIFFERENCE"
    mod.object = obstacle
    volume._set_first(mod, ("solver",), "EXACT")  # robust to imperfect input
    surface._bake_modifiers(context, domain_obj)

    domain_obj.name = surface.PROXY_NAME
    return domain_obj


# --- domain construction ---------------------------------------------------

def _domain_object(context, sources, props, work_coll):
    """The flow-domain solid: a user domain object if set, else a padded bbox."""
    if props.domain_object is not None:
        return _world_mesh_copy(props.domain_object, work_coll)
    lo, hi = _world_bounds(sources)
    pad = props.domain_padding
    lo = [c - pad for c in lo]
    hi = [c + pad for c in hi]
    return _box_object(lo, hi, work_coll)


def _world_bounds(sources):
    """Axis-aligned world-space bounds over the source meshes."""
    inf = float("inf")
    lo = [inf, inf, inf]
    hi = [-inf, -inf, -inf]
    for obj in sources:
        if obj.type != "MESH":
            continue
        mw = obj.matrix_world
        for v in obj.data.vertices:
            co = mw @ v.co
            for i in range(3):
                lo[i] = min(lo[i], co[i])
                hi[i] = max(hi[i], co[i])
    return lo, hi


def _box_object(lo, hi, work_coll):
    """A closed box mesh object spanning lo..hi."""
    center = mathutils.Vector(((lo[0] + hi[0]) * 0.5,
                               (lo[1] + hi[1]) * 0.5,
                               (lo[2] + hi[2]) * 0.5))
    size = (hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2])
    bm = bmesh.new()
    bmesh.ops.create_cube(bm, size=1.0)  # unit cube centered at origin
    bm.transform(mathutils.Matrix.Translation(center)
                 @ mathutils.Matrix.Diagonal((size[0], size[1], size[2], 1.0)))
    mesh = bpy.data.meshes.new(DOMAIN_NAME)
    bm.to_mesh(mesh)
    bm.free()
    obj = bpy.data.objects.new(DOMAIN_NAME, mesh)
    work_coll.objects.link(obj)
    return obj


def _world_mesh_copy(src_obj, work_coll):
    """A world-space mesh-object copy of a user-supplied domain object."""
    mesh = src_obj.data.copy()
    mesh.transform(src_obj.matrix_world)
    obj = bpy.data.objects.new(DOMAIN_NAME, mesh)
    work_coll.objects.link(obj)
    return obj
