# -*- coding: utf-8 -*-
"""Phase 2 — input processing.

Gather source objects, realize instances, apply modifiers, and convert
supported object types to mesh — **without touching the originals**. The result
is a set of freshly created working mesh objects living in a dedicated work
collection, ready for cleanup (Phase 3) and voxelization (Phase 4).

The heavy lifting is done through the dependency graph: iterating
``depsgraph.object_instances`` and calling ``bpy.data.meshes.new_from_object``
realizes collection/vertex/face instances, bakes modifier stacks, and converts
curves, text, surfaces and metaballs to mesh — all in one pass. We never modify
the source objects, so originals are preserved by construction.
"""

import bpy
import bmesh

from .errors import CollectError

WORK_COLLECTION = "CADCFD_work"

# Object types ``new_from_object`` can evaluate to a mesh.
_CONVERTIBLE = {"MESH", "CURVE", "SURFACE", "FONT", "META"}

# Custom property recording which source an object/patch came from.
SOURCE_KEY_PROP = "cadcfd_source"


def gather(context, props):
    """Return a list of working mesh objects realized from the source set.

    Raises
    ------
    CollectError
        If no valid source geometry can be resolved.
    """
    sources = set(_resolve_sources(context, props))
    if not sources:
        raise CollectError("No source objects found for the selected source mode")

    depsgraph = context.evaluated_depsgraph_get()
    work_coll = _ensure_work_collection(context)

    # Group realized meshes by patch key so export (Phase 11) can keep identity.
    groups = {}
    for inst in depsgraph.object_instances:
        eval_obj, key = _classify_instance(inst, sources, props)
        if eval_obj is None:
            continue
        if eval_obj.type not in _CONVERTIBLE:
            continue
        mesh = _mesh_from_eval(eval_obj, depsgraph, inst.matrix_world)
        if mesh is None:
            continue
        groups.setdefault(key, []).append(mesh)

    if not groups:
        raise CollectError("Source objects produced no convertible geometry")

    produced = []
    for key, meshes in groups.items():
        obj = _join_meshes(key, meshes, work_coll)
        obj[SOURCE_KEY_PROP] = key
        produced.append(obj)
    return produced


# --- source resolution -----------------------------------------------------

def _resolve_sources(context, props):
    """Resolve the set of source objects per the chosen source mode.

    ``Collection.all_objects`` is recursive, so nested collections are covered.
    """
    mode = props.source
    if mode == "SELECTED":
        return list(context.selected_objects)
    if mode == "ACTIVE_COLLECTION":
        layer = context.view_layer.active_layer_collection
        if layer is None:
            raise CollectError("No active collection")
        return list(layer.collection.all_objects)
    if mode == "NAMED_COLLECTION":
        coll = bpy.data.collections.get(props.source_collection)
        if coll is None:
            raise CollectError("Collection %r not found" % props.source_collection)
        return list(coll.all_objects)
    raise CollectError("Unknown source mode %r" % mode)


def _classify_instance(inst, sources, props):
    """Return ``(eval_obj, patch_key)`` for an instance we should keep, else
    ``(None, None)``.

    Real (non-instance) objects are kept when they're in the source set;
    generated instances are kept when their instancer is in the source set.
    """
    if inst.is_instance:
        instancer = inst.parent.original if inst.parent else None
        if instancer not in sources:
            return None, None
        return inst.object, _patch_key(instancer, props)
    orig = inst.object.original
    if orig not in sources:
        return None, None
    return inst.object, _patch_key(orig, props)


def _patch_key(obj, props):
    """Patch/region name for *obj* under the chosen naming scheme (Phase 11)."""
    scheme = props.patch_naming
    if scheme == "SINGLE":
        return "proxy"
    if scheme == "PER_COLLECTION":
        colls = obj.users_collection
        return colls[0].name if colls else "proxy"
    return obj.name  # PER_OBJECT


# --- mesh realization -------------------------------------------------------

def _mesh_from_eval(eval_obj, depsgraph, matrix_world):
    """Evaluated object → new world-space mesh datablock, or ``None`` if empty.

    Applies modifiers and converts to mesh via ``new_from_object``; the caller
    owns the returned datablock and must free it (``_join_meshes`` does).
    """
    try:
        mesh = bpy.data.meshes.new_from_object(
            eval_obj, preserve_all_data_layers=False, depsgraph=depsgraph)
    except RuntimeError:
        return None
    if mesh is None:
        return None
    if not mesh.polygons:
        bpy.data.meshes.remove(mesh)
        return None
    # Instances carry their transform in matrix_world; bake it into the copy.
    mesh.transform(matrix_world)
    return mesh


def _join_meshes(key, meshes, work_coll):
    """Merge realized *meshes* into one working object linked to *work_coll*.

    Consumes (frees) the input mesh datablocks.
    """
    bm = bmesh.new()
    for mesh in meshes:
        bm.from_mesh(mesh)      # appends; geometry is already world-space
        bpy.data.meshes.remove(mesh)
    name = "%s_proxy_src" % key
    joined = bpy.data.meshes.new(name)
    bm.to_mesh(joined)
    bm.free()
    obj = bpy.data.objects.new(name, joined)
    work_coll.objects.link(obj)
    return obj


# --- work collection --------------------------------------------------------

def _ensure_work_collection(context):
    """Get or create the scene-linked collection holding working geometry."""
    coll = bpy.data.collections.get(WORK_COLLECTION)
    if coll is None:
        coll = bpy.data.collections.new(WORK_COLLECTION)
    if coll.name not in context.scene.collection.children:
        context.scene.collection.children.link(coll)
    return coll
