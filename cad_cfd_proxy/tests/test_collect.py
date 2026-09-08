# -*- coding: utf-8 -*-
"""Headless functional test for Phase 2 (input processing).

Run with Blender::

    blender --background --factory-startup \
        --python cad_cfd_proxy/tests/test_collect.py

Builds a small scene exercising the three things Phase 2 must handle —
modifier application, curve→mesh conversion, and collection-instance
realization — then asserts that ``collect.gather`` produces real geometry and
leaves the originals untouched. Exits non-zero on failure for CI.
"""

import sys


def _clear_scene():
    import bpy
    bpy.ops.wm.read_factory_settings(use_empty=True)


def _build_scene():
    import bpy

    # 1) A mesh with a modifier (must be applied on realize).
    bpy.ops.mesh.primitive_cube_add(location=(0, 0, 0))
    cube = bpy.context.active_object
    cube.name = "Cube"
    cube.modifiers.new("bev", "BEVEL")

    # 2) A curve (must convert to mesh).
    bpy.ops.curve.primitive_bezier_circle_add(location=(4, 0, 0))
    curve = bpy.context.active_object
    curve.name = "Ring"
    curve.data.extrude = 0.1  # give it surface area so it meshes to polys

    # 3) A collection instance (must realize).
    src_coll = bpy.data.collections.new("InstSrc")
    bpy.context.scene.collection.children.link(src_coll)
    bpy.ops.mesh.primitive_uv_sphere_add(location=(0, 0, 0))
    sphere = bpy.context.active_object
    sphere.name = "Sphere"
    # move sphere out of the scene master collection into src_coll
    for c in list(sphere.users_collection):
        c.objects.unlink(sphere)
    src_coll.objects.link(sphere)
    inst = bpy.data.objects.new("SphereInstance", None)
    inst.instance_type = "COLLECTION"
    inst.instance_collection = src_coll
    inst.location = (-4, 0, 0)
    bpy.context.scene.collection.objects.link(inst)

    return [cube, curve, inst]


def main():
    import bpy
    import cad_cfd_proxy
    from cad_cfd_proxy.core import collect

    _clear_scene()
    cad_cfd_proxy.register()

    sources = _build_scene()
    n_objects_before = len(bpy.data.objects)

    # Select the top-level sources and gather.
    for obj in bpy.data.objects:
        obj.select_set(False)
    for obj in sources:
        obj.select_set(True)
    bpy.context.view_layer.objects.active = sources[0]

    props = bpy.context.scene.cad_cfd_proxy
    props.source = "SELECTED"
    props.patch_naming = "PER_OBJECT"

    produced = collect.gather(bpy.context, props)

    # Working geometry was produced with real faces.
    assert produced, "gather produced nothing"
    for obj in produced:
        assert obj.type == "MESH" and obj.data.polygons, \
            "working object %r has no faces" % obj.name
        assert collect.SOURCE_KEY_PROP in obj, "missing source key on %r" % obj.name
    print("  produced %d working object(s): %s"
          % (len(produced), ", ".join(o.name for o in produced)))

    # Originals preserved: the cube still has its unapplied modifier, the curve
    # is still a curve, and the sphere source still lives in its collection.
    cube = bpy.data.objects["Cube"]
    assert cube.modifiers, "original cube modifier was consumed"
    assert bpy.data.objects["Ring"].type == "CURVE", "original curve was converted"
    assert bpy.data.objects["Sphere"].type == "MESH", "instance source disturbed"

    # Only the new work objects (+ their collection) were added.
    added = len(bpy.data.objects) - n_objects_before
    assert added == len(produced), \
        "expected %d new objects, got %d" % (len(produced), added)

    cad_cfd_proxy.unregister()
    print("collect test: PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print("collect test: FAIL: %s" % exc)
        sys.exit(1)
