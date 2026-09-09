# -*- coding: utf-8 -*-
"""Headless functional test for Phase 3 (geometry cleanup).

Run with Blender::

    blender --background --factory-startup \
        --python cad_cfd_proxy/tests/test_cleanup.py

Builds a mesh with a coincident (double) vertex and a stray loose vertex, then
asserts cleanup welds the double, drops the loose vert, and leaves every
remaining vertex attached to a face. Exits non-zero on failure for CI.
"""

import sys


def main():
    import bpy
    import bmesh
    import cad_cfd_proxy
    from cad_cfd_proxy.core import cleanup

    bpy.ops.wm.read_factory_settings(use_empty=True)
    cad_cfd_proxy.register()

    # One quad + a coincident double on a corner + a far-away loose vertex.
    bm = bmesh.new()
    v1 = bm.verts.new((0.0, 0.0, 0.0))
    v2 = bm.verts.new((1.0, 0.0, 0.0))
    v3 = bm.verts.new((1.0, 1.0, 0.0))
    v4 = bm.verts.new((0.0, 1.0, 0.0))
    bm.faces.new((v1, v2, v3, v4))
    bm.verts.new((0.0, 0.0, 0.0))   # double, coincident with v1 (loose)
    bm.verts.new((5.0, 5.0, 5.0))   # clearly loose
    mesh = bpy.data.meshes.new("dirty")
    bm.to_mesh(mesh)
    bm.free()
    obj = bpy.data.objects.new("dirty", mesh)
    bpy.context.scene.collection.objects.link(obj)

    props = bpy.context.scene.cad_cfd_proxy
    props.merge_distance = 1e-4
    props.delete_loose = True
    props.fix_normals = True

    cleanup.clean([obj], props)

    assert len(mesh.vertices) == 4, \
        "expected 4 verts after cleanup, got %d" % len(mesh.vertices)

    check = bmesh.new()
    check.from_mesh(mesh)
    orphans = [v for v in check.verts if not v.link_faces]
    check.free()
    assert not orphans, "%d loose vertices survived cleanup" % len(orphans)

    cad_cfd_proxy.unregister()
    print("cleanup test: PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        # Re-raise so Blender's --python-exit-code produces a non-zero exit.
        print("cleanup test: FAIL: %s" % exc)
        raise
