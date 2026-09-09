# -*- coding: utf-8 -*-
"""Headless functional test for Phase 12 (validation).

Run with Blender::

    blender --background --factory-startup \
        --python cad_cfd_proxy/tests/test_validate.py

Checks the validator on three meshes: a clean cube (watertight/manifold/
FDS-ready), an open cube (not watertight), and two crossing quads (a genuine
self-intersection). Exits non-zero on failure for CI.
"""

import sys


def _cube_object(name="c"):
    import bpy
    bpy.ops.mesh.primitive_cube_add(size=2.0)
    obj = bpy.context.active_object
    obj.name = name
    return obj


def main():
    import bpy
    import bmesh
    import cad_cfd_proxy
    from cad_cfd_proxy.core import validate

    bpy.ops.wm.read_factory_settings(use_empty=True)
    cad_cfd_proxy.register()

    # 1) Clean cube: watertight, manifold, oriented, no self-intersections.
    cube = _cube_object("clean")
    rep = validate.check(cube)
    print("  clean cube: %s" % rep.summary())
    assert rep.watertight and rep.manifold and rep.normals_consistent, \
        "clean cube should be watertight/manifold/oriented"
    assert rep.self_intersections == 0, "clean cube has no self-intersections"
    assert rep.fds_ready(), "clean cube should be FDS-ready"

    # 2) Open cube: remove one face → not watertight, not FDS-ready.
    holed = _cube_object("holed")
    bm = bmesh.new()
    bm.from_mesh(holed.data)
    bm.faces.ensure_lookup_table()
    bmesh.ops.delete(bm, geom=[bm.faces[0]], context="FACES")
    bm.to_mesh(holed.data)
    bm.free()
    rep = validate.check(holed)
    print("  holed cube: %s" % rep.summary())
    assert not rep.watertight, "holed cube should not be watertight"
    assert not rep.fds_ready(), "holed cube should not be FDS-ready"

    # 3) Two crossing quads in one mesh: a real self-intersection.
    bm = bmesh.new()
    z = [bm.verts.new(v) for v in
         [(-1, -1, 0), (1, -1, 0), (1, 1, 0), (-1, 1, 0)]]       # quad in z=0
    x = [bm.verts.new(v) for v in
         [(0, -1, -1), (0, 1, -1), (0, 1, 1), (0, -1, 1)]]        # quad in x=0
    bm.faces.new(z)
    bm.faces.new(x)
    cross_mesh = bpy.data.meshes.new("cross")
    bm.to_mesh(cross_mesh)
    bm.free()
    cross = bpy.data.objects.new("cross", cross_mesh)
    bpy.context.scene.collection.objects.link(cross)
    rep = validate.check(cross)
    print("  crossing quads: %s" % rep.summary())
    assert rep.self_intersections >= 1, "crossing quads should self-intersect"

    cad_cfd_proxy.unregister()
    print("validate test: PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print("validate test: FAIL: %s" % exc)
        raise
