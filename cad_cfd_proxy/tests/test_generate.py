# -*- coding: utf-8 -*-
"""Headless end-to-end test: full generate pipeline (Phases 2-4, 8).

Run with Blender::

    blender --background --factory-startup \
        --python cad_cfd_proxy/tests/test_generate.py

Runs ``core.generate_proxy`` on a cube in default EXTERNAL mode and asserts a
real proxy mesh comes out. Clearance/feature-size default to 0 (dilate/close
no-op) and smoothing/decimation are pending no-ops, so this exercises the
collect → cleanup → mesh_to_sdf → volume_to_mesh path. Exits non-zero on
failure for CI.
"""

import sys


def main():
    import bpy
    import cad_cfd_proxy
    from cad_cfd_proxy import core

    bpy.ops.wm.read_factory_settings(use_empty=True)
    cad_cfd_proxy.register()

    bpy.ops.mesh.primitive_cube_add(size=2.0, location=(0, 0, 0))
    cube = bpy.context.active_object
    cube.select_set(True)
    bpy.context.view_layer.objects.active = cube

    props = bpy.context.scene.cad_cfd_proxy
    props.source = "SELECTED"
    props.mode = "EXTERNAL"
    props.voxel_size = 0.1
    props.fill_volume = True

    proxy = core.generate_proxy(bpy.context, props)

    assert proxy is not None, "generate_proxy returned None"
    assert proxy.type == "MESH", "proxy is not a mesh"
    n_faces = len(proxy.data.polygons)
    assert n_faces > 0, "proxy mesh has no faces"
    # A voxelized 2 m cube at 0.1 m should be a substantial closed surface.
    assert len(proxy.data.vertices) > 8, "proxy suspiciously coarse"
    print("  proxy: %d verts, %d faces" % (len(proxy.data.vertices), n_faces))

    cad_cfd_proxy.unregister()
    print("generate test: PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        # Re-raise so Blender's --python-exit-code produces a non-zero exit.
        print("generate test: FAIL: %s" % exc)
        raise
