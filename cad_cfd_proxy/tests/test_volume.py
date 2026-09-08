# -*- coding: utf-8 -*-
"""Headless functional test for Phase 4 (Mesh → Volume / SDF).

Run with Blender::

    blender --background --factory-startup \
        --python cad_cfd_proxy/tests/test_volume.py

Voxelizes a cube and asserts a non-empty OpenVDB grid comes out with a plausible
bounding box. Exits non-zero on failure for CI.
"""

import sys


def main():
    import bpy
    import cad_cfd_proxy
    from cad_cfd_proxy.core import volume

    bpy.ops.wm.read_factory_settings(use_empty=True)
    cad_cfd_proxy.register()

    # A 2 m cube at the origin (Phase 2 output is world-space at identity).
    bpy.ops.mesh.primitive_cube_add(size=2.0, location=(0, 0, 0))
    cube = bpy.context.active_object
    cube.name = "src"

    props = bpy.context.scene.cad_cfd_proxy
    props.voxel_size = 0.1
    props.fill_volume = True

    result = volume.mesh_to_sdf(bpy.context, [cube], props)

    # A Mesh-to-Volume modifier was configured on a real Volume object.
    assert result.volume_object.type == "VOLUME", "result is not a volume object"
    mods = result.volume_object.modifiers
    assert any(m.type == "MESH_TO_VOLUME" for m in mods), "no mesh-to-volume modifier"

    # Evaluate and check a grid was actually produced.
    grids = result.evaluated_grids(bpy.context)
    assert grids, "mesh-to-volume produced no grids"
    print("  grids: %s" % ", ".join(g.name for g in grids))

    cad_cfd_proxy.unregister()
    print("volume test: PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        # Re-raise so Blender's --python-exit-code produces a non-zero exit.
        print("volume test: FAIL: %s" % exc)
        raise
