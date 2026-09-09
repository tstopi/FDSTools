# -*- coding: utf-8 -*-
"""Headless functional test for Phase 10 (adaptive decimation).

Run with Blender::

    blender --background --factory-startup \
        --python cad_cfd_proxy/tests/test_decimate.py

Decimates a dense icosphere to a target triangle count and asserts the result
lands within tolerance. Exits non-zero on failure for CI.
"""

import sys


def main():
    import bpy
    import cad_cfd_proxy
    from cad_cfd_proxy.core import surface

    bpy.ops.wm.read_factory_settings(use_empty=True)
    cad_cfd_proxy.register()

    bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=6)  # ~20k tris
    obj = bpy.context.active_object
    start = len(obj.data.polygons)

    target = 2000
    surface.decimate_to_target(bpy.context, obj, target, tolerance=0.05)
    result = len(obj.data.polygons)
    print("  decimate: %d -> %d (target %d)" % (start, result, target))

    # Within 10% (function aims for 5%; allow slack for collapse non-linearity).
    assert abs(result - target) <= 0.10 * target, \
        "decimation missed target: %d vs %d" % (result, target)
    # Modifier was baked, not left live.
    assert not obj.modifiers, "decimate modifier was not baked"

    cad_cfd_proxy.unregister()
    print("decimate test: PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print("decimate test: FAIL: %s" % exc)
        raise
