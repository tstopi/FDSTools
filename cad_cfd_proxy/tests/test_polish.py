# -*- coding: utf-8 -*-
"""Headless functional test for Phases 9, 14, 15 (smooth / estimate / cleanup).

Run with Blender::

    blender --background --factory-startup \
        --python cad_cfd_proxy/tests/test_polish.py

Exits non-zero on failure for CI.
"""

import sys


def _test_smooth():
    import bpy
    from cad_cfd_proxy.core import surface

    bpy.ops.mesh.primitive_cube_add(size=2.0)
    obj = bpy.context.active_object
    n_faces = len(obj.data.polygons)
    before = [v.co.copy() for v in obj.data.vertices]

    props = bpy.context.scene.cad_cfd_proxy
    props.smooth_iterations = 3
    props.smooth_strength = 0.5
    surface.smooth(obj, props)

    assert len(obj.data.polygons) == n_faces, "smoothing changed topology"
    moved = any((v.co - b).length > 1e-6 for v, b in zip(obj.data.vertices, before))
    assert moved, "smoothing moved no vertices"
    print("  smooth: topology preserved, vertices moved")


def _select_only(obj):
    import bpy
    for o in bpy.data.objects:
        o.select_set(o is obj)
    bpy.context.view_layer.objects.active = obj


def _test_estimate():
    import bpy
    from cad_cfd_proxy.core import estimate

    bpy.ops.mesh.primitive_cube_add(size=2.0)
    cube = bpy.context.active_object
    _select_only(cube)

    props = bpy.context.scene.cad_cfd_proxy
    props.source = "SELECTED"
    props.voxel_size = 0.1

    est = estimate.estimate(bpy.context, props)
    print("  estimate: %s" % "; ".join(estimate.format_lines(est)))
    assert est["surface_area"] > 20.0, "cube area ~24 expected, got %.1f" % est["surface_area"]
    assert est["sparse_voxels"] > 0 and est["est_faces"] > 0, "empty estimate"
    assert 1.9 < est["bbox"][0] < 2.1, "bbox x ~2 expected, got %.2f" % est["bbox"][0]


def _test_finalize():
    import bpy
    from cad_cfd_proxy import core
    from cad_cfd_proxy.core import collect

    # Default run: intermediates cleaned, proxy re-homed in the scene.
    bpy.ops.mesh.primitive_cube_add(size=2.0)
    cube = bpy.context.active_object
    _select_only(cube)
    props = bpy.context.scene.cad_cfd_proxy
    props.source = "SELECTED"
    props.mode = "EXTERNAL"
    props.voxel_size = 0.15
    props.keep_intermediates = False

    proxy = core.generate_proxy(bpy.context, props)
    work = bpy.data.collections.get(collect.WORK_COLLECTION)
    assert work is None or len(work.objects) == 0, "intermediates were not cleaned"
    assert proxy.name in bpy.context.scene.collection.objects, "proxy not re-homed"
    print("  finalize: intermediates cleaned, proxy in scene")

    # keep_intermediates=True leaves the work collection populated (no reset:
    # re-registering would fail; a fresh selection is enough).
    bpy.ops.mesh.primitive_cube_add(size=2.0, location=(10, 0, 0))
    cube2 = bpy.context.active_object
    _select_only(cube2)
    props.voxel_size = 0.15
    props.keep_intermediates = True
    core.generate_proxy(bpy.context, props)
    work = bpy.data.collections.get(collect.WORK_COLLECTION)
    assert work is not None and len(work.objects) > 0, "intermediates were dropped"
    print("  finalize: keep_intermediates retained work objects")


def main():
    import bpy
    import cad_cfd_proxy

    bpy.ops.wm.read_factory_settings(use_empty=True)
    cad_cfd_proxy.register()

    _test_smooth()
    _test_estimate()
    _test_finalize()

    print("polish test: PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print("polish test: FAIL: %s" % exc)
        raise
