# -*- coding: utf-8 -*-
"""Headless test for the progress generator and cancel cleanup.

Run with Blender::

    blender --background --factory-startup \
        --python cad_cfd_proxy/tests/test_progress.py

The modal operator and panel progress bar need a window, so they aren't tested
here; this covers the generate_job generator (monotonic progress, proxy result)
and discard_work (removes partial intermediates on cancel). Exits non-zero on
failure for CI.
"""

import sys


def main():
    import bpy
    import cad_cfd_proxy
    from cad_cfd_proxy import core
    from cad_cfd_proxy.core import collect, volume

    bpy.ops.wm.read_factory_settings(use_empty=True)
    cad_cfd_proxy.register()

    bpy.ops.mesh.primitive_cube_add(size=2.0)
    cube = bpy.context.active_object
    for o in bpy.data.objects:
        o.select_set(o is cube)
    bpy.context.view_layer.objects.active = cube

    props = bpy.context.scene.cad_cfd_proxy
    props.source = "SELECTED"
    props.mode = "EXTERNAL"
    props.voxel_size = 0.15

    # Drive the generator like the modal operator does.
    gen = core.generate_job(bpy.context, props)
    fractions = []
    proxy = None
    try:
        while True:
            label, frac = next(gen)
            fractions.append(frac)
    except StopIteration as stop:
        proxy = stop.value

    print("  progress steps: %s" % ", ".join("%.2f" % f for f in fractions))
    assert proxy is not None and len(proxy.data.polygons) > 0, "no proxy produced"
    assert fractions == sorted(fractions), "progress not monotonic"
    assert abs(fractions[-1] - 1.0) < 1e-9, "progress did not reach 1.0"
    assert 0.0 < fractions[0] < 1.0, "first step out of range"

    # discard_work removes a partially-built work collection (cancel path).
    grid = volume.mesh_to_sdf(bpy.context, [cube], props)
    assert bpy.data.collections.get(collect.WORK_COLLECTION) is not None, \
        "mesh_to_sdf did not create the work collection"
    core.discard_work(bpy.context)
    assert bpy.data.collections.get(collect.WORK_COLLECTION) is None, \
        "discard_work left the work collection behind"
    print("  discard_work: cleaned partial intermediates")

    cad_cfd_proxy.unregister()
    print("progress test: PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print("progress test: FAIL: %s" % exc)
        raise
