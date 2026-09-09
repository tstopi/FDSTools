# -*- coding: utf-8 -*-
"""Headless functional test for Phase 11 (export).

Run with Blender::

    blender --background --factory-startup \
        --python cad_cfd_proxy/tests/test_export.py

Exports a cube to snappyHexMesh STL and FDS &GEOM and checks the files have the
expected structure and metre scaling. Exits non-zero on failure for CI.
"""

import os
import sys
import tempfile


def main():
    import bpy
    import cad_cfd_proxy
    from cad_cfd_proxy import core

    bpy.ops.wm.read_factory_settings(use_empty=True)
    cad_cfd_proxy.register()

    # A 1000 mm cube; export scale 0.001 should yield a 1 m cube.
    bpy.ops.mesh.primitive_cube_add(size=1000.0)
    proxy = bpy.context.active_object
    props = bpy.context.scene.cad_cfd_proxy
    props.export_scale = 0.001

    tmp = tempfile.mkdtemp(prefix="cadcfd_")

    # --- STL ---
    stl_path = os.path.join(tmp, "proxy.stl")
    props.export_target = "SNAPPY"
    props.export_path = stl_path
    props.patch_naming = "SINGLE"
    core.export_proxy(bpy.context, proxy, props)
    stl = open(stl_path).read()
    assert stl.startswith("solid proxy"), "STL missing solid header"
    assert "facet normal" in stl and "endsolid" in stl, "STL malformed"
    # 1000 mm * 0.001 = 0.5 m half-extent → coordinates ~0.5, not ~500.
    assert "5.000000e-01" in stl, "STL not scaled to metres"
    print("  STL: %d bytes" % len(stl))

    # --- FDS &GEOM ---
    fds_path = os.path.join(tmp, "proxy.geom")
    props.export_target = "FDS"
    props.export_path = fds_path
    core.export_proxy(bpy.context, proxy, props)
    fds = open(fds_path).read()
    assert "&GEOM" in fds and "VERTS=" in fds and "FACES=" in fds, "GEOM malformed"
    assert fds.rstrip().endswith("/"), "GEOM namelist not terminated"
    print("  GEOM: %d bytes" % len(fds))

    cad_cfd_proxy.unregister()
    print("export test: PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print("export test: FAIL: %s" % exc)
        raise
