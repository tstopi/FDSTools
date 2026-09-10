# -*- coding: utf-8 -*-
"""Headless tests for the code-review fixes.

Run with Blender::

    blender --background --factory-startup \
        --python cad_cfd_proxy/tests/test_review.py

Covers: preset application (the enum now writes parameters) and the export
error contract (write raises PipelineError on no path / non-solver-ready FDS
geometry, which the Export operator turns into a clean report). Exits non-zero
on failure for CI.
"""

import os
import sys
import tempfile


def _test_preset():
    import bpy
    from cad_cfd_proxy.properties import PRESET_VALUES

    props = bpy.context.scene.cad_cfd_proxy
    props.preset = "CUSTOM"
    props.voxel_size = 0.123  # sentinel

    props.preset = "VEHICLE"
    vs, cl, ft = PRESET_VALUES["VEHICLE"]
    assert abs(props.voxel_size - vs) < 1e-9, "preset did not set voxel_size"
    assert abs(props.clearance - cl) < 1e-9, "preset did not set clearance"
    assert abs(props.feature_size - ft) < 1e-9, "preset did not set feature_size"

    # CUSTOM must not overwrite the current values.
    props.preset = "CUSTOM"
    assert abs(props.voxel_size - vs) < 1e-9, "CUSTOM preset clobbered voxel_size"
    print("  preset: VEHICLE applied, CUSTOM left values intact")


def _test_export_errors():
    import bpy
    import bmesh
    from cad_cfd_proxy.core import export
    from cad_cfd_proxy.core.errors import PipelineError

    bpy.ops.mesh.primitive_cube_add(size=2.0)
    cube = bpy.context.active_object
    props = bpy.context.scene.cad_cfd_proxy

    # No path → PipelineError (Export operator turns this into a clean report).
    props.export_target = "SNAPPY"
    props.export_path = ""
    raised = False
    try:
        export.write(bpy.context, cube, props)
    except PipelineError:
        raised = True
    assert raised, "export with no path should raise PipelineError"

    # FDS export of non-watertight geometry is hard-blocked.
    bm = bmesh.new()
    bm.from_mesh(cube.data)
    bm.faces.ensure_lookup_table()
    bmesh.ops.delete(bm, geom=[bm.faces[0]], context="FACES")
    bm.to_mesh(cube.data)
    bm.free()
    props.export_target = "FDS"
    props.export_path = os.path.join(tempfile.mkdtemp(prefix="cadcfd_"), "x.geom")
    raised = False
    try:
        export.write(bpy.context, cube, props)
    except PipelineError:
        raised = True
    assert raised, "FDS export of non-watertight geometry should be blocked"
    print("  export errors: no-path and non-watertight-FDS both raise PipelineError")


def main():
    import bpy
    import cad_cfd_proxy

    bpy.ops.wm.read_factory_settings(use_empty=True)
    cad_cfd_proxy.register()

    _test_preset()
    _test_export_errors()

    cad_cfd_proxy.unregister()
    print("review test: PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print("review test: FAIL: %s" % exc)
        raise
