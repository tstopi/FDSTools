# -*- coding: utf-8 -*-
"""Headless functional test for Phases 5-6 (clearance + morphological close).

Run with Blender::

    blender --background --factory-startup \
        --python cad_cfd_proxy/tests/test_morphology.py

* Clearance: dilating a cube grows its bounding box by ~2·distance.
* Close: two cubes separated by a gap smaller than 2·radius merge into one
  connected component.

Exits non-zero on failure for CI.
"""

import sys


def _bbox_dim(obj):
    xs = [v.co.x for v in obj.data.vertices]
    return max(xs) - min(xs)


def _component_count(mesh):
    import bmesh
    bm = bmesh.new()
    bm.from_mesh(mesh)
    bm.verts.ensure_lookup_table()
    seen = set()
    components = 0
    for seed in bm.verts:
        if seed.index in seen:
            continue
        components += 1
        stack = [seed]
        seen.add(seed.index)
        while stack:
            v = stack.pop()
            for e in v.link_edges:
                o = e.other_vert(v)
                if o.index not in seen:
                    seen.add(o.index)
                    stack.append(o)
    bm.free()
    return components


def _test_clearance():
    import bpy
    from cad_cfd_proxy.core import volume, surface

    bpy.ops.mesh.primitive_cube_add(size=2.0)
    cube = bpy.context.active_object
    props = bpy.context.scene.cad_cfd_proxy
    props.voxel_size = 0.1
    props.fill_volume = True

    grid = volume.mesh_to_sdf(bpy.context, [cube], props)
    base = surface.volume_to_mesh(bpy.context, grid, props)
    base_dim = _bbox_dim(base)

    grid2 = volume.mesh_to_sdf(bpy.context, [cube], props)
    grid2 = volume.dilate(bpy.context, grid2, 0.3, props)
    grown = surface.volume_to_mesh(bpy.context, grid2, props)
    grown_dim = _bbox_dim(grown)

    print("  clearance: bbox %.3f -> %.3f (+%.3f, expect ~+0.6)"
          % (base_dim, grown_dim, grown_dim - base_dim))
    delta = grown_dim - base_dim
    assert 0.3 < delta < 0.9, "dilation did not grow bbox as expected: +%.3f" % delta


def _test_close():
    import bpy
    import bmesh
    from cad_cfd_proxy.core import volume, surface

    # Two unit cubes with a 0.3 gap (right face of A at x=0.5, left of B at 0.8).
    def add_cube(cx):
        bpy.ops.mesh.primitive_cube_add(size=1.0, location=(cx, 0, 0))
        return bpy.context.active_object

    a = add_cube(0.0)
    b = add_cube(1.3)
    props = bpy.context.scene.cad_cfd_proxy
    props.voxel_size = 0.05
    props.fill_volume = True

    grid = volume.mesh_to_sdf(bpy.context, [a, b], props)
    before = surface.volume_to_mesh(bpy.context, grid, props)
    n_before = _component_count(before.data)

    grid = volume.mesh_to_sdf(bpy.context, [a, b], props)
    grid = volume.morphological_close(bpy.context, grid, 0.2, props)  # 2r=0.4 > 0.3
    after = surface.volume_to_mesh(bpy.context, grid, props)
    n_after = _component_count(after.data)

    print("  close: components %d -> %d (expect 2 -> 1)" % (n_before, n_after))
    assert n_before == 2, "expected two separate cubes before close, got %d" % n_before
    assert n_after == 1, "close did not bridge the gap: %d components" % n_after


def main():
    import bpy
    import cad_cfd_proxy

    bpy.ops.wm.read_factory_settings(use_empty=True)
    cad_cfd_proxy.register()

    _test_clearance()
    _test_close()

    cad_cfd_proxy.unregister()
    print("morphology test: PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print("morphology test: FAIL: %s" % exc)
        raise
