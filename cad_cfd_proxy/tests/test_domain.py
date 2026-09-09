# -*- coding: utf-8 -*-
"""Headless functional test for Phase 7 (internal fluid volume).

Run with Blender::

    blender --background --factory-startup \
        --python cad_cfd_proxy/tests/test_domain.py

A cube obstacle inside a padded bounding box should yield a flow domain with two
closed shells (outer box + inner obstacle cavity) spanning the padded box.
Exits non-zero on failure for CI.
"""

import sys


def _component_count(mesh):
    import bmesh
    bm = bmesh.new()
    bm.from_mesh(mesh)
    bm.verts.ensure_lookup_table()
    bm.verts.index_update()
    seen = set()
    n = 0
    for seed in bm.verts:
        if seed.index in seen:
            continue
        n += 1
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
    return n


def _bbox_dim(obj):
    xs = [v.co.x for v in obj.data.vertices]
    return max(xs) - min(xs)


def main():
    import bpy
    import cad_cfd_proxy
    from cad_cfd_proxy.core import volume, domain

    bpy.ops.wm.read_factory_settings(use_empty=True)
    cad_cfd_proxy.register()

    bpy.ops.mesh.primitive_cube_add(size=2.0, location=(0, 0, 0))  # obstacle, span 2
    obstacle = bpy.context.active_object

    props = bpy.context.scene.cad_cfd_proxy
    props.mode = "INTERNAL"
    props.voxel_size = 0.1
    props.fill_volume = True
    props.domain_padding = 1.0   # box spans ~4 m

    grid = volume.mesh_to_sdf(bpy.context, [obstacle], props)
    fluid = domain.extract_fluid_volume(bpy.context, [obstacle], grid, props)

    assert fluid.type == "MESH", "fluid domain is not a mesh"
    n_faces = len(fluid.data.polygons)
    assert n_faces > 0, "fluid domain has no faces"

    # Two closed shells: outer box boundary + inner obstacle cavity.
    n_comp = _component_count(fluid.data)
    dim = _bbox_dim(fluid)
    print("  fluid domain: %d faces, %d shells, bbox %.2f m" % (n_faces, n_comp, dim))
    assert n_comp == 2, "expected 2 shells (box + cavity), got %d" % n_comp
    # Domain spans the padded box (~4 m), clearly larger than the 2 m obstacle.
    assert dim > 3.0, "domain bbox too small (%.2f), padding not applied?" % dim

    # Full pipeline in INTERNAL mode (orchestrator branch + decimate + validate).
    from cad_cfd_proxy import core
    for o in list(bpy.data.objects):
        o.select_set(o is obstacle)
    bpy.context.view_layer.objects.active = obstacle
    props.source = "SELECTED"
    proxy = core.generate_proxy(bpy.context, props)
    assert proxy is not None and len(proxy.data.polygons) > 0, \
        "internal generate_proxy produced no geometry"
    print("  internal generate_proxy: %d faces" % len(proxy.data.polygons))

    cad_cfd_proxy.unregister()
    print("domain test: PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print("domain test: FAIL: %s" % exc)
        raise
