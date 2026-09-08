# -*- coding: utf-8 -*-
"""Phase 3 — geometry cleanup.

Preprocess the working meshes so voxelization behaves: weld coincident verts,
collapse degenerate geometry, drop loose (wire/point) geometry, and make face
normals consistent — the last matters for a correct SDF sign. Curve/text/etc.
conversion already happened in Phase 2 via the dependency graph, so there's
nothing to do here for that.

Operates in place on the working objects; originals were already left untouched
by Phase 2.
"""

import bmesh


def clean(objects, props):
    """Clean each working mesh object in *objects* in place."""
    for obj in objects:
        if obj.type == "MESH":
            _clean_object(obj, props)


def _clean_object(obj, props):
    mesh = obj.data
    bm = bmesh.new()
    bm.from_mesh(mesh)

    # 1) Weld coincident vertices. Do this first: it's what turns a pile of
    #    disconnected CAD triangles into a connected surface.
    if props.merge_distance > 0.0:
        bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=props.merge_distance)

    # 2) Collapse zero-length edges / degenerate faces left by welding.
    bmesh.ops.dissolve_degenerate(
        bm, dist=max(props.merge_distance, 1e-9), edges=bm.edges)

    # 3) Remove loose geometry (verts/edges bounding no face).
    if props.delete_loose:
        loose_verts = [v for v in bm.verts if not v.link_faces]
        if loose_verts:
            bmesh.ops.delete(bm, geom=loose_verts, context="VERTS")
        loose_edges = [e for e in bm.edges if not e.link_faces]
        if loose_edges:
            bmesh.ops.delete(bm, geom=loose_edges, context="EDGES")

    # 4) Drop any remaining zero-area faces.
    zero_area = [f for f in bm.faces if f.calc_area() <= 0.0]
    if zero_area:
        bmesh.ops.delete(bm, geom=zero_area, context="FACES")

    # 5) Make normals consistent (outward for closed islands).
    if props.fix_normals:
        bmesh.ops.recalc_face_normals(bm, faces=bm.faces)

    bm.to_mesh(mesh)
    bm.free()
    mesh.update()
