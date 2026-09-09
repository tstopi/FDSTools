# -*- coding: utf-8 -*-
"""Phase 12 — validation.

Watertight + manifold usually pass by construction (the VDB isosurface); the
checks that actually catch problems are self-intersection and normal
consistency, introduced by smoothing/decimation. FDS ``&GEOM`` hard-requires
watertight + manifold + consistent orientation with no self-intersections;
snappyHexMesh is more tolerant.
"""

import bmesh
from mathutils.bvhtree import BVHTree

# BVH overlap epsilon (metres, pre-scale) for self-intersection detection.
_BVH_EPSILON = 1e-6
# Faces below this area are treated as degenerate.
_MIN_FACE_AREA = 1e-12


class ValidationReport(object):
    """Collected validation results for a proxy mesh."""

    def __init__(self):
        self.watertight = None
        self.manifold = None
        self.self_intersections = None
        self.normals_consistent = None
        self.degenerate_faces = None
        self.face_count = None

    @property
    def ok(self):
        checks = (self.watertight, self.manifold, self.normals_consistent)
        return all(c for c in checks if c is not None) and not self.self_intersections

    def fds_ready(self):
        """FDS &GEOM requires watertight + manifold + consistent orientation."""
        return bool(self.watertight and self.manifold and self.normals_consistent
                    and not self.self_intersections)

    def summary(self):
        parts = [
            "watertight" if self.watertight else "NOT watertight",
            "manifold" if self.manifold else "NOT manifold",
            "normals ok" if self.normals_consistent else "normals inconsistent",
        ]
        if self.self_intersections:
            parts.append("%d self-intersection(s)" % self.self_intersections)
        if self.degenerate_faces:
            parts.append("%d degenerate face(s)" % self.degenerate_faces)
        parts.append("%s faces" % self.face_count)
        return ", ".join(parts)


def check(obj, props=None):
    """Run validation on *obj* and return a :class:`ValidationReport`."""
    rep = ValidationReport()
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    bm.normal_update()
    bm.faces.ensure_lookup_table()
    bm.verts.ensure_lookup_table()

    rep.face_count = len(bm.faces)

    # Edge topology: every edge of a closed manifold surface bounds exactly two
    # faces. 1 = open boundary; 0 or >2 = non-manifold.
    boundary = 0
    nonmanifold_edges = 0
    for e in bm.edges:
        n = len(e.link_faces)
        if n == 1:
            boundary += 1
        elif n == 0 or n > 2:
            nonmanifold_edges += 1

    rep.watertight = (boundary == 0 and nonmanifold_edges == 0)
    rep.manifold = (nonmanifold_edges == 0
                    and all(v.is_manifold for v in bm.verts))

    # Consistent winding: every two-face (manifold) edge must be contiguous.
    rep.normals_consistent = all(
        e.is_contiguous for e in bm.edges if len(e.link_faces) == 2)

    rep.degenerate_faces = sum(
        1 for f in bm.faces if f.calc_area() <= _MIN_FACE_AREA)

    rep.self_intersections = _count_self_intersections(bm)

    bm.free()
    return rep


def _count_self_intersections(bm):
    """Count intersecting face pairs via BVH overlap, excluding adjacent faces.

    Faces that merely share a vertex/edge legitimately overlap within epsilon,
    so they're filtered out; what remains are genuine crossings.
    """
    tree = BVHTree.FromBMesh(bm, epsilon=_BVH_EPSILON)
    seen = set()
    count = 0
    for i, j in tree.overlap(tree):
        if i == j:
            continue
        key = (i, j) if i < j else (j, i)
        if key in seen:
            continue
        seen.add(key)
        vi = {v.index for v in bm.faces[i].verts}
        vj = {v.index for v in bm.faces[j].verts}
        if vi & vj:  # shared vertex/edge → adjacency, not a real intersection
            continue
        count += 1
    return count
