# -*- coding: utf-8 -*-
"""Phase 12 — validation.

Watertight/manifold should pass by construction (assert it); the meaningful
checks are self-intersection, normal consistency, and degenerate faces. FDS
export hard-requires all of these; snappyHexMesh is more tolerant.
"""


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


def check(obj, props):
    """Run validation on *obj*, return a :class:`ValidationReport`. TODO(phase-12)."""
    raise NotImplementedError("validate.check is a Phase 12 stub")
