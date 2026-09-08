# -*- coding: utf-8 -*-
"""Phase 11 — export (the CFD hand-off).

Common to all targets: convert to metres, enforce consistent orientation,
triangulate. Then dispatch to the chosen writer.
"""


def write(context, proxy, props, report=None):
    """Export *proxy* to the configured target. TODO(phase-11)."""
    if props.export_target == "SNAPPY":
        return write_snappy_stl(context, proxy, props, report=report)
    if props.export_target == "FDS":
        return write_fds_geom(context, proxy, props, report=report)
    raise ValueError("unknown export target: %r" % props.export_target)


def write_snappy_stl(context, proxy, props, report=None):
    """STL with named solids/patches, metre-scaled. TODO(phase-11).

    snappyHexMesh tolerates small leaks, so patch naming and orientation matter
    more than strict watertightness.
    """
    raise NotImplementedError("export.write_snappy_stl is a Phase 11 stub")


def write_fds_geom(context, proxy, props, report=None):
    """FDS ``&GEOM`` geometry, metre-scaled. TODO(phase-11).

    Hard-enforce watertight + manifold + consistent orientation before writing;
    fail loudly otherwise (FDS otherwise rejects or misbehaves silently). Map
    source object/collection → SURF_ID.
    """
    raise NotImplementedError("export.write_fds_geom is a Phase 11 stub")
