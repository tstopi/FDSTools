# -*- coding: utf-8 -*-
"""Phase 11 — export (the CFD hand-off).

Common to all targets: triangulate, scale to metres, and rely on consistent
normals from cleanup. Then dispatch to the chosen writer.

Note on patch identity: voxelization merges everything into one surface, so the
per-source patch mapping from Phase 2 does not survive to here — the proxy is
written as a single named solid/region. Multi-patch output would need a
post-hoc spatial assignment, which is out of scope for v1.
"""

import bmesh
import bpy

from . import validate
from .errors import PipelineError


def write(context, proxy, props, report=None):
    """Export *proxy* to the configured target. Returns the output path."""
    path = bpy.path.abspath(props.export_path) if props.export_path else ""
    if not path:
        raise PipelineError("No export path set")

    name = _solid_name(proxy, props)
    verts, faces, normals = _triangulated(proxy, props.export_scale)
    if not faces:
        raise PipelineError("Nothing to export: proxy has no faces")

    if props.export_target == "SNAPPY":
        _write_stl_ascii(path, verts, faces, normals, name)
    elif props.export_target == "FDS":
        _check_fds_ready(proxy, props, report)
        _write_fds_geom(path, verts, faces, name)
    else:
        raise PipelineError("Unknown export target: %r" % props.export_target)

    if report:
        report({"INFO"}, "Exported %d triangles to %s" % (len(faces), path))
    return path


# --- shared prep -----------------------------------------------------------

def _triangulated(proxy, scale):
    """Return (verts, faces, normals): metre-scaled verts, 0-based triangle
    index triples, and per-face unit normals."""
    bm = bmesh.new()
    bm.from_mesh(proxy.data)
    bmesh.ops.triangulate(bm, faces=bm.faces)
    bm.normal_update()
    bm.verts.index_update()

    verts = [(v.co.x * scale, v.co.y * scale, v.co.z * scale) for v in bm.verts]
    faces, normals = [], []
    for f in bm.faces:
        faces.append(tuple(v.index for v in f.verts))
        normals.append((f.normal.x, f.normal.y, f.normal.z))
    bm.free()
    return verts, faces, normals


def _solid_name(proxy, props):
    if props.patch_naming == "SINGLE":
        return "proxy"
    # Per-object/collection identity is lost through voxelization; use the
    # proxy's own base name so the solid at least has a stable label.
    return proxy.name.replace(" ", "_") or "proxy"


def _check_fds_ready(proxy, props, report):
    """FDS &GEOM requires watertight+manifold+oriented geometry with no
    self-intersections. Fail loudly rather than write geometry FDS will reject.
    """
    rep = validate.check(proxy, props)
    if rep.fds_ready():
        return
    raise PipelineError(
        "FDS &GEOM export blocked — geometry is not solver-ready: %s. "
        "Fix the proxy (or export to snappyHexMesh STL, which is tolerant)."
        % rep.summary())


# --- writers ---------------------------------------------------------------

def _write_stl_ascii(path, verts, faces, normals, name):
    """ASCII STL — snappyHexMesh-ready. One named solid."""
    with open(path, "w") as fh:
        fh.write("solid %s\n" % name)
        for (a, b, c), n in zip(faces, normals):
            fh.write(" facet normal %e %e %e\n" % n)
            fh.write("  outer loop\n")
            for i in (a, b, c):
                fh.write("   vertex %e %e %e\n" % verts[i])
            fh.write("  endloop\n")
            fh.write(" endfacet\n")
        fh.write("endsolid %s\n" % name)


def _write_fds_geom(path, verts, faces, name, surf_id="INERT"):
    """FDS ``&GEOM`` namelist: flat VERTS (metres) + 1-based FACES with a
    trailing surface index referencing the single SURF_ID."""
    verts_rows = ",\n            ".join(
        "%.6f,%.6f,%.6f" % (x, y, z) for x, y, z in verts)
    faces_rows = ",\n            ".join(
        "%d,%d,%d,1" % (a + 1, b + 1, c + 1) for a, b, c in faces)
    with open(path, "w") as fh:
        fh.write("&GEOM ID='%s',\n" % name)
        fh.write("      SURF_ID='%s',\n" % surf_id)
        fh.write("      VERTS=%s,\n" % verts_rows)
        fh.write("      FACES=%s /\n" % faces_rows)
