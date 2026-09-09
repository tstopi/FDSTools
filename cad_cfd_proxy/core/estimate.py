# -*- coding: utf-8 -*-
"""Phase 14 — pre-generation estimates.

Estimates voxel count, memory, and face count before the (potentially heavy)
generate runs, and raises warnings for large jobs.

The memory model uses the **sparse narrow band**: OpenVDB only stores voxels
near the surface, so cost scales with surface area, not bounding-box volume.
Estimating from bbox volume (dense) overstates memory by orders of magnitude and
would scare users off jobs that are actually fine — so we report the sparse
figure and only use the dense figure to flag genuinely huge jobs.
"""

import bmesh

# Rough bytes per active voxel (grid value + tree/topology overhead).
_BYTES_PER_VOXEL = 8.0
# Warn above these thresholds.
_WARN_VOXELS = 5.0e7
_WARN_MEMORY_MB = 2048.0


def estimate(context, props):
    """Return an estimate dict for the current source selection and settings."""
    sources = _resolve_source_objects(context, props)
    area, bbox = _area_and_bounds(context, sources)
    voxel = max(props.voxel_size, 1e-9)

    band = max(getattr(props, "interior_band_width", 3.0), 1.0)
    # Sparse: (surface area / voxel²) surface voxels, thickened by the band.
    surface_voxels = area / (voxel * voxel)
    sparse_voxels = surface_voxels * band
    # Dense (bbox volume / voxel³) — only for the "this is huge" warning.
    bbox_vol = bbox[0] * bbox[1] * bbox[2]
    dense_voxels = bbox_vol / (voxel ** 3)

    est_faces = surface_voxels * 2.0  # ~2 triangles per surface voxel face
    est_memory_mb = sparse_voxels * _BYTES_PER_VOXEL / 1.0e6

    warnings = []
    if sparse_voxels > _WARN_VOXELS:
        warnings.append("high voxel count (~%.0f M) — consider a larger voxel size"
                        % (sparse_voxels / 1e6))
    if est_memory_mb > _WARN_MEMORY_MB:
        warnings.append("high memory (~%.0f MB)" % est_memory_mb)
    if area <= 0.0:
        warnings.append("no surface area found in the selected source")

    return {
        "voxel_size": voxel,
        "surface_area": area,
        "bbox": bbox,
        "sparse_voxels": sparse_voxels,
        "dense_voxels": dense_voxels,
        "est_faces": est_faces,
        "est_memory_mb": est_memory_mb,
        "warnings": warnings,
    }


def format_lines(est):
    """Human-readable summary lines for the preview operator/panel."""
    lines = [
        "Voxel size: %.4g m" % est["voxel_size"],
        "Bounding box: %.3g x %.3g x %.3g m" % tuple(est["bbox"]),
        "Est. voxels (sparse band): ~%.2g" % est["sparse_voxels"],
        "Est. memory: ~%.0f MB" % est["est_memory_mb"],
        "Est. faces: ~%.2g" % est["est_faces"],
    ]
    lines += ["Warning: " + w for w in est["warnings"]]
    return lines


# --- helpers ---------------------------------------------------------------

def _resolve_source_objects(context, props):
    # Reuse the same source resolution as Phase 2 without duplicating geometry.
    from . import collect
    try:
        return collect._resolve_sources(context, props)
    except Exception:  # noqa: BLE001 — estimate must never hard-fail
        return list(context.selected_objects)


def _area_and_bounds(context, sources):
    """World-space surface area and bounding-box dimensions of the evaluated
    source meshes. Reads evaluated meshes (modifiers applied) without leaving
    datablocks behind."""
    depsgraph = context.evaluated_depsgraph_get()
    inf = float("inf")
    lo = [inf, inf, inf]
    hi = [-inf, -inf, -inf]
    area = 0.0
    for obj in sources:
        if obj.type not in {"MESH", "CURVE", "SURFACE", "FONT", "META"}:
            continue
        eval_obj = obj.evaluated_get(depsgraph)
        try:
            mesh = eval_obj.to_mesh()
        except RuntimeError:
            continue
        if mesh is not None and mesh.polygons:
            bm = bmesh.new()
            bm.from_mesh(mesh)
            bm.transform(obj.matrix_world)
            area += sum(f.calc_area() for f in bm.faces)
            for v in bm.verts:
                for i in range(3):
                    lo[i] = min(lo[i], v.co[i])
                    hi[i] = max(hi[i], v.co[i])
            bm.free()
        eval_obj.to_mesh_clear()

    if lo[0] == inf:
        return 0.0, (0.0, 0.0, 0.0)
    return area, (hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2])
