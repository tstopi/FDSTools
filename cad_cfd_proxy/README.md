# CAD → CFD Proxy Mesh (Blender addon)

Turns complex CAD assemblies into simplified, solver-ready proxy meshes for CFD.
Supports **external envelopes** (solid obstacle) and **internal fluid volumes**,
with export to **OpenFOAM / snappyHexMesh** and **FDS `&GEOM`**.

See [`../CAD_Proxy_Mesh_Development_Plan.md`](../CAD_Proxy_Mesh_Development_Plan.md)
for the full development plan.

## Status

- **Phase 1 — scaffolding.** Done. Registers, shows panel + diagnostics, wires
  the pipeline together.
- **Phase 2 — input processing.** Done (`core/collect.py`). Realizes instances,
  applies modifiers, and converts curves/text/surfaces/metaballs to mesh via the
  dependency graph, into a `CADCFD_work` collection — originals untouched.
- **Phase 3 — geometry cleanup.** Done (`core/cleanup.py`). bmesh: merge by
  distance, dissolve degenerate, delete loose, drop zero-area faces, recalculate
  normals.
- **Phase 4 — Mesh→Volume (SDF).** Done (`core/volume.py`). Joins the cleaned
  meshes and runs a Mesh-to-Volume modifier into an OpenVDB grid. Note: defaults
  to **fill** (not band) so Volume→Mesh yields a solid envelope — a deliberate
  deviation from the plan's band default, since band-only gives double-walled
  shells. Phase 6 close makes fill robust on open shells.
- **Phase 8 — Volume→Mesh.** Done (`core/surface.py`). Volume-to-Mesh modifier
  on the SDF grid, baked to a real watertight/manifold mesh. **The EXTERNAL
  pipeline now runs end-to-end** (collect → cleanup → volume → surface).
- **Phase 10 — adaptive decimation.** Done (`core/surface.py`). Collapse-decimate
  with ratio bisection to hit a target triangle count within ±5%.
- **Phase 11 — export.** Done (`core/export.py`). snappyHexMesh **ASCII STL**
  (named solid) and **FDS `&GEOM`** (VERTS/FACES namelist), both triangulated and
  metre-scaled. Patch identity is lost through voxelization, so output is a
  single named solid/region. FDS watertight/manifold enforcement is a warning
  until Phase 12.
- **Phase 12 — validation.** Done (`core/validate.py`). bmesh watertight/
  manifold/normal-consistency/degenerate checks + BVH-tree self-intersection
  detection. `generate` reports the summary; **FDS export now hard-blocks** on
  non-solver-ready geometry (snappyHexMesh stays tolerant).
- **Phases 5–6 — clearance + morphological close.** Done (`core/volume.py`),
  API-only route: each op is a volume→mesh → offset-along-normals → mesh→volume
  round-trip (re-voxelization discards the self-intersections a raw offset
  creates). Clearance dilates; close = dilate then erode, bridging gaps and
  filling holes smaller than ~2·radius. Both no-op at zero distance.
- **Phase 7 — internal fluid volume.** Done (`core/domain.py`). Boolean
  DIFFERENCE of a padded bounding box (or a user domain object) minus the
  simplified solid → the flow domain around the object (wind-tunnel style, as
  snappyHexMesh meshes it). Enclosed hollow-part interiors (duct internals) are
  out of scope for v1.
- **Phase 9 — smoothing.** Done (`core/surface.py`). Laplacian vertex smoothing
  (iterations/strength); no-op at zero.
- **Phase 13 — compat hardening.** Done. `generate` refuses to run and operators
  disable with a poll message on an unsupported build (missing nodes/modifiers).
- **Phase 14 — preview/estimates.** Done (`core/estimate.py`). Voxel/memory/face
  estimates using the **sparse narrow-band** model (not bbox-volume), with
  large-job warnings; surfaced by the Preview operator.
- **Phase 15 — performance/cleanup.** Done. Temporary volume/mesh datablocks are
  freed after a run and the proxy is re-homed in the scene (`keep_intermediates`
  overrides).

All 15 phases implemented.

## Layout

| Path | Role |
|------|------|
| `__init__.py` | registration, `bl_info` |
| `blender_manifest.toml` | extension manifest (Blender 4.2+/5.x) |
| `compat.py` | version + node/modifier capability layer |
| `properties.py` | `Scene.cad_cfd_proxy` parameter group |
| `operators/` | generate · preview · export |
| `core/` | pipeline: collect · cleanup · volume · domain · surface · validate · export |
| `ui/` | sidebar panel + diagnostics |
| `utils/` | logging, temp-datablock cleanup |
| `tests/` | headless smoke test |

## Install

This ships with a `blender_manifest.toml`, so Blender 4.2+/5.x treats it as an
**extension** (not a legacy add-on). It has no external Python dependencies, so
there are no wheels to bundle. After enabling, the panel appears in the 3D
Viewport **N-panel ▸ CAD→CFD** tab, with a **Diagnostics** sub-panel showing the
capability checks.

**Option A — build a zip, then Install from Disk (recommended)**

```sh
# from the repo root; produces cad_cfd_proxy-<version>.zip
blender --command extension build --source-dir cad_cfd_proxy
```

Then in Blender: **Edit ▸ Preferences ▸ Get Extensions ▸ ▾ (top-right) ▸
Install from Disk…**, pick the zip, and enable it.

**Option B — drag-and-drop**

Zip the folder so `blender_manifest.toml` sits at the **root** of the archive
(zip the *contents* of `cad_cfd_proxy/`, not the parent folder), then drag the
`.zip` onto the Blender window and confirm.

**Option C — dev symlink (no rebuild while iterating)**

Link or copy the `cad_cfd_proxy` folder into your user extensions directory,
e.g. `~/.config/blender/<version>/extensions/user_default/cad_cfd_proxy` (paths
differ per OS), then refresh local extensions in Preferences.

> **Blender 5.2 note:** the manifest's `blender_version_min = "4.2.0"` permits
> 5.x, but that floor is conservative and has only been CI-tested on 4.2 LTS.
> If a 5.x release renamed a Mesh-to-Volume / Volume-to-Mesh node or modifier
> property, the **Diagnostics** panel will flag it and `Generate` refuses to run
> with a clear message rather than crashing — that's the place to reconcile.

## Test

```sh
# all tests in one process
blender --background --factory-startup --python-exit-code 1 \
    --python cad_cfd_proxy/tests/run_all.py

# or a single test
blender --background --factory-startup --python cad_cfd_proxy/tests/test_smoke.py
```

CI (`.github/workflows/ci.yml`) byte-compiles the addon and runs `run_all.py`
against Blender on every push.
