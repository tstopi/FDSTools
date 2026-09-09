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
- **Phases 5–7, 9, 13–15** — pending. Dilate/close no-op at zero distance;
  smoothing is a graceful no-op so a default run completes.
  Internal-mode (Phase 7) still raises until implemented.

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

## Install (dev)

Blender ≥ 4.2: **Edit ▸ Preferences ▸ Get Extensions ▸ Install from Disk** and
point at this folder (or a zip of it). Panel appears in the 3D viewport sidebar
under **CAD→CFD**.

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
