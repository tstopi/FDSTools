# CAD → CFD Proxy Mesh (Blender addon)

Turns complex CAD assemblies into simplified, solver-ready proxy meshes for CFD.
Supports **external envelopes** (solid obstacle) and **internal fluid volumes**,
with export to **OpenFOAM / snappyHexMesh** and **FDS `&GEOM`**.

See [`../CAD_Proxy_Mesh_Development_Plan.md`](../CAD_Proxy_Mesh_Development_Plan.md)
for the full development plan.

## Status

**Phase 1 — scaffolding.** The addon registers, shows its panel and diagnostics,
and wires the pipeline together. The pipeline phases themselves are stubs that
raise `NotImplementedError` until their phase lands.

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
blender --background --factory-startup --python cad_cfd_proxy/tests/test_smoke.py
```
