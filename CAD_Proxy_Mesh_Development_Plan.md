# CAD → CFD Proxy Mesh Generator for Blender

Blender addon that turns complex CAD assemblies into simplified, solver-ready proxy
meshes for CFD — supporting both **external envelopes** (solid obstacle) and **internal
fluid volumes** (the cavity the fluid flows through).

Export targets:
- **OpenFOAM / snappyHexMesh** — named-patch STL, metre scale, tolerant of small leaks.
- **FDS `&GEOM`** — unstructured triangulated geometry for the cut-cell solver. FDS is
  **strict**: geometry must be **watertight, manifold, and consistently outward-oriented**
  in metres, or the cut-cell mesher rejects it. The VDB isosurface pipeline satisfies this
  by construction, which is the main reason this approach suits FDS well.

(Part of the FDSTools monorepo of CFD-workflow tooling.)

---

## 0. Guiding principles

- **The deliverable is an STL a CFD solver accepts**, not just a nice mesh in Blender.
  Units, normals, and patch names are first-class, not an afterthought.
- **Non-destructive core.** Originals are never modified. The generator runs on a
  duplicate and, where possible, as a re-runnable parametric node group.
- **The SDF is the workhorse.** Mesh-to-Volume produces a signed distance field. Almost
  every processing step (clearance, feature suppression, hole filling) is an SDF offset
  or morphological operation — cheaper and more predictable than repeated remeshing.
- **Be honest about accuracy.** Level-set remeshing rounds sharp edges and drops
  sub-voxel features. Surface the trade-off; don't hide it.

---

## Pipeline overview

```text
                    ┌─────────────────────────────┐
Input assembly ───► │ Collect · realize · to-mesh │
                    └──────────────┬──────────────┘
                                   ▼
                    ┌─────────────────────────────┐
                    │ Cleanup (merge, normals,     │
                    │ curves→mesh, loose geom)     │
                    └──────────────┬──────────────┘
                                   ▼
                    ┌─────────────────────────────┐
                    │ Mesh → Volume  (SDF)         │
                    └──────────────┬──────────────┘
                                   ▼
              ┌────────────────────┴────────────────────┐
   EXTERNAL   ▼                                          ▼   INTERNAL
 ┌────────────────────────┐              ┌──────────────────────────────────┐
 │ Dilate (clearance)     │              │ Domain SDF − solid SDF            │
 │ Morph-close (suppress  │              │ (extract negative space)          │
 │ internals + fill holes)│              │ + inlet/outlet caps               │
 └───────────┬────────────┘              └────────────────┬─────────────────┘
             └───────────────────┬────────────────────────┘
                                 ▼
                    ┌─────────────────────────────┐
                    │ Volume → Mesh (adaptive)     │
                    └──────────────┬──────────────┘
                                   ▼
                    ┌─────────────────────────────┐
                    │ Smooth · decimate to target  │
                    └──────────────┬──────────────┘
                                   ▼
                    ┌─────────────────────────────┐
                    │ Validate · orient · scale ·  │
                    │ name patches · EXPORT (STL)  │
                    └─────────────────────────────┘
```

---

## Phase 1 — Core architecture

```text
cad_cfd_proxy/
├── __init__.py            # registration, bl_info
├── preferences.py
├── compat.py             # version + API capability layer (see Phase 13)
├── core/
│   ├── collect.py        # input gathering / realize / convert
│   ├── cleanup.py
│   ├── volume.py         # mesh→volume, SDF offsets, morphology
│   ├── domain.py         # internal-mode domain + boolean
│   ├── surface.py        # volume→mesh, smooth, decimate
│   ├── validate.py
│   └── export.py         # units, orientation, patch naming, STL/OBJ
├── nodes/                # parametric Geometry Nodes group builders
├── operators/
├── ui/
└── tests/                # headless bpy tests
```

Goals: modular, re-runnable, version-guarded.

**Recommendation:** implement the volume→surface core as a single **Geometry Nodes node
group** driven by the panel, rather than a baked stack of modifiers. It's
non-destructive, re-runnable on parameter change, and reuses Blender's own
Mesh-to-Volume / Volume-to-Mesh nodes. Fall back to modifiers only where a needed
operation has no node equivalent.

---

## Phase 2 — Input processing

Sources: selected mesh objects · active collection · named collection · nested
collections.

- Duplicate source; **never touch originals**.
- Realize instances / collection instances.
- Convert curves, text, metaballs, surfaces to mesh; apply modifiers.
- Skip/flag unsupported types with a report, don't fail the whole run.

---

## Phase 3 — Geometry cleanup

- Merge by distance, delete loose verts/edges.
- Recalculate normals outward (external mode) — needed for a correct SDF sign.
- Triangulate is **not** required here (VDB handles it); defer triangulation to export.
- Remove zero-area faces / degenerate geometry.

---

## Phase 4 — Volume creation (SDF)

`Mesh → Volume` producing a signed distance field.

**Density mode matters:** CAD assemblies are full of open shells. Default to
**band/exterior** density (offset around surfaces), not **interior fill** — interior
fill silently under-fills on non-watertight input. Offer interior fill only when the
source is known-closed.

User controls: **voxel size** (drives everything downstream), interior band width.
Presets: mechanical · product · vehicle · factory (each sets sensible voxel/clearance
defaults).

---

## Phase 5 — Envelope clearance (external mode)

Add a configurable offset around the source = **dilate the SDF** (shift the isosurface
outward by N metres). Cheaper and smoother than re-voxelizing.

Parameters: clearance distance, optional exterior band width.

---

## Phase 6 — Internal feature suppression + hole filling (morphological close)

One operation, not two: **dilate then erode** the SDF (morphological *close*) by a
"feature size" radius.

- Fills openings and gaps **smaller** than the radius (bolt holes, cable ports, vents).
- Removes thin internal detail (bolts, brackets, hidden internals) that the outer
  surface never sees.
- Preserves openings/features **larger** than the radius.

This replaces the plan's separate "internal suppression" and "hole filling" phases and
the double mesh↔volume round-trip — same result, one SDF pass.

Control: single **feature size** slider (radius). Optional separate open/close radii for
advanced use.

---

## Phase 7 — Internal fluid volume (internal mode)

For extracting the domain the fluid flows through:

1. **Domain source:** bounding box, padded bbox, or a user-supplied domain object.
2. Compute domain SDF and **subtract** the solid SDF → negative space.
3. **Cap** inlets/outlets: user marks faces/planes where the fluid domain is cut open, so
   the result is closed at boundaries (these become inlet/outlet patches).
4. Same smooth/decimate/export tail.

This is a genuinely separate branch from the external envelope; the UI exposes a mode
switch (External / Internal).

---

## Phase 8 — Surface reconstruction

`Volume → Mesh` (adaptive). Output is **watertight and manifold by construction** (it's
an isosurface of a level set). Controls: adaptivity, surface detail (grid-relative).

---

## Phase 9 — Surface smoothing

Smooth position → optional Laplacian smooth. Controls: iterations, strength.

**Caveat to surface in UI:** smoothing + decimation are where **self-intersections** can
appear (not leaks). Keep smoothing conservative and re-check in Phase 12.

---

## Phase 10 — Adaptive decimation to target face count

Targets: 5k · 10k · 25k · 50k · 100k (and custom). Bisection on collapse ratio →
face count (monotonic, converges in a few passes). Target within ±5%.

Note: snappyHexMesh resamples the surface, so low poly is about **STL size / speed**, not
solver necessity. Planar/collapse decimation; preserve boundary/feature edges where
possible.

---

## Phase 11 — Export (the CFD hand-off) **[new, critical]**

Common to all targets:
- **Units:** convert scene units → **metres** (explicit scale factor, shown in UI).
- **Orientation:** enforce consistent normals (outward for solids / toward fluid for
  internal); flag if inconsistent.
- **Triangulate** at export (VDB output is already tri, but guarantee it).

**OpenFOAM / snappyHexMesh:**
- ASCII or binary STL; OBJ as secondary.
- Name solids per object/collection so snappy sees named patches; option for multi-region
  STL or **one file per patch**. Inlet/outlet caps (Phase 7) become named patches.
- Optional starter `snappyHexMeshDict` snippet + suggested `locationInMesh` point.

**FDS `&GEOM`:**
- Emit geometry FDS accepts — inline `VERTS`/`FACES`, or a companion binary/`bingeom`
  file referenced from the namelist, depending on target FDS version.
- **Hard-enforce watertight + manifold + consistent orientation** before writing (fail
  the export with a clear report if validation from Phase 12 doesn't pass — FDS will
  otherwise reject or misbehave silently).
- Map object/collection → `SURF_ID` so materials/boundary conditions can be assigned.
- Optional: stub `&GEOM` and matching `&SURF` lines for the input file.

---

## Phase 12 — Validation suite

- Watertight / manifold check (should pass by construction — assert it).
- **Self-intersection** check after smoothing/decimation (the real failure mode).
- Face-count compliance (±5%).
- Normal consistency.
- Non-zero, non-degenerate faces.
- **Headless harness:** `blender --background --python tests/…` over mechanical /
  electronics / vehicle / plant fixtures, runnable in CI.

---

## Phase 13 — Compatibility layer

- **Verify the actual target Blender version and pin its API** (don't assume a version
  number; confirm Mesh-to-Volume / Volume-to-Mesh node identifiers and socket names on
  the real build).
- Runtime capability checks: GeometryNodeTree, Mesh-to-Volume node, Volume-to-Mesh node,
  Nodes modifier, Decimate.
- Dynamic socket lookup by name, not index. Startup diagnostics panel.

---

## Phase 14 — Preview & estimates

Before generation: **voxel count**, memory, face-count, warnings.

**Fix the memory model:** VDB is a *sparse narrow band*, so estimate from
`surface_area × band_width / voxel³`, **not** `bbox_volume / voxel³`. The volume estimate
overstates memory by orders of magnitude and will scare users off viable jobs.

---

## Phase 15 — Performance

Large-assembly support · temporary datablock cleanup · node-group reuse · memory recovery
· progress reporting (modal operator with cancel).

---

## UI

```text
CAD → CFD Proxy
─────────────────────────────
Mode:  ( ) External envelope
       ( ) Internal fluid volume
Source:            [ ... ]
Domain (internal): [ ... ]
Preset:            [ Mechanical ▾ ]
Voxel size
Clearance / feature size
Target faces
Smooth iterations
─────────────────────────────
Export units:  [ mm → m ]
Patch naming:  [ per-object ▾ ]
─────────────────────────────
[ Preview ]   [ Generate ]   [ Export STL ]
```

---

## v1.0 deliverables

- External envelope **and** internal fluid-volume modes
- SDF-based clearance + morphological close (suppression + hole fill)
- Adaptive Volume→Mesh, smoothing, target-face decimation (±5%)
- **Export to snappyHexMesh STL (named patches) and FDS `&GEOM`**, metre scaling,
  consistent orientation, watertight/manifold enforced for FDS
- Preview with corrected (sparse) memory estimate
- Version-pinned compatibility + diagnostics
- Presets, progress, watertight/manifold + self-intersection validation
- Headless test harness
