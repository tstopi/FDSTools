# -*- coding: utf-8 -*-
"""CAD → CFD Proxy Mesh generator for Blender.

Turns complex CAD assemblies into simplified, solver-ready proxy meshes for
CFD, supporting external envelopes and internal fluid volumes, with export to
OpenFOAM/snappyHexMesh and FDS ``&GEOM``.

Phase 1 scaffolding: modular package, registration wiring, and a version /
capability compatibility layer. Pipeline phases live in :mod:`.core` and are
filled in by later phases.
"""

# bl_info is retained for the legacy add-on install path. On Blender 4.2+/5.x
# the extension system reads blender_manifest.toml instead and ignores this.
bl_info = {
    "name": "CAD → CFD Proxy Mesh",
    "author": "tstopi",
    "version": (0, 1, 0),
    "blender": (4, 2, 0),
    "location": "View3D > Sidebar > CAD→CFD",
    "description": "Simplified proxy meshes from CAD assemblies for OpenFOAM and FDS",
    "category": "Mesh",
}

# Submodules are registered in this order; unregistered in reverse.
from . import compat, properties, operators, ui

_MODULES = (properties, operators, ui)


def register():
    ok, detail = compat.is_supported()
    if not ok:
        # Register nothing meaningful on an unsupported build; surface why.
        print("[cad_cfd_proxy] unsupported Blender build: %s" % detail)
    compat.run_diagnostics(log=True)
    for mod in _MODULES:
        mod.register()


def unregister():
    for mod in reversed(_MODULES):
        mod.unregister()
