# -*- coding: utf-8 -*-
"""Phase 7 — internal fluid volume.

Extract the negative space: domain SDF minus solid SDF, then cap inlets/outlets
so the result is closed at the boundaries (caps become inlet/outlet patches).
"""


def build_domain_grid(context, source, props):
    """SDF of the flow domain (user object, or padded bbox). TODO(phase-7)."""
    raise NotImplementedError("domain.build_domain_grid is a Phase 7 stub")


def extract_fluid_volume(context, source, solid_grid, props):
    """Return the fluid-volume grid = domain − solid, capped. TODO(phase-7)."""
    raise NotImplementedError("domain.extract_fluid_volume is a Phase 7 stub")
