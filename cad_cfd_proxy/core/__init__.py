# -*- coding: utf-8 -*-
"""Pipeline orchestration.

The individual phases live in sibling modules (:mod:`.collect`, :mod:`.cleanup`,
:mod:`.volume`, :mod:`.domain`, :mod:`.surface`, :mod:`.validate`,
:mod:`.export`). :func:`generate_proxy` chains them; each is a no-op stub in the
Phase 1 scaffold and raises :class:`NotImplementedError` until its phase lands.
"""

from . import collect, cleanup, volume, domain, surface, validate, export
from .errors import PipelineError  # noqa: F401 (re-exported)


def generate_proxy(context, props, report=None):
    """Run the full generation pipeline and return the proxy object.

    Parameters
    ----------
    context : bpy.types.Context
    props : CADCFDProxyProperties
    report : callable, optional
        ``operator.report``-style callback for user-facing progress/warnings.

    Returns
    -------
    bpy.types.Object
        The generated proxy mesh (not yet exported).
    """
    def _say(msg):
        if report:
            report({"INFO"}, msg)

    _say("Collecting source geometry")
    sources = collect.gather(context, props)

    _say("Cleaning geometry")
    cleanup.clean(sources, props)

    _say("Building volume (SDF)")
    grid = volume.mesh_to_sdf(context, sources, props)

    if props.mode == "EXTERNAL":
        volume.dilate(grid, props.clearance)
        volume.morphological_close(grid, props.feature_size)
    else:  # INTERNAL
        grid = domain.extract_fluid_volume(context, source, grid, props)

    _say("Reconstructing surface")
    proxy = surface.volume_to_mesh(context, grid, props)
    surface.smooth(proxy, props)
    surface.decimate_to_target(proxy, props.resolve_target_faces())

    _say("Validating")
    validate.check(proxy, props)

    return proxy


def export_proxy(context, proxy, props, report=None):
    """Export *proxy* to the configured CFD target."""
    return export.write(context, proxy, props, report=report)
