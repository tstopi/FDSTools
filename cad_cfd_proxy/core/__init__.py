# -*- coding: utf-8 -*-
"""Pipeline orchestration.

The individual phases live in sibling modules (:mod:`.collect`, :mod:`.cleanup`,
:mod:`.volume`, :mod:`.domain`, :mod:`.surface`, :mod:`.validate`,
:mod:`.export`). :func:`generate_proxy` chains them; each is a no-op stub in the
Phase 1 scaffold and raises :class:`NotImplementedError` until its phase lands.
"""

from .. import compat
from . import collect, cleanup, volume, domain, surface, validate, export, estimate
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

    missing = compat.missing_capabilities()
    if missing:
        raise PipelineError(
            "Unsupported Blender build; missing: %s" % ", ".join(missing))

    _say("Collecting source geometry")
    sources = collect.gather(context, props)

    _say("Cleaning geometry")
    cleanup.clean(sources, props)

    _say("Building volume (SDF)")
    grid = volume.mesh_to_sdf(context, sources, props)

    if props.mode == "EXTERNAL":
        grid = volume.dilate(context, grid, props.clearance, props)
        grid = volume.morphological_close(context, grid, props.feature_size, props)
        _say("Reconstructing surface")
        proxy = surface.volume_to_mesh(context, grid, props)
    else:  # INTERNAL
        _say("Extracting fluid volume")
        proxy = domain.extract_fluid_volume(context, sources, grid, props)

    surface.smooth(proxy, props)
    surface.decimate_to_target(context, proxy, props.resolve_target_faces())

    _say("Validating")
    rep = validate.check(proxy, props)
    if report:
        level = "INFO" if rep.ok else "WARNING"
        report({level}, "Validation: %s" % rep.summary())

    _say("Cleaning up")
    finalize_proxy(context, proxy, props)
    return proxy


def finalize_proxy(context, proxy, props):
    """Phase 15 — move *proxy* to the scene and free intermediates.

    The pipeline stages leave temporary objects (joined mesh, volumes, obstacle,
    domain) in the ``CADCFD_work`` collection. Unless ``keep_intermediates`` is
    set, remove them and their datablocks to recover memory, and re-home the
    proxy in the scene's master collection.
    """
    import bpy

    work = bpy.data.collections.get(collect.WORK_COLLECTION)
    if work is None:
        return

    # Re-home the proxy in the scene master collection.
    if proxy.name in work.objects:
        work.objects.unlink(proxy)
    if proxy.name not in context.scene.collection.objects:
        context.scene.collection.objects.link(proxy)

    if props.keep_intermediates:
        return

    for obj in list(work.objects):
        data = obj.data
        bpy.data.objects.remove(obj, do_unlink=True)
        _remove_orphan_data(data)

    if not work.objects and not work.children:
        bpy.data.collections.remove(work)


def _remove_orphan_data(data):
    """Remove a mesh/volume datablock if nothing references it any more."""
    import bpy

    if data is None or data.users:
        return
    if isinstance(data, bpy.types.Mesh):
        bpy.data.meshes.remove(data)
    elif isinstance(data, bpy.types.Volume):
        bpy.data.volumes.remove(data)


def export_proxy(context, proxy, props, report=None):
    """Export *proxy* to the configured CFD target."""
    return export.write(context, proxy, props, report=report)
