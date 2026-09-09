# -*- coding: utf-8 -*-
"""Pipeline orchestration.

The individual phases live in sibling modules (:mod:`.collect`, :mod:`.cleanup`,
:mod:`.volume`, :mod:`.domain`, :mod:`.surface`, :mod:`.validate`,
:mod:`.export`).

Two entry points chain them:

* :func:`generate_job` is a generator that performs one phase per ``next()`` and
  yields ``(label, fraction)`` progress. It lets the modal operator drive the
  pipeline a step at a time, drawing a progress bar and honouring cancel between
  phases (the heavy per-phase C ops can't be interrupted mid-step).
* :func:`generate_proxy` exhausts that generator synchronously — used by scripts
  and the headless tests.
"""

from .. import compat
from . import collect, cleanup, volume, domain, surface, validate, export, estimate
from .errors import PipelineError  # noqa: F401 (re-exported)


def generate_job(context, props):
    """Generator running the pipeline one phase at a time.

    Yields ``(label, fraction)`` after each phase completes; returns the proxy
    object (via ``StopIteration.value``). Raises :class:`PipelineError` /
    :class:`NotImplementedError` like the phase functions do.
    """
    missing = compat.missing_capabilities()
    if missing:
        raise PipelineError(
            "Unsupported Blender build; missing: %s" % ", ".join(missing))

    sources = collect.gather(context, props)
    yield ("Collected source geometry", 0.15)

    cleanup.clean(sources, props)
    yield ("Cleaned geometry", 0.30)

    grid = volume.mesh_to_sdf(context, sources, props)
    yield ("Built volume (SDF)", 0.45)

    if props.mode == "EXTERNAL":
        grid = volume.dilate(context, grid, props.clearance, props)
        yield ("Applied clearance", 0.55)
        grid = volume.morphological_close(context, grid, props.feature_size, props)
        yield ("Closed features", 0.65)
        proxy = surface.volume_to_mesh(context, grid, props)
        yield ("Reconstructed surface", 0.78)
    else:  # INTERNAL
        proxy = domain.extract_fluid_volume(context, sources, grid, props)
        yield ("Extracted fluid volume", 0.78)

    surface.smooth(proxy, props)
    yield ("Smoothed", 0.85)

    surface.decimate_to_target(context, proxy, props.resolve_target_faces())
    yield ("Decimated to target", 0.92)

    rep = validate.check(proxy, props)
    yield ("Validated: %s" % rep.summary(), 0.97)

    finalize_proxy(context, proxy, props)
    yield ("Done", 1.0)
    return proxy


def generate_proxy(context, props, report=None):
    """Run the full pipeline synchronously and return the proxy object.

    Exhausts :func:`generate_job`, forwarding each phase label to *report*.
    """
    gen = generate_job(context, props)
    try:
        while True:
            label, _fraction = next(gen)
            if report:
                report({"INFO"}, label)
    except StopIteration as stop:
        return stop.value


def discard_work(context):
    """Remove the whole ``CADCFD_work`` collection and its datablocks.

    Used to clean up partially-built intermediates when a run is cancelled.
    """
    import bpy

    work = bpy.data.collections.get(collect.WORK_COLLECTION)
    if work is None:
        return
    for obj in list(work.objects):
        data = obj.data
        bpy.data.objects.remove(obj, do_unlink=True)
        _remove_orphan_data(data)
    if not work.objects and not work.children:
        bpy.data.collections.remove(work)


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
