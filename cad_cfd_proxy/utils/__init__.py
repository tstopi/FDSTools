# -*- coding: utf-8 -*-
"""Shared helpers: logging, temporary-datablock tracking/cleanup (Phase 15)."""

import contextlib

PREFIX = "[cad_cfd_proxy]"


def log(msg):
    print("%s %s" % (PREFIX, msg))


@contextlib.contextmanager
def temp_datablocks(keep=False):
    """Track datablocks created during the block and remove them on exit.

    Usage::

        with temp_datablocks() as tmp:
            tmp.append(some_object)

    TODO(phase-15): wire into the pipeline for reliable cleanup / memory
    recovery. When *keep* is True, intermediates are left in the scene.
    """
    created = []
    try:
        yield created
    finally:
        if not keep:
            # TODO(phase-15): remove tracked datablocks via bpy.data.*.remove
            pass
