# -*- coding: utf-8 -*-
"""Operator registration."""

import bpy

from . import generate, preview, export, cancel

_MODULE_CLASSES = (
    generate.CADCFD_OT_generate,
    preview.CADCFD_OT_preview,
    export.CADCFD_OT_export,
    cancel.CADCFD_OT_cancel,
)


def register():
    for cls in _MODULE_CLASSES:
        bpy.utils.register_class(cls)


def unregister():
    for cls in reversed(_MODULE_CLASSES):
        bpy.utils.unregister_class(cls)
