# -*- coding: utf-8 -*-
"""UI registration."""

import bpy

from . import panel

_MODULE_CLASSES = (
    panel.CADCFD_PT_main,
    panel.CADCFD_PT_diagnostics,
)


def register():
    for cls in _MODULE_CLASSES:
        bpy.utils.register_class(cls)


def unregister():
    for cls in reversed(_MODULE_CLASSES):
        bpy.utils.unregister_class(cls)
