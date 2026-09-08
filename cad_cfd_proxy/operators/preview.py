# -*- coding: utf-8 -*-
"""Preview operator — estimates before generation (plan Phase 14)."""

from bpy.types import Operator

from .. import compat


class CADCFD_OT_preview(Operator):
    bl_idname = "cadcfd.preview"
    bl_label = "Preview"
    bl_description = ("Estimate voxel count, memory and face count before "
                      "generating (sparse-band estimate)")
    bl_options = {"REGISTER"}

    @classmethod
    def poll(cls, context):
        return compat.all_capabilities_ok()

    def execute(self, context):
        # TODO(phase-14): estimate from surface_area * band / voxel**3
        # (sparse), NOT bbox_volume / voxel**3 which grossly overstates memory.
        self.report({"WARNING"}, "Preview estimates are a Phase 14 stub")
        return {"CANCELLED"}
