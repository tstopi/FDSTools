# -*- coding: utf-8 -*-
"""Preview operator — pre-generation estimates (Phase 14)."""

from bpy.types import Operator

from .. import compat, core


class CADCFD_OT_preview(Operator):
    bl_idname = "cadcfd.preview"
    bl_label = "Preview"
    bl_description = ("Estimate voxel count, memory and face count before "
                      "generating (sparse-band estimate)")
    bl_options = {"REGISTER"}

    @classmethod
    def poll(cls, context):
        missing = compat.missing_capabilities()
        if missing:
            cls.poll_message_set("Unsupported Blender build: " + ", ".join(missing))
            return False
        return True

    def execute(self, context):
        props = context.scene.cad_cfd_proxy
        est = core.estimate.estimate(context, props)
        lines = core.estimate.format_lines(est)
        for line in lines:
            print("[cad_cfd_proxy] %s" % line)
        # Surface the headline (+ first warning) in the status bar.
        headline = "~%.2g faces, ~%.0f MB" % (est["est_faces"], est["est_memory_mb"])
        if est["warnings"]:
            self.report({"WARNING"}, "Preview: %s — %s"
                        % (headline, est["warnings"][0]))
        else:
            self.report({"INFO"}, "Preview: " + headline)
        return {"FINISHED"}
