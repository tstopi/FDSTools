# -*- coding: utf-8 -*-
"""Generate operator — runs the full pipeline."""

import bpy
from bpy.types import Operator

from .. import compat, core


class CADCFD_OT_generate(Operator):
    bl_idname = "cadcfd.generate"
    bl_label = "Generate Envelope"
    bl_description = "Build a simplified CFD proxy mesh from the source assembly"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        missing = compat.missing_capabilities()
        if missing:
            cls.poll_message_set("Unsupported Blender build: " + ", ".join(missing))
            return False
        return True

    def execute(self, context):
        props = context.scene.cad_cfd_proxy
        try:
            core.generate_proxy(context, props, report=self.report)
        except NotImplementedError as exc:
            self.report({"WARNING"}, "Not yet implemented: %s" % exc)
            return {"CANCELLED"}
        except core.PipelineError as exc:
            self.report({"ERROR"}, str(exc))
            return {"CANCELLED"}
        self.report({"INFO"}, "Proxy generated")
        return {"FINISHED"}
