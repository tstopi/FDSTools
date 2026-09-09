# -*- coding: utf-8 -*-
"""Cancel operator — requests the running Generate modal to stop.

Sets a flag the Generate modal checks each tick; cancellation takes effect at
the next phase boundary. (ESC in the viewport also cancels.)
"""

from bpy.types import Operator


class CADCFD_OT_cancel(Operator):
    bl_idname = "cadcfd.cancel"
    bl_label = "Cancel"
    bl_description = "Cancel the running generation (takes effect between phases)"

    @classmethod
    def poll(cls, context):
        return context.scene.cad_cfd_proxy.is_running

    def execute(self, context):
        context.scene.cad_cfd_proxy.cancel_requested = True
        return {"FINISHED"}
