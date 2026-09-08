# -*- coding: utf-8 -*-
"""Export operator — writes the proxy to the configured CFD target."""

from bpy.types import Operator

from .. import core


class CADCFD_OT_export(Operator):
    bl_idname = "cadcfd.export"
    bl_label = "Export"
    bl_description = "Export the active proxy mesh to snappyHexMesh STL or FDS &GEOM"
    bl_options = {"REGISTER"}

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and obj.type == "MESH"

    def execute(self, context):
        props = context.scene.cad_cfd_proxy
        try:
            core.export_proxy(context, context.active_object, props, report=self.report)
        except NotImplementedError as exc:
            self.report({"WARNING"}, "Not yet implemented: %s" % exc)
            return {"CANCELLED"}
        self.report({"INFO"}, "Exported")
        return {"FINISHED"}
