# -*- coding: utf-8 -*-
"""Sidebar panels (plan Phase 13 layout)."""

from bpy.types import Panel

from .. import compat

_CATEGORY = "CAD→CFD"


class CADCFD_PT_main(Panel):
    bl_label = "CAD → CFD Proxy"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = _CATEGORY

    def draw(self, context):
        layout = self.layout
        props = context.scene.cad_cfd_proxy

        layout.prop(props, "mode")

        box = layout.box()
        box.label(text="Source")
        box.prop(props, "source")
        if props.source == "NAMED_COLLECTION":
            box.prop(props, "source_collection")
        if props.mode == "INTERNAL":
            box.prop(props, "domain_object")
            if props.domain_object is None:
                box.prop(props, "domain_padding")

        box = layout.box()
        box.label(text="Volume")
        box.prop(props, "preset")
        col = box.column(align=True)
        col.prop(props, "voxel_size")
        if props.mode == "EXTERNAL":
            col.prop(props, "clearance")
        col.prop(props, "feature_size")

        box = layout.box()
        box.label(text="Surface")
        col = box.column(align=True)
        col.prop(props, "adaptivity")
        col.prop(props, "smooth_iterations")
        col.prop(props, "smooth_strength")
        col.prop(props, "target_faces_preset")
        if props.target_faces_preset == "CUSTOM":
            col.prop(props, "target_faces_custom")

        box = layout.box()
        box.label(text="Export")
        col = box.column(align=True)
        col.prop(props, "export_target")
        col.prop(props, "export_scale")
        col.prop(props, "patch_naming")
        col.prop(props, "export_path")

        row = layout.row(align=True)
        row.operator("cadcfd.preview", icon="ZOOM_ALL")
        row.operator("cadcfd.generate", icon="MOD_REMESH")
        layout.operator("cadcfd.export", icon="EXPORT")


class CADCFD_PT_diagnostics(Panel):
    bl_label = "Diagnostics"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = _CATEGORY
    bl_options = {"DEFAULT_CLOSED"}

    def draw(self, context):
        layout = self.layout
        for label, ok, detail in compat.check_capabilities():
            row = layout.row()
            row.label(text=label, icon="CHECKMARK" if ok else "ERROR")
            row.label(text=detail)
