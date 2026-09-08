# -*- coding: utf-8 -*-
"""Scene-level property group holding all generator parameters.

Attached to ``Scene.cad_cfd_proxy`` so settings persist in the .blend and drive
both the UI panel and the operators. Grouped to mirror the pipeline phases.
"""

import bpy
from bpy.props import (
    BoolProperty,
    EnumProperty,
    FloatProperty,
    IntProperty,
    PointerProperty,
    StringProperty,
)
from bpy.types import PropertyGroup


MODE_ITEMS = (
    ("EXTERNAL", "External envelope", "Wrap the outer surface (solid obstacle)"),
    ("INTERNAL", "Internal fluid volume", "Extract the cavity the fluid flows through"),
)

SOURCE_ITEMS = (
    ("SELECTED", "Selected objects", "Use the current selection"),
    ("ACTIVE_COLLECTION", "Active collection", "Use the active collection"),
    ("NAMED_COLLECTION", "Named collection", "Use a named collection"),
)

PRESET_ITEMS = (
    ("CUSTOM", "Custom", "Use the values below as-is"),
    ("MECHANICAL", "Mechanical", "Fine voxel, small clearance"),
    ("PRODUCT", "Product / electronics", "Medium voxel"),
    ("VEHICLE", "Vehicle", "Coarser voxel, larger clearance"),
    ("FACTORY", "Factory / plant", "Coarse voxel, large features"),
)

TARGET_FACE_ITEMS = (
    ("5000", "5k", ""),
    ("10000", "10k", ""),
    ("25000", "25k", ""),
    ("50000", "50k", ""),
    ("100000", "100k", ""),
    ("CUSTOM", "Custom", "Use the custom face count below"),
)

EXPORT_TARGET_ITEMS = (
    ("SNAPPY", "OpenFOAM / snappyHexMesh", "Named-patch STL, tolerant of small leaks"),
    ("FDS", "FDS &GEOM", "Watertight/manifold triangulated geometry, enforced"),
)

PATCH_NAMING_ITEMS = (
    ("PER_OBJECT", "Per object", "One patch/region per source object"),
    ("PER_COLLECTION", "Per collection", "One patch/region per source collection"),
    ("SINGLE", "Single", "One patch for the whole proxy"),
)


class CADCFDProxyProperties(PropertyGroup):
    # --- mode / input (phases 2, 7) ---
    mode: EnumProperty(name="Mode", items=MODE_ITEMS, default="EXTERNAL")
    source: EnumProperty(name="Source", items=SOURCE_ITEMS, default="SELECTED")
    source_collection: StringProperty(name="Collection")
    domain_object: PointerProperty(
        name="Domain", type=bpy.types.Object,
        description="Domain bounds for internal fluid-volume mode (empty = padded bbox)")
    domain_padding: FloatProperty(
        name="Domain padding", default=0.1, min=0.0, unit="LENGTH",
        description="Padding around the bounding box when no domain object is set")

    # --- cleanup (phase 3) ---
    merge_distance: FloatProperty(
        name="Merge distance", default=1e-4, min=0.0, soft_max=0.01, unit="LENGTH",
        description="Weld vertices closer than this (0 disables)")
    delete_loose: BoolProperty(
        name="Delete loose", default=True,
        description="Remove vertices/edges not bounding any face")
    fix_normals: BoolProperty(
        name="Recalculate normals", default=True,
        description="Make face normals consistent/outward (needed for correct SDF sign)")

    # --- volume / SDF (phases 4, 5, 6) ---
    preset: EnumProperty(name="Preset", items=PRESET_ITEMS, default="CUSTOM")
    voxel_size: FloatProperty(
        name="Voxel size", default=0.02, min=1e-5, soft_max=1.0, unit="LENGTH")
    fill_volume: BoolProperty(
        name="Fill volume", default=True,
        description="Fill the interior so Volume→Mesh yields a solid envelope "
                    "(needs reasonably closed input; Phase 6 close helps). "
                    "Disable for a thin band around open surfaces")
    interior_band_width: FloatProperty(
        name="Interior band", default=3.0, min=0.5, soft_max=10.0,
        description="SDF band thickness inside the surface, in voxels")
    clearance: FloatProperty(
        name="Clearance", default=0.0, min=0.0, unit="LENGTH",
        description="Outward SDF offset around the source (external mode)")
    feature_size: FloatProperty(
        name="Feature size", default=0.0, min=0.0, unit="LENGTH",
        description="Morphological close radius: fills holes / removes internals "
                    "smaller than this")

    # --- surface (phases 8, 9, 10) ---
    surface_threshold: FloatProperty(
        name="Iso threshold", default=0.5, min=0.0, max=1.0,
        description="Grid value at the extracted isosurface (Volume→Mesh)")
    adaptivity: FloatProperty(name="Adaptivity", default=0.0, min=0.0, max=1.0)
    smooth_iterations: IntProperty(name="Smooth iterations", default=2, min=0, max=100)
    smooth_strength: FloatProperty(name="Smooth strength", default=0.5, min=0.0, max=1.0)
    target_faces_preset: EnumProperty(
        name="Target faces", items=TARGET_FACE_ITEMS, default="25000")
    target_faces_custom: IntProperty(name="Custom faces", default=25000, min=100)

    # --- export (phase 11) ---
    export_target: EnumProperty(
        name="Export target", items=EXPORT_TARGET_ITEMS, default="SNAPPY")
    export_scale: FloatProperty(
        name="Scene → metres", default=1.0, min=1e-9,
        description="Scale factor applied on export (e.g. 0.001 for mm → m)")
    patch_naming: EnumProperty(
        name="Patch naming", items=PATCH_NAMING_ITEMS, default="PER_OBJECT")
    export_path: StringProperty(name="Output", subtype="FILE_PATH")

    # --- housekeeping ---
    keep_intermediates: BoolProperty(
        name="Keep intermediates", default=False,
        description="Leave temporary volume/mesh datablocks in the scene for debugging")

    def resolve_target_faces(self):
        """Effective target face count from the preset / custom fields."""
        if self.target_faces_preset == "CUSTOM":
            return self.target_faces_custom
        return int(self.target_faces_preset)


classes = (CADCFDProxyProperties,)


def register():
    for cls in classes:
        bpy.utils.register_class(cls)
    bpy.types.Scene.cad_cfd_proxy = PointerProperty(type=CADCFDProxyProperties)


def unregister():
    del bpy.types.Scene.cad_cfd_proxy
    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)
