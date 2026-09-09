# -*- coding: utf-8 -*-
"""Blender version and API capability layer (plan Phase 13).

Everything version- or API-sensitive is funnelled through here so the rest of
the addon can stay declarative. Two rules:

* **Never look sockets up by index** — identifiers and positions drift between
  Blender releases. Use :func:`find_socket`.
* **Never assume a node type exists** — probe it via :func:`check_capabilities`
  and degrade with a clear message rather than a traceback.
"""

import bpy

# Minimum supported Blender. NOTE: this is a conservative floor (extension
# system + Geometry Nodes Mesh/Volume conversion nodes present). Confirm the
# real target build and tighten if needed.
MIN_VERSION = (4, 2, 0)

# Node bl_idnames the pipeline depends on.
REQUIRED_NODES = (
    "GeometryNodeMeshToVolume",
    "GeometryNodeVolumeToMesh",
)

# Modifier types the pipeline depends on.
REQUIRED_MODIFIERS = (
    "NODES",
    "DECIMATE",
)


def get_version():
    """Return the running Blender version as a ``(major, minor, patch)`` tuple."""
    return bpy.app.version


def is_supported():
    """``(ok, detail)`` — whether the running build meets :data:`MIN_VERSION`."""
    version = get_version()
    if version < MIN_VERSION:
        return False, "Blender %s < required %s" % (
            ".".join(map(str, version)),
            ".".join(map(str, MIN_VERSION)),
        )
    return True, "Blender %s" % ".".join(map(str, version))


def _node_available(bl_idname):
    return hasattr(bpy.types, bl_idname)


def _modifier_available(mod_type):
    # Modifier types are enum items on the Modifier.type property.
    rna = bpy.types.Modifier.bl_rna.properties.get("type")
    if rna is None:
        return False
    return mod_type in {item.identifier for item in rna.enum_items}


def check_capabilities():
    """Return an ordered list of ``(label, ok, detail)`` capability probes."""
    ok, detail = is_supported()
    results = [("Blender version", ok, detail)]
    results.append(("Geometry Nodes", hasattr(bpy.types, "GeometryNodeTree"),
                    "GeometryNodeTree"))
    for node in REQUIRED_NODES:
        results.append((node, _node_available(node), "node type"))
    for mod in REQUIRED_MODIFIERS:
        results.append(("%s modifier" % mod, _modifier_available(mod), "modifier type"))
    return results


def run_diagnostics(log=False):
    """Run :func:`check_capabilities`; optionally print. Returns the results."""
    results = check_capabilities()
    if log:
        for label, ok, detail in results:
            print("[cad_cfd_proxy] %-24s %s  (%s)"
                  % (label, "OK" if ok else "MISSING", detail))
    return results


def all_capabilities_ok():
    return all(ok for _, ok, _ in check_capabilities())


def missing_capabilities():
    """Return the labels of capability probes that are failing (empty = all ok)."""
    return [label for label, ok, _ in check_capabilities() if not ok]


def find_socket(node, name, in_out="INPUT"):
    """Look a node socket up by *name*, not index.

    Returns the :class:`bpy.types.NodeSocket` or ``None``. Works for both nodes
    (``node.inputs`` / ``node.outputs``) and node-group interface items.
    """
    collection = node.inputs if in_out == "INPUT" else node.outputs
    for socket in collection:
        if socket.name == name:
            return socket
    return None
