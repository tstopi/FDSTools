# -*- coding: utf-8 -*-
"""Headless registration smoke test (plan Phase 12 harness).

Run with Blender::

    blender --background --factory-startup \
        --python cad_cfd_proxy/tests/test_smoke.py

Exits non-zero on failure so it can gate CI. It verifies the addon
registers/unregisters cleanly and that the compat probes run — it does *not*
exercise the (still stubbed) pipeline.
"""

import sys


def main():
    import bpy  # noqa: F401 (only importable inside Blender)
    import cad_cfd_proxy

    cad_cfd_proxy.register()
    assert hasattr(bpy.types.Scene, "cad_cfd_proxy"), "properties not attached"

    results = cad_cfd_proxy.compat.check_capabilities()
    assert results, "no capability probes ran"
    for label, ok, detail in results:
        print("  %-24s %s (%s)" % (label, "OK" if ok else "MISSING", detail))

    # Operators registered?
    for op in ("cadcfd.generate", "cadcfd.preview", "cadcfd.export"):
        module, name = op.split(".")
        assert hasattr(getattr(bpy.ops, module), name), "missing operator %s" % op

    cad_cfd_proxy.unregister()
    assert not hasattr(bpy.types.Scene, "cad_cfd_proxy"), "properties not detached"
    print("smoke test: PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        print("smoke test: FAIL: %s" % exc)
        sys.exit(1)
