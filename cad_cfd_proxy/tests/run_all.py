# -*- coding: utf-8 -*-
"""Run every headless test in one Blender process.

    blender --background --factory-startup --python-exit-code 1 \
        --python cad_cfd_proxy/tests/run_all.py

Each test module resets to factory settings in its own ``main()``, so they are
independent. Any uncaught exception propagates so ``--python-exit-code`` fails
the run.
"""

import os
import sys

# Make the repo root importable regardless of where Blender was launched from.
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from cad_cfd_proxy.tests import (
    test_smoke, test_collect, test_cleanup, test_volume, test_generate,
    test_decimate, test_export, test_validate, test_morphology,
)

for module in (test_smoke, test_collect, test_cleanup, test_volume,
               test_generate, test_decimate, test_export, test_validate,
               test_morphology):
    print("=== %s ===" % module.__name__)
    module.main()

print("all headless tests: PASS")
