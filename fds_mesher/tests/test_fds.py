import os
import shutil
import subprocess
import tempfile
import unittest

import fds_mesher as fm
from fds_mesher.tests import cases


@unittest.skipUnless(shutil.which("fds"), "fds not on PATH")
class FdsStarts(unittest.TestCase):
    def test_fds_accepts_generated_meshes(self):
        res = fm.run(cases.rock_tunnel(), 0.2, block=(10, 5, 5))
        txt = res.output.replace("&HEAD CHID='t'", "&HEAD CHID='mesher_t'") \
            + "&TIME T_END=0 /\n&TAIL /\n"
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "mesher_t.fds")
            with open(path, "w") as f:
                f.write(txt)
            p = subprocess.run(["fds", path], cwd=d, capture_output=True,
                               text=True, timeout=120)
            out = p.stdout + p.stderr
            self.assertNotIn("ERROR", out.upper().replace("STOP: FDS COMPLETED", ""))
