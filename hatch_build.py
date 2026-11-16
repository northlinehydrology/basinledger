"""A wheel with MODFLOW 6 inside: set BASINLEDGER_LIBMF6 to USGS's libmf6 for this platform (``get-modflow``) and the
wheel carries it in ``basinledger/_lib/``, tagged for the platform (``tools/lib_tag.py``), so ``pip install
basinledger`` needs no ``get-modflow``. Without it the wheel is the pure one."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from hatchling.builders.hooks.plugin.interface import BuildHookInterface


class BundleLibmf6(BuildHookInterface):
    PLUGIN_NAME = "custom"

    def initialize(self, version, build_data):
        lib = os.environ.get("BASINLEDGER_LIBMF6")
        if self.target_name != "wheel" or not lib:
            return
        lib = Path(lib).resolve()
        if not lib.is_file():
            raise FileNotFoundError(f"BASINLEDGER_LIBMF6: {lib} isn't a file")
        r = subprocess.run([sys.executable, str(Path(self.root) / "tools" / "lib_tag.py"), str(lib)],
                           capture_output=True, text=True)
        if r.returncode:
            raise RuntimeError(f"{lib.name} can't go in a wheel as it is: {(r.stderr or r.stdout).strip()}")
        build_data["force_include"][str(lib)] = f"basinledger/_lib/{lib.name}"
        build_data["force_include"][str(Path(self.root) / "LICENSES" / "MODFLOW6.md")] = "basinledger/_lib/MODFLOW6.md"
        build_data["pure_python"] = False
        build_data["tag"] = f"py3-none-{r.stdout.strip()}"
