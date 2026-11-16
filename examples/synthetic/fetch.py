"""Get what the synthetic example needs. It downloads no data of its own.

- SWAT+ with the exchange mode, and SWAT+'s example projects (refdata/Osu_1hru): `swatplus-bmi build`.
- MODFLOW 6's shared library (libmf6): `get-modflow`, from flopy, takes it from USGS's MODFLOW 6 release.

    uv run python examples/synthetic/fetch.py          # checks, and says what's missing
    uv run python examples/synthetic/fetch.py --get    # also gets what's missing (SWAT+ is compiled with 2 jobs)
"""

from __future__ import annotations

import subprocess
import sys

from basinledger.config import swatplus_refdata
from basinledger.coupling.modflow6 import LIBNAME, MF6, find_library


def main(get: bool) -> int:
    missing = 0
    try:
        from swatplus_bmi.bmi import find_executable
        print("SWAT+ (exchange mode):", find_executable())
        ok = (swatplus_refdata() / "Osu_1hru" / "file.cio").exists()
    except FileNotFoundError:
        ok = False
    if not ok:
        missing += 1
        print("SWAT+ or its refdata missing: uv run swatplus-bmi build --jobs 2")
        if get:
            subprocess.run([sys.executable, "-m", "swatplus_bmi.build", "build", "--jobs", "2"], check=True)
    else:
        print("Osu_1hru:", swatplus_refdata() / "Osu_1hru")
    try:
        print("libmf6:", find_library())
    except FileNotFoundError:
        missing += 1
        print(f"libmf6 missing: uv run get-modflow :flopy --repo modflow6 --subset {MF6},{LIBNAME}")
        if get:
            subprocess.run(["get-modflow", ":flopy", "--repo", "modflow6", "--subset", f"{MF6},{LIBNAME}"], check=True)
    return 0 if (missing == 0 or get) else 1


if __name__ == "__main__":
    sys.exit(main("--get" in sys.argv))
