"""Build the wheel with MODFLOW 6 inside and check it the way a user gets it: installed with nothing else into a new
environment, no get-modflow, no source checkout, and run.

    uv run python tools/check_wheel.py ~/.local/share/flopy/bin/libmf6.dylib [--out dist]

The wheel is built with BASINLEDGER_LIBMF6 (hatch_build.py), tagged for this platform (tools/lib_tag.py), and
installed by uv into a new virtual environment. There, with a home folder of its own (so no libmf6 that get-modflow
put there is found) and from a folder outside the checkout, Basin Ledger must find the wheel's own libmf6, and a
small MODFLOW 6 model (flopy, 10 by 10 cells, a well and a constant head) is run through it for 3 days. Prints the
wheel's path.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

RUN = """
import sys
from pathlib import Path
import flopy
from modflowapi import ModflowApi
from basinledger.coupling.modflow6 import find_library
lib = find_library()
assert "_lib" in lib.parts, f"not the wheel's own libmf6: {lib}"
ws = Path(sys.argv[1])
sim = flopy.mf6.MFSimulation(sim_ws=str(ws))
flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(3.0, 3, 1.0)])
flopy.mf6.ModflowIms(sim)
gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf")
flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=10, ncol=10, delr=100.0, delc=100.0, top=10.0, botm=0.0)
flopy.mf6.ModflowGwfic(gwf, strt=9.0)
flopy.mf6.ModflowGwfnpf(gwf, k=5.0)
flopy.mf6.ModflowGwfsto(gwf, sy=0.2, ss=1e-5, iconvert=1, transient={0: True})
flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, r, 0), 9.0] for r in range(10)])
flopy.mf6.ModflowGwfwel(gwf, stress_period_data=[[(0, 5, 7), -50.0]])
sim.write_simulation(silent=True)
m = ModflowApi(str(lib), working_directory=str(ws))
m.initialize()
while m.get_current_time() < m.get_end_time():
    m.update()
head = m.get_value(m.get_var_address("X", "GWF"))
m.finalize()
assert head[5 * 10 + 7] < 9.0, "the well drew no water down"
print(f"ran 3 days with {lib.name}; head at the well {head[5 * 10 + 7]:.3f}")
"""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("lib", type=Path)
    ap.add_argument("--out", type=Path, default=ROOT / "dist")
    a = ap.parse_args(argv)
    uv = shutil.which("uv")
    if uv is None:
        raise SystemExit("uv not found")

    a.out.mkdir(parents=True, exist_ok=True)
    before = set(a.out.glob("*.whl"))
    env = dict(os.environ, BASINLEDGER_LIBMF6=str(a.lib.resolve()))
    subprocess.run([uv, "build", "--wheel", "-o", str(a.out), str(ROOT)], check=True, env=env)
    wheel, = set(a.out.glob("*.whl")) - before
    print(f"built {wheel.name}")

    with tempfile.TemporaryDirectory(prefix="basinledger-wheel-") as tmp:
        tmp = Path(tmp)
        venv = tmp / "venv"
        subprocess.run([uv, "venv", "-q", "--python", f"{sys.version_info.major}.{sys.version_info.minor}", str(venv)],
                       check=True)
        py = venv / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")
        subprocess.run([uv, "pip", "install", "-q", "--python", str(py), str(wheel)], check=True)
        home = tmp / "home"
        home.mkdir()
        clean = {k: v for k, v in os.environ.items() if not k.startswith(("LIBMF6", "VIRTUAL_ENV", "PYTHON", "UV_"))}
        clean.update(HOME=str(home), USERPROFILE=str(home), LOCALAPPDATA=str(home / "AppData" / "Local"))
        (tmp / "model").mkdir()
        r = subprocess.run([str(py), "-I", "-c", RUN, str(tmp / "model")], cwd=tmp, env=clean,
                           capture_output=True, text=True)
        print(r.stdout, end="")
        if r.returncode:
            print(r.stderr, file=sys.stderr)
            return r.returncode
    print(wheel)
    return 0


if __name__ == "__main__":
    sys.exit(main())
