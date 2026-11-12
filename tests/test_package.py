"""The package as a user gets it: the wheel installs in a clean environment, its extras import, and the command
line runs an example from start to finish."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys

import pytest

from conftest import ROOT, RUNS


def test_extras_import():
    """The modules each extra exists for: users (waterledger-models), nextgen (netCDF4)."""
    pytest.importorskip("netCDF4")
    from waterledger_models.datacenter.bmi import CoolingWaterBmi  # users

    import basinledger
    from basinledger.coupling import nextgen  # noqa: F401
    from basinledger import forcing  # noqa: F401

    assert basinledger.__version__
    CoolingWaterBmi().initialize(None)


@pytest.mark.slow
def test_wheel_installs_in_a_clean_environment(tmp_path):
    uv = shutil.which("uv")
    if uv is None:
        pytest.skip("uv not found")
    subprocess.run([uv, "build", "--wheel", "--out-dir", str(tmp_path / "dist"), str(ROOT)], check=True,
                   capture_output=True)
    wheel = next((tmp_path / "dist").glob("basinledger-*.whl"))
    venv = tmp_path / "venv"
    subprocess.run([uv, "venv", "--python", f"{sys.version_info.major}.{sys.version_info.minor}", str(venv)],
                   check=True, capture_output=True)
    py = venv / "bin" / "python"
    subprocess.run([uv, "pip", "install", "--python", str(py), f"{wheel}[users,nextgen]"], check=True,
                   capture_output=True)
    out = subprocess.run([str(py), "-c", "import basinledger, basinledger.driver, basinledger.coupling.nextgen, "
                          "netCDF4, waterledger_models.datacenter.bmi; print(basinledger.__version__)"],
                         check=True, capture_output=True, text=True, cwd=tmp_path)
    import basinledger
    assert out.stdout.strip() == basinledger.__version__
    h = subprocess.run([str(venv / "bin" / "basinledger"), "--help"], check=True, capture_output=True, text=True)
    assert "run" in h.stdout and "status" in h.stdout


def test_command_line_end_to_end(tmp_path, swatplus):
    """basinledger run on the surface example: SWAT+, MODFLOW 6, a utility, two users and three rules, 30 days."""
    out = tmp_path / "out"
    r = subprocess.run([sys.executable, "-m", "basinledger", "run", str(ROOT / "examples" / "synthetic" / "surface.toml"),
                        "--days", "30", "--run", "cli-e2e", "--out", str(out)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "ledger closed every day" in r.stdout and "screening" in r.stdout
    s = json.loads((out / "summary.json").read_text())
    assert s["days"] == 30 and s["totals"]["route:well:pw1:delivered"] > 0
    assert (out / "ledger.csv").exists() and (out / "series.csv").exists()


def test_a_swatplus_built_by_another_wrapper_says_to_rebuild(swatplus, monkeypatch):
    """A SWAT+ built by an earlier swatplus-bmi speaks an older exchange protocol; the error says how to rebuild."""
    from swatplus_bmi.exchange import ExchangeError

    from basinledger import driver
    from conftest import synthetic_config

    def old_build(self, path):
        raise ExchangeError("SWAT+ speaks exchange protocol 4, this wrapper 5")

    monkeypatch.setattr(driver.SwatPlusBmi, "initialize", old_build)
    with pytest.raises(ExchangeError, match=r"protocol 4, this wrapper 5: .*swatplus-bmi build --jobs 2"):
        driver.Basin(synthetic_config(), run="old-build", days=2)
