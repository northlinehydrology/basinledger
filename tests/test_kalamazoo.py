"""The Kalamazoo River basin model (Michigan; USGS doi:10.5066/F76D5R5V; examples/kalamazoo), MODFLOW-NWT converted to
MODFLOW 6 by mf5to6: the second test model, as given. Skipped until fetched and converted (model.zip is downloaded by
hand, behind ScienceBase's captcha; ``uv run python examples/kalamazoo/fetch.py``).

What these check is what Basin Ledger would meet in it, read from the converted files' headers, not by loading its
7.2 million cells: feet, a steady single period with no storage, its own boundaries, its multi-node wells as MAW, and
SFR with a Manning's constant for feet and days, which the exchange's check takes once the config states feet
(Basin Ledger converts at each link; coupling/modflow6.py, Units). Tier 0: no model runs.
"""

from __future__ import annotations

import importlib.util
import re

import pytest

from conftest import ROOT

spec = importlib.util.spec_from_file_location("kalamazoo_fetch", ROOT / "examples" / "kalamazoo" / "fetch.py")
F = importlib.util.module_from_spec(spec)
spec.loader.exec_module(F)
WS = F.mf6_folder()
pytestmark = pytest.mark.skipif(not (WS / "mfsim.nam").exists(),
                                reason="Kalamazoo not fetched: see examples/kalamazoo/fetch.py")


def head(name: str, lines: int = 400) -> str:
    with open(WS / name) as f:
        return "".join(next(f, "") for _ in range(lines))


def test_packages_as_converted():
    nam = (WS / "kala.nam").read_text()
    kinds = re.findall(r"^\s*(\w+)6\s", nam.split("BEGIN Packages")[1], flags=re.M)
    assert kinds == ["DIS", "IC", "NPF", "OC", "MAW", "SFR", "DRN", "RIV", "GHB", "RCH", "CHD"]
    assert "STO" not in kinds  # steady: no storage, so daily steps need specific yield and storage added


def test_feet_and_days():
    assert "LENGTH_UNITS  UNDEFINED" in head("kala.dis", 10)
    assert "length_units feet" in (F.nwt_folder() / "usgs.model.reference").read_text()
    assert "TIME_UNITS DAYS" in head("kala.tdis")
    assert "NPER 1" in head("kala.tdis")
    # Manning's constant for feet and days: 1.486 x 86,400, which the exchange's check takes for feet (1.4859 to five
    # figures; a thousandth allowed), once [modflow6] length_units = "feet" states what the DIS leaves undefined
    unitconv = float(re.search(r"UNIT_CONVERSION\s+([\d.]+)", head("kala.sfr")).group(1))
    assert unitconv == pytest.approx(1.486 * 86400)
    from basinledger.coupling.modflow6 import DAYS, METERS
    assert abs(unitconv / (86400.0 * DAYS["days"] / METERS["feet"] ** (1 / 3)) - 1.0) < 1e-3


def test_its_size_and_its_boundaries():
    dims = head("kala.dis", 12)
    assert [int(re.search(rf"{k}\s+(\d+)", dims).group(1)) for k in ("NLAY", "NROW", "NCOL")] == [8, 730, 1230]
    count = lambda name, key: int(re.search(rf"{key}\s+(\d+)", head(name, 40)).group(1))  # noqa: E731
    assert count("kala.sfr", "NREACHES") == 20010
    assert count("kala.maw", "NMAWWELLS") == 710
    assert count("kala.riv6", "MAXBOUND") == 52969
    assert count("kala.drn", "MAXBOUND") == 29897
    assert count("kala.ghb", "MAXBOUND") == 137100
    assert count("kala.ib.ibound.chd", "MAXBOUND") == 28424


def test_storage_added_to_the_steady_model(tmp_path):
    """Kalamazoo has no storage package (steady): ``[modflow6] storage`` adds one to a run's copy, its SS written per
    foot, the units the study states (its DIS leaves them undefined). Text only: the name file copied, the grid
    linked, nothing run."""
    import shutil

    from basinledger.coupling.existing import add_storage, steady_periods
    for name in ("mfsim.nam", "kala.nam"):
        shutil.copy(WS / name, tmp_path / name)
    (tmp_path / "kala.dis").symlink_to(WS / "kala.dis")
    assert steady_periods(tmp_path, "kala") == "no STO"
    with pytest.raises(ValueError, match="length units are undefined"):
        add_storage(tmp_path, "kala", 0.15, 1e-5)
    add_storage(tmp_path, "kala", 0.15, 1e-5, length_units="feet")
    assert "STO6  kala.sto_basinledger.sto  sto_basinledger" in (tmp_path / "kala.nam").read_text()
    assert "CONSTANT 3.048e-06" in (tmp_path / "kala.sto_basinledger.sto").read_text()
    assert steady_periods(tmp_path, "kala") is None
