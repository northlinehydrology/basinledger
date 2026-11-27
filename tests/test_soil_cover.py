"""The soil cover SWAT+'s soil evaporation sees (``[soil_cover]``, coupling/soilcover.py; docs/design.md, Soil cover
for soil evaporation): the config's section, checked on read, and the HRUs it names by name or land use."""

from __future__ import annotations

import pytest
from conftest import SYNTHETIC

from basinledger import config as C
from basinledger.coupling.soilcover import Settings, selected


def _load(tmp_path, section: str):
    (tmp_path / "basin.toml").write_text(SYNTHETIC.read_text() + "\n" + section)
    return C.load(tmp_path / "basin.toml")


def test_off_without_the_section():
    assert C.load(SYNTHETIC).soil_cover is None


def test_the_section_is_read(tmp_path):
    cfg = _load(tmp_path, '[soil_cover]\nrule = "leaves"\nland_uses = ["frst_lum", "frse_lum"]\nhrus = ["hru0003"]\n')
    assert cfg.soil_cover == Settings(rule="leaves", hrus=("hru0003",), land_uses=("frst_lum", "frse_lum"))


@pytest.mark.parametrize("section, match", [
    ('[soil_cover]\nland_uses = ["frst_lum"]\n', "rule is one of"),
    ('[soil_cover]\nrule = "wood"\nland_uses = ["frst_lum"]\n', "rule is one of"),
    ('[soil_cover]\nrule = "leaves"\n', "names no HRUs or land uses"),
    ('[soil_cover]\nrule = "leaves"\nland_uses = "frst_lum"\n', "list of names"),
    ('[soil_cover]\nrule = "leaves"\nhrus = ["hru1", "hru1"]\n', "names one twice"),
    ('[soil_cover]\nrule = "leaves"\nland_uses = ["frst_lum"]\nlai = 2\n', "unknown settings"),
])
def test_a_bad_section_is_refused(tmp_path, section, match):
    with pytest.raises(ValueError, match=match):
        _load(tmp_path, section)


def test_the_hrus_named_directly_or_by_land_use():
    hrus = ["hru1", "hru2", "hru3", "hru4"]
    landuse = ["frst_lum", "past_lum", "frst_lum", "corn_lum"]
    assert selected(Settings("leaves", land_uses=("frst_lum",)), hrus, landuse).tolist() == [0, 2]
    assert selected(Settings("leaves", hrus=("hru4",), land_uses=("frst_lum",)), hrus, landuse).tolist() == [0, 2, 3]
    assert selected(Settings("leaves", hrus=("hru2", "hru1")), hrus, landuse).tolist() == [0, 1]


def test_names_the_project_does_not_have_are_refused():
    with pytest.raises(ValueError, match="hrus .*hru9"):
        selected(Settings("leaves", hrus=("hru9",)), ["hru1"], ["frst_lum"])
    with pytest.raises(ValueError, match="land_uses .*frsd_lum"):
        selected(Settings("leaves", land_uses=("frsd_lum",)), ["hru1"], ["frst_lum"])
