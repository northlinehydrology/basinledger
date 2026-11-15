"""Irrigation pumped from MODFLOW 6's wells (coupling/irrigation.py; docs/design.md, "Large users and irrigation,
MODFLOW first").

- each HRU's irrigation sources read from the project: Osu_1hru's paddy irrigates from an outside source (its
  ``pad_irr_r`` table's irrigate from null); the same table irrigating from a channel is refused for a mapped HRU,
  and an irrigated HRU left unmapped is noted as water from outside (tier 1);
- the map read, from the config or a CSV, its shares checked (tier 1);
- coupled on the synthetic existing model through Osu_1hru's first irrigation (its 128th day): the water SWAT+ put
  on the paddy each day is pumped that day from two added wells in their shares, MODFLOW 6 delivers it, and the
  ledger closes every day (tier 2).
"""

from __future__ import annotations

import shutil
import warnings

import pytest
from test_existing import CONFIG, SYNTHETIC, config

from basinledger import config as C
from basinledger.coupling import irrigation as IR

FIRST = 127  # Osu_1hru's first day with irrigation (0-based): 8 May 2010, 60 mm on the paddy


@pytest.fixture
def osu(swatplus, tmp_path):
    p = tmp_path / "Osu_1hru"
    shutil.copytree(swatplus / "Osu_1hru", p)
    return p


def test_the_paddy_irrigates_from_outside(osu):
    assert IR.sources(osu) == {"hru0001": {IR.OUTSIDE}}
    assert IR.problems(osu, {"hru0001": [("irr1", 1.0)]}, {"irr1"}) == ([], [])


def test_irrigation_from_a_channel_isnt_pumped(osu):
    dtl = osu / "lum.dtl"
    lines = dtl.read_text().splitlines()
    start = next(i for i, x in enumerate(lines) if x.split()[:1] == ["pad_irr_r"])
    end = next(i for i in range(start, len(lines)) if lines[i].split()[:2] == ["name", "conds"])
    for i in range(start, end):  # its two irrigate rows (on and off), from a channel instead
        w = lines[i].split()
        if w[:2] == ["irrigate", "null"]:
            lines[i] = "  ".join([w[0], "sdc"] + w[2:])
    dtl.write_text("\n".join(lines) + "\n")
    assert IR.sources(osu)["hru0001"] == {"a SWAT+ channel"}
    bad, _ = IR.problems(osu, {"hru0001": [("irr1", 1.0)]}, {"irr1"})
    assert len(bad) == 1 and "HRU hru0001 irrigates from a SWAT+ channel" in bad[0]


def test_what_the_check_refuses_and_notes(osu):
    bad, notes = IR.problems(osu, {"hru0002": [("irr1", 1.0)], "hru0001": [("irr9", 1.0)]}, {"irr1"})
    assert "[irrigation] HRU hru0002: no such HRU in the SWAT+ project (hru.con)" in bad
    assert "[irrigation] HRU hru0001: well 'irr9' isn't one of the [[wells]]" in bad
    _, notes = IR.problems(osu, {}, set())
    assert notes == [("irrigation on 1 HRU(s) not in [irrigation] (hru0001) stays SWAT+'s water from outside "
                      "the basin")]


def test_the_map_read(tmp_path):
    assert IR.read_map(tmp_path, [{"hru": "h1", "well": "a"}]) == {"h1": [("a", 1.0)]}
    (tmp_path / "m.csv").write_text("hru,well,share\nh1,a,0.25\nh1,b,0.75\nh2,a,1\n")
    assert IR.read_map(tmp_path, "m.csv") == {"h1": [("a", 0.25), ("b", 0.75)], "h2": [("a", 1.0)]}
    with pytest.raises(ValueError, match="HRU h1's shares add up to 0.5, not 1"):
        IR.read_map(tmp_path, [{"hru": "h1", "well": "a", "share": 0.5}])
    with pytest.raises(ValueError, match="HRU h1 names a well twice"):
        IR.read_map(tmp_path, [{"hru": "h1", "well": "a", "share": 0.5}, {"hru": "h1", "well": "a", "share": 0.5}])


def test_the_irrigation_settings_read(tmp_path):
    base = CONFIG.read_text().replace('"maps/', f'"{SYNTHETIC}/maps/')
    extra = ('\n[[wells]]\nname = "irr1"\npackage = "wel_irr"\ncell = [1, 6, 6]\nadd = true\n'
             '\n[irrigation]\nwells = [{hru = "hru0001", well = "irr1"}]\n')
    (tmp_path / "x.toml").write_text(base + extra)
    assert C.load(tmp_path / "x.toml").irrigation == {"hru0001": [("irr1", 1.0)]}
    (tmp_path / "x.toml").write_text(base + extra.replace('well = "irr1"', 'well = "irr2"'))
    with pytest.raises(ValueError, match=r"\[irrigation\] wells \['irr2'\] aren't in \[\[wells\]\]"):
        C.load(tmp_path / "x.toml")


def irrigated(cfg):
    """Two wells the model doesn't have, in its bottom (confined) layer, pumping the paddy's irrigation 60:40."""
    cfg.wells += [C.WellConfig("irr1", "wel_irr", (1, 6, 6), add=True),
                  C.WellConfig("irr2", "wel_irr", (1, 8, 4), add=True)]
    cfg.irrigation = {"hru0001": [("irr1", 0.6), ("irr2", 0.4)]}
    return cfg


def test_irrigation_pumped_from_wells(swatplus):
    from basinledger.driver import Basin
    days = FIRST + 8
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            b = Basin(irrigated(config()), run="existing-run-irrigation", days=days)
        except ValueError as e:
            if "swatplus-bmi 0.11.0" in str(e):
                pytest.skip("swatplus-bmi before 0.11.0: no HRU irrigation to read")
            raise
        try:
            b.run()
            t = b.ledger.table()
        finally:
            b.finalize()
    pumped = t["irrigation:pumped"]
    assert (pumped[:FIRST] == 0).all() and pumped[FIRST] == pytest.approx(60.0 / 1000 * 10 * 1e4, rel=1e-6)
    assert pumped.sum() > 6000.0
    assert (t["irrigation:from_outside"] == 0).all()
    for w, share in (("irr1", 0.6), ("irr2", 0.4)):
        assert t[f"irrigation:well:{w}"] == pytest.approx(share * pumped, rel=1e-12)
        assert t[f"route:well:{w}:delivered"] == pytest.approx(t[f"route:well:{w}:requested"], rel=1e-9)
        assert t[f"route:well:{w}:requested"] == pytest.approx(share * pumped, rel=1e-12)
    assert abs(t["groundwater:closure_share"]).max() < 1e-6
