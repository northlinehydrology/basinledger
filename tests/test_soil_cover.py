"""The soil cover SWAT+'s soil evaporation sees (``[soil_cover]``, coupling/soilcover.py; docs/design.md, Soil cover
for soil evaporation):
- the config's section, checked on read, and the HRUs it names by name or land use;
- the wrapper, on a stand-in for swatplus-bmi;
- on SWAT+, with the synthetic example's HRU made a forest:
  - off, the surface model isn't wrapped;
  - SWAT+'s own cover set back is the same run to the bit;
  - leaves raise soil evaporation, and the basin still closes."""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest
from conftest import ROOT, RUNS, SYNTHETIC

from basinledger import config as C
from basinledger.balance import basin_balance
from basinledger.coupling import soilcover as SC
from basinledger.coupling.soilcover import Settings, selected
from basinledger.driver import Basin


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


# ---------- the wrapper, on a stand-in for swatplus-bmi ----------

class _Stub:
    """swatplus-bmi's cover variables on three HRUs, recording what is set."""

    def __init__(self):
        self.values = {SC.LEAF: np.array([0.3, 0.0, 0.2]), SC.ABOVE: np.array([9.0, 0.5, 6.0]),
                       SC.RESIDUE: np.array([1.0, 0.1, 0.4])}
        self.set, self.updates = [], 0

    def get_input_var_names(self):
        return (SC.COVER,)

    def get_output_var_names(self):
        return (SC.COVER, SC.LEAF, SC.ABOVE, SC.RESIDUE)

    def get_var_grid(self, var):
        return 0

    def get_grid_size(self, grid):
        return 3

    def get_var_type(self, var):
        return "float64"

    def get_value(self, var, dest):
        dest[:] = self.values[var]
        return dest

    def set_value_at_indices(self, var, inds, src):
        self.set.append((var, list(inds), list(src)))

    def update(self):
        self.updates += 1

    def get_current_time(self):
        return 7.0


def test_leaves_set_the_leaf_mass_and_residue_on_the_named_hrus_before_each_day():
    m = _Stub()
    w = SC.SoilCover(m, "leaves", np.array([0, 2]))
    w.update()
    assert m.set == [(SC.COVER, [0, 2], [pytest.approx(1.3), pytest.approx(0.6)])] and m.updates == 1
    assert w.days == 1 and w.get_current_time() == 7.0                  # everything else passes through


def test_swat_sets_back_its_own_cover():
    m = _Stub()
    SC.SoilCover(m, "swat", np.array([1])).update()
    assert m.set == [(SC.COVER, [1], [pytest.approx(0.6)])]


def test_the_sum_is_swatplus_own_4_byte_sum():
    """SWAT+ adds the masses in kg/ha as 4-byte reals. A float64 sum in kg m-2 rounds about 1 in 20 of them otherwise;
    the rule's sum never does."""
    rng = np.random.default_rng(0)
    n = 20000
    above = rng.uniform(1e3, 6e4, n).astype(np.float32)        # kg/ha, as SWAT+ holds them
    residue = rng.uniform(1e2, 2e4, n).astype(np.float32)
    m = _Stub()
    m.get_grid_size = lambda grid: n
    m.values = {SC.LEAF: np.zeros(n), SC.ABOVE: above.astype(np.float64) * 1e-4,
                SC.RESIDUE: residue.astype(np.float64) * 1e-4}
    cover = SC.SoilCover(m, "swat", np.arange(n)).cover()
    assert np.array_equal((cover * 1e4).astype(np.float32), above + residue)
    naive = ((m.values[SC.ABOVE] + m.values[SC.RESIDUE]) * 1e4).astype(np.float32)
    assert (naive != above + residue).sum() > n // 50                 # what it guards against


def test_a_swatplus_bmi_without_the_cover_is_refused():
    m = _Stub()
    m.get_output_var_names = lambda: (SC.LEAF,)
    m.get_input_var_names = lambda: ()
    with pytest.raises(ValueError, match="needs swatplus-bmi's soil cover"):
        SC.SoilCover(m, "leaves", np.array([0]))


# ---------- on SWAT+: Osu_1hru's HRU made a forest ----------

DAYS = 730


def forest(model, folder):
    """The synthetic example's HRU as a temperate forest, growing from the start (50 t/ha, LAI 2, 10 t/ha of
    residue, as Ames_sub1's frst), with no paddy and no management schedule: no decision-table actions."""
    if model != "swatplus":
        return
    p = folder / "hru-data.hru"
    p.write_text(p.read_text().replace("rice140_lum", "frst_tems_lum").replace("paddy0001", "null"))
    p = folder / "plant.ini"
    lines = p.read_text().splitlines()
    i = next(k for k, ln in enumerate(lines) if ln.split()[:1] == ["frst_tems"])
    lines[i] = (f"{'':40s}frst_tems{'':12s} y{2.0:14.5f}{50000.0:14.5f}{0.0:14.5f}{0.0:14.5f}{1.0:14.5f}"
                f"{10000.0:14.5f}")
    p.write_text("\n".join(lines) + "\n")


@pytest.fixture(scope="module")
def cover(swatplus):
    """SWAT+ through a swatplus-bmi with the soil cover."""
    import swatplus_bmi
    src = Path(swatplus_bmi.__file__).parent
    if not any(SC.COVER in f.read_text(errors="ignore") for f in src.rglob("*") if f.suffix in (".py", ".f90")):
        pytest.skip(f"the soil cover needs a swatplus-bmi with {SC.COVER} (this is {swatplus_bmi.__version__})")
    return swatplus


def _run(name, soil_cover=None, two_way=False, days=DAYS):
    cfg = C.load(ROOT / "examples" / "synthetic" / ("two_way.toml" if two_way else "basin.toml"))
    cfg.run_folder = RUNS
    cfg.soil_cover = soil_cover
    if two_way:
        cfg.outputs = dict(cfg.outputs, balance="cha01")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(cfg, run=name, days=days, prepare=forest)
    try:
        b.run()
        b.write_outputs()
        return b.ledger.table(), b.swat
    finally:
        b.finalize()


def _same(a: dict, b: dict) -> bool:
    return a.keys() == b.keys() and all(np.asarray(a[k]).tobytes() == np.asarray(b[k]).tobytes() for k in a)


def test_off_the_surface_model_is_not_wrapped(cover):
    t, swat = _run("cover-off", days=30)
    assert not isinstance(swat, SC.SoilCover)


def test_swatplus_own_cover_set_back_is_the_same_run_to_the_bit(cover):
    """On a forest with no decision-table actions the masses read before a day are those et_act uses, so setting
    above-ground + residue back changes nothing, every day."""
    off, _ = _run("cover-off")
    same, swat = _run("cover-swat", Settings("swat", land_uses=("frst_tems_lum",)))
    assert swat.days == DAYS
    assert _same(off, same)


def test_leaves_raise_the_forests_soil_evaporation_and_the_basin_still_closes(cover):
    _run("cover-off-2way", two_way=True)
    _, swat = _run("cover-leaves-2way", Settings("leaves", land_uses=("frst_tems_lum",)), two_way=True)
    assert swat.days == DAYS
    off = basin_balance(RUNS / "cover-off-2way", "cha01")
    on = basin_balance(RUNS / "cover-leaves-2way", "cha01")
    assert on["et_land"] > off["et_land"]                     # the wood no longer shades the soil
    p = on["precipitation_land"]
    assert abs(on["land_residual"]) <= 0.002 * (p + on["irrigation"])
    assert abs(on["inflow_from_elsewhere"]) <= 0.001 * on["outflow"]


def test_a_land_use_no_hru_has_is_refused(cover):
    with pytest.raises(ValueError, match="land_uses .*frst_lum"):
        _run("cover-bad", Settings("leaves", land_uses=("frst_lum",)), days=1)
