"""Mass on the returns to the ground and the wells over several cells (docs/design.md, Mass through the human-system
routes, part 1).

On the synthetic transport example (chloride at stated concentrations, coupled both ways), a house pumps a well
over two cells (3:1) and returns 80% of it to a drainfield over three (1:2:1), at a stated septic concentration
(its ReportedUser reports no chloride); the example's plant takes in the aquifer's chloride from the same well:

- each drainfield cell takes its weight's share of the route's mass the same day as its water, and the transport
  model receives what was sent, every day;
- the well's mass is its cells' (each at its own concentration, by what it gave), and the plant's intake is
  yesterday's concentration of that water;
- the transport model's mass closes every day with both, and neither is booked as a boundary;
- with salt (chloride and weathering's ions), the whole system closes with the returns' mass from outside.
"""

from __future__ import annotations

import datetime as dt
import warnings

import numpy as np
import pytest
from conftest import ROOT, RUNS

from basinledger import config as C
from basinledger.config import UserConfig
from basinledger.coupling import swatplus as S
from basinledger.driver import Basin

pytestmark = pytest.mark.slow

DAY = 86400.0
DAYS = 60
PUMP = 0.0004             # m3 s-1, the house
RETURNED = 0.8            # to the drainfield; the rest consumed
SEPTIC = 0.09             # kg m-3 of chloride reaching the water table under the drainfield [A]
WELL = (((0, 5, 6), 0.75), ((0, 7, 6), 0.25))
FIELD = (((0, 2, 7), 0.25), ((0, 3, 7), 0.5), ((0, 4, 7), 0.25))


def _house(cfg, tmp_path, days):
    start = S._period(cfg.swat_project)[0]
    w = PUMP * DAY
    rows = ["date,withdrawal_m3,consumption_m3,return_m3"]
    rows += [f"{start + dt.timedelta(days=d)},{w},{(1 - RETURNED) * w},{RETURNED * w}" for d in range(days)]
    (tmp_path / "house.csv").write_text("\n".join(rows) + "\n")
    (tmp_path / "house.toml").write_text(f'series = "house.csv"\nstart = "{start}"\n')
    return UserConfig(name="house", component="basinledger.generic:ReportedUser", config=tmp_path / "house.toml",
                      withdrawal="well:pw2", returns="ground:field")


def _config(example, tmp_path, days, concentration):
    cfg = C.load(ROOT / "examples" / "synthetic" / example)
    cfg.run_folder = RUNS
    cfg.wells = [C.WellConfig("pw2", "wel_house", None, True, cells=WELL)]
    cfg.ground = [C.WellConfig("field", "wel_ground", None, True, cells=FIELD, concentration=concentration)]
    plant = next(u for u in cfg.users if u.name == "plant")
    plant.withdrawal = "well:pw2"  # the plant takes in the aquifer's chloride from the well over two cells
    cfg.users = [plant, _house(cfg, tmp_path, days)]
    cfg.outputs = {}
    return cfg


@pytest.fixture(scope="module")
def transport(swatplus, tmp_path_factory):
    cfg = _config("transport.toml", tmp_path_factory.mktemp("ground-mass"), DAYS, {"chloride": SEPTIC})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(cfg, run="ground-mass", days=DAYS)
    t = b.transport["chloride"]
    cells = {}  # each entry's mass and water as the transport model took them, every day
    for _ in range(DAYS):
        b.step()
        for name, (pkg, i) in list(b._ground_links.items()) + list(b._mass_links[1].items()):
            q, m = t.last[pkg.upper()]
            cells.setdefault(name, []).append((float(q[i]), float(m[i]), float(t.concentration()[t.entry_nodes(pkg)[i]])))
    out = b.ledger.table(), cells
    b.finalize()
    return out


def test_each_drainfield_cell_takes_its_weights_share_the_same_day(transport):
    t, cells = transport
    k = "transport:chloride:ground:field:"
    due = SEPTIC * RETURNED * PUMP * DAY
    np.testing.assert_allclose(t[k + "due"], due, rtol=1e-12)
    np.testing.assert_allclose(t[k + "sent"], t[k + "due"], rtol=1e-12)
    np.testing.assert_allclose(t[k + "no_water"], 0.0, atol=0.0)
    np.testing.assert_allclose(t[k + "received"], t[k + "sent"], rtol=1e-9)
    for i, (_, w) in enumerate(FIELD, 1):
        q, m = np.array([r[0] for r in cells[f"field#{i}"]]), np.array([r[1] for r in cells[f"field#{i}"]])
        np.testing.assert_allclose(q, w * RETURNED * PUMP * DAY, rtol=1e-9)   # its water
        np.testing.assert_allclose(m, w * due, rtol=1e-9)                      # and its mass, the same day


def test_the_well_takes_each_cells_concentration_by_what_it_gave(transport):
    t, cells = transport
    k = "transport:chloride:well:pw2"
    parts = [cells[f"pw2#{i}"] for i in (1, 2)]
    water = -sum(np.array([r[0] for r in p]) for p in parts)
    for i, p in enumerate(parts, 1):  # each cell's mass leaves at its own concentration
        np.testing.assert_allclose(t[f"{k}#{i}:withdrawn"], [-r[0] * r[2] for r in p], rtol=1e-9)
    np.testing.assert_allclose(t[k + ":withdrawn"], t[f"{k}#1:withdrawn"] + t[f"{k}#2:withdrawn"], rtol=1e-12)
    np.testing.assert_allclose(t[k + ":concentration"], t[k + ":withdrawn"] / water, rtol=1e-12)
    np.testing.assert_allclose(water, t["route:well:pw2:delivered"], rtol=1e-9)
    # the cells differ, so the weighting matters
    c1, c2 = np.array([r[2] for r in parts[0]]), np.array([r[2] for r in parts[1]])
    assert np.abs(c1 - c2).max() > 1e-6 * c1.max()
    # the plant takes in yesterday's water from the well
    np.testing.assert_allclose(t["user:plant:intake:chloride"][1:], t[k + ":concentration"][:-1], rtol=1e-12)


def test_the_transport_model_closes_with_them(transport):
    t, _ = transport
    assert np.nanmax(np.abs(t["transport:chloride:closure_share"])) < 1e-6
    for p in ("wel_ground", "wel_house"):
        assert f"transport:chloride:boundary:{p}" not in t
    # the drainfield raises the aquifer's chloride under it: 0.09 kg m-3 into water starting at 0.01
    _, cells = transport
    assert cells["field#2"][-1][2] > 0.011


def test_the_ssm_takes_the_ground_package(transport, swatplus):
    text = (RUNS / "ground-mass" / "mf6" / "gwt_chloride.ssm").read_text().lower()
    assert "wel_ground  aux  chloride" in text
    assert "auxiliary chloride" in (RUNS / "ground-mass" / "mf6" / "gwf.wel_ground.wel").read_text().lower()


def test_a_return_with_no_mass_for_a_species_is_refused(swatplus, tmp_path):
    cfg = _config("transport.toml", tmp_path, 3, {})
    with pytest.raises(ValueError, match=r"returns to ground:field reporting no \['chloride'\]"):
        Basin(cfg, run="ground-mass-refused", days=3)
    cfg = _config("transport.toml", tmp_path, 3, {"nitrate": 0.01, "chloride": SEPTIC})
    with pytest.raises(ValueError, match=r"concentration for \['nitrate'\], which aren't transported"):
        Basin(cfg, run="ground-mass-refused", days=3)


def test_with_salt_the_whole_system_closes(swatplus, tmp_path):
    cfg = _config("salt.toml", tmp_path, 90, {"chloride": SEPTIC, "weathering": 0.2})
    for s in cfg.salt["chloride"].sources:  # the example's septic land source covers homes the drainfield isn't
        if s.group == "septic":
            s.outside_routes = True
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(cfg, run="ground-mass-salt", days=90)
    b.run()
    t = b.ledger.table()
    b.finalize()
    for sp, c in (("chloride", SEPTIC), ("weathering", 0.2)):
        k = f"salt:{sp}:"
        assert np.nanmax(np.abs(t[k + "closure_share"])) < 1e-9
        assert np.nanmax(np.abs(t[f"transport:{sp}:closure_share"])) < 1e-6
        np.testing.assert_allclose(t[k + "aquifer:from_ground"], c * RETURNED * PUMP * DAY, rtol=1e-9)
        np.testing.assert_allclose(t[k + "aquifer:from_ground"], t[f"transport:{sp}:ground:field:received"],
                                   rtol=1e-12)
