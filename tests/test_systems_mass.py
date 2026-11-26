"""Constituents through a water system, coupled (docs/design.md, Mass through the human-system routes, parts 2 and 3).

On the synthetic salt example (chloride and weathering's ions through land, aquifer and streams), a town draws from
an intake on the channel first and then from a well over two cells. Its homes are sewered, its edge on septic over
three cells, half its losses leak to the ground, and a cooling user on its water returns to its sewer. Its plant's
chloride leaves at a monthly effluent concentration and its weathering ions pass through; it discharges at the
point source the example's plant also uses. Then on the nitrate example, the same town's plant loads SWAT+'s point
source with nitrate and soluble P.

- the sources' concentrations are the models' (the channel's mixed water, the well's cells by what each gave);
- every pair and each system's closure hold every day, with the salt balance's and the transport models';
- the septic water and leaks reach the aquifer the next day with their mass, the outfall's mass the stream;
- the cooling user takes in the town's water a day later;
- salt's septic land source and the routes' septic are refused together unless it says outside_routes.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from conftest import ROOT, RUNS

from basinledger import config as C
from basinledger.config import UserConfig
from basinledger.driver import Basin

pytestmark = pytest.mark.slow

DAY = 86400.0
DAYS = 120
WELL = (((0, 5, 6), 0.75), ((0, 7, 6), 0.25))
EDGE = (((0, 2, 7), 0.25), ((0, 3, 7), 0.5), ((0, 4, 7), 0.25))
MAINS = (((0, 6, 3), 1.0),)
EFFLUENT_CL = [0.10] * 6 + [0.12] * 6   # kg m-3 by month, as a plant's DMR gives it [A]


def _town(effluent, septic):
    return {"name": "town", "outfall": "swat:point_source:wwtp", "losses": 0.1,
            "leaks": {"share": 0.5, "to": "ground:mains"},
            "sources": [{"route": "swat:place_of_use:town_intake", "priority": 1, "capacity": 0.0008},
                        {"route": "well:pw2", "priority": 2}],
            "customers": [{"name": "homes", "demand": 0.0008, "sewered": 0.7, "consumed": 0.3},
                          {"name": "edge", "demand": 0.0003, "septic": 0.9, "septic_to": "ground:edge",
                           "consumed": 0.1}],
            "effluent": effluent, "septic_effluent": septic}


def _config(example, effluent, septic, outside_routes=True):
    cfg = C.load(ROOT / "examples" / "synthetic" / example)
    cfg.run_folder = RUNS
    cfg.places_of_use = [("town_intake", "cha01")]
    cfg.point_sources = [("wwtp", "cha01")]
    cfg.wells = cfg.wells + [C.WellConfig("pw2", "wel_town", None, True, cells=WELL)]
    cfg.ground = [C.WellConfig("edge", "wel_edge", None, True, cells=EDGE),
                  C.WellConfig("mains", "wel_mains", None, True, cells=MAINS)]
    cfg.utilities = [_town(effluent, septic)]
    cool = UserConfig(name="cooling", component="basinledger.testing:SyntheticUser",
                      config=ROOT / "examples" / "synthetic" / "intake_user.toml",
                      withdrawal="utility:town", returns="utility:town")
    cfg.users = [u for u in cfg.users] + [cool]
    for s in (cfg.salt.get("chloride").sources if "chloride" in cfg.salt else []):
        if s.group == "septic":
            s.outside_routes = outside_routes
    cfg.outputs = {"channels": ["cha01"]}
    return cfg


def _run(cfg, name, days, each_day=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(cfg, run=name, days=days)
    rows = []
    for _ in range(days):
        b.step()
        if each_day is not None:
            rows.append(each_day(b))
    t = b.ledger.table()
    b.finalize()
    return t, rows


@pytest.fixture(scope="module")
def salted(swatplus):
    cfg = _config("salt.toml", {"chloride": EFFLUENT_CL, "weathering": "sewer"},
                  {"chloride": 0.08, "weathering": "tap"})
    return _run(cfg, "systems-mass-salt", DAYS)[0]


def test_every_pair_and_closure_holds(salted):
    t = salted
    for sp in ("chloride", "weathering"):
        assert np.nanmax(np.abs(t[f"utility:town:{sp}:closure_share"])) < 1e-12
        assert np.nanmax(np.abs(t[f"salt:{sp}:closure_share"])) < 1e-9
        assert np.nanmax(np.abs(t[f"transport:{sp}:closure_share"])) < 1e-6
        assert np.nansum(t[f"utility:town:{sp}:produced:kg"]) > 0
    # the ledger's pairs held every day, or the run would have stopped; the routes delivered
    assert np.nansum(t["utility:town:source:swat:place_of_use:town_intake:delivered"]) > 0
    assert np.nansum(t["utility:town:source:well:pw2:delivered"]) > 0


def test_the_sources_concentrations_are_the_models(salted):
    t = salted
    k = "utility:town:chloride:"
    # the intake's water left the channel at salt's mixed concentration: the town is its only taker
    np.testing.assert_allclose(t[k + "source:swat:place_of_use:town_intake:kg"], t["salt:chloride:streams:withdrawn"],
                               rtol=1e-6, atol=1e-12)
    # the well's, as the transport model took it out of its two cells
    np.testing.assert_allclose(t[k + "source:well:pw2:kg"], t["transport:chloride:well:pw2:withdrawn"], rtol=1e-9,
                               atol=1e-12)


def test_septic_and_leaks_reach_the_aquifer_the_next_day(salted):
    t = salted
    for sp in ("chloride", "weathering"):
        k = f"utility:town:{sp}:"
        for route in ("edge", "mains"):
            # day 0 releases nothing: nothing was held
            np.testing.assert_allclose(t[f"transport:{sp}:ground:{route}:due"],
                                       np.nan_to_num(t[k + f"to_ground:ground:{route}:kg"]), rtol=1e-12, atol=1e-15)
            np.testing.assert_allclose(t[f"transport:{sp}:ground:{route}:received"],
                                       t[f"transport:{sp}:ground:{route}:sent"], rtol=1e-9, atol=1e-12)
        np.testing.assert_allclose(t[k + "to_ground:ground:edge:kg"][1:], t[k + "septic:edge:kg"][:-1], rtol=1e-12)
        np.testing.assert_allclose(t[f"salt:{sp}:aquifer:from_ground"],
                                   np.nan_to_num(t[k + "to_ground:ground:edge:kg"] + t[k + "to_ground:ground:mains:kg"]),
                                   rtol=1e-9, atol=1e-12)
    # septic chloride at 0.08; weathering's at the tap water's ("tap"), so nothing added
    k = "utility:town:"
    np.testing.assert_allclose(t[k + "chloride:septic:edge:kg"], 0.08 * t[k + "septic:edge"], rtol=1e-12)
    np.testing.assert_allclose(t[k + "weathering:septic_added:edge"], 0.0, atol=1e-12)


def test_the_plants_mass_reaches_the_stream(salted):
    t = salted
    k = "utility:town:"
    # chloride at the month's effluent concentration (Jan to Apr here); weathering's ions pass the plant
    np.testing.assert_allclose(t[k + "chloride:outfall:kg"], 0.10 * t[k + "outfall"], rtol=1e-12)
    np.testing.assert_allclose(t[k + "weathering:plant:added"], 0.0, atol=1e-12)
    # the point source carries the plant user's reported chloride and the town's outfall, and for weathering's ions,
    # which the plant user doesn't report, its water at the stated 0.15
    np.testing.assert_allclose(t["salt:chloride:streams:point_sources"],
                               t["user:plant:chloride:kg"] + t[k + "chloride:outfall:kg"], rtol=1e-6, atol=1e-12)
    np.testing.assert_allclose(t["salt:weathering:streams:point_sources"],
                               0.15 * t["user:plant:return"] + t[k + "weathering:outfall:kg"], rtol=1e-5, atol=1e-12)


def test_the_cooling_user_takes_in_the_towns_water_a_day_later(salted):
    t = salted
    np.testing.assert_allclose(t["user:cooling:intake:chloride"][1:], t["utility:town:chloride:concentration"][:-1],
                               rtol=1e-12)
    # and returns it to the sewer, reported (SyntheticUser returns what it takes in)
    assert np.nansum(t["utility:town:chloride:users_return:kg"]) > 0


def test_septic_counted_once(swatplus):
    cfg = _config("salt.toml", {"chloride": 0.1, "weathering": "sewer"}, {"chloride": 0.08, "weathering": "tap"},
                  outside_routes=False)
    with pytest.raises(ValueError, match="would count the same septic chloride twice"):
        Basin(cfg, run="systems-mass-refused", days=2)


def test_missing_concentrations_are_refused(swatplus):
    cfg = _config("salt.toml", {"chloride": 0.1}, {"weathering": "tap"})
    with pytest.raises(ValueError, match=r"utility town: septic_effluent chloride \(class edge\)\n- utility town: effluent weathering"):
        Basin(cfg, run="systems-mass-refused", days=2)


def test_the_plant_loads_swatplus_point_source_with_nitrate(swatplus):
    """On the nitrate example: SWAT+ carries nitrate and soluble P in its streams, so the outfall's mass is set as
    the point source's load before SWAT+'s day; with the plant's effluent the channel carries more nitrate."""
    def conc(b):
        return float(b.sw.get("channel_water_no3__mass_concentration")[0])

    days = 40
    runs = {}
    for name, no3 in (("with", 0.02), ("without", 0.0)):
        cfg = _config("nitrate.toml", {"nitrate": no3, "phosphorus": 0.001},
                      {"nitrate": 0.04, "phosphorus": 0.0})
        runs[name] = _run(cfg, f"systems-mass-no3-{name}", days, conc)
    t, c = runs["with"]
    _, c0 = runs["without"]
    k = "utility:town:nitrate:"
    np.testing.assert_allclose(t["transport:nitrate:point_source:wwtp:sent"], t[k + "outfall:kg"], rtol=1e-12)
    np.testing.assert_allclose(t[k + "outfall:kg"], 0.02 * t["utility:town:outfall"], rtol=1e-12)
    assert np.nansum(t[k + "outfall:kg"]) > 0 and np.nansum(t["transport:nitrate:point_source:wwtp:unstated_m3"]) == 0
    c, c0 = np.array(c), np.array(c0)  # the channel is dry on some days, and carries nothing either way
    assert (c >= c0 - 1e-12).all() and c[1:].sum() > 2 * c0[1:].sum()
    for sp in ("nitrate", "phosphorus"):
        assert np.nanmax(np.abs(t[f"utility:town:{sp}:closure_share"])) < 1e-12
        assert np.nanmax(np.abs(t[f"transport:{sp}:closure_share"])) < 1e-6
        np.testing.assert_allclose(t[f"transport:{sp}:ground:edge:received"], t[f"transport:{sp}:ground:edge:sent"],
                                   rtol=1e-9, atol=1e-12)


def test_a_restart_carries_the_systems_mass(swatplus):
    n, at = 30, 12

    def run(name, save_at=None, restart=None):
        cfg = _config("salt.toml", {"chloride": EFFLUENT_CL, "weathering": "sewer"},
                      {"chloride": 0.08, "weathering": "tap"})
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            b = Basin(cfg, run=name, days=n, restart=restart)
        if save_at:
            b.run(save_at - b.day)
            b.save(b.folder / "state")
        else:
            b.run()
        t = b.ledger.table()
        b.finalize()
        return t

    a = run("systems-mass-unbroken")
    run("systems-mass-warm", save_at=at)
    b = run("systems-mass-restarted", restart=RUNS / "systems-mass-warm" / "state")
    for key in a:
        if key.startswith(("utility:town:chloride:", "utility:town:weathering:", "salt:chloride:", "user:cooling:")):
            np.testing.assert_array_equal(a[key], b[key], err_msg=key)
