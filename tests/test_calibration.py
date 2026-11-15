"""Calibration: tables, statistics, simulated values at observations, the forward model, parameters applied to
SWAT+ and MODFLOW 6, and PEST++ IES files."""

from __future__ import annotations

import csv

import numpy as np
import pytest

from basinledger import calibration as K
from conftest import ROOT, RUNS


def test_statistics_on_known_values():
    obs = np.array([1.0, 2.0, 3.0, 4.0])
    assert K.nse(obs, obs) == 1.0
    assert K.nse(np.full(4, obs.mean()), obs) == pytest.approx(0.0)
    assert K.kge(obs, obs) == pytest.approx(1.0)
    assert K.pbias(obs * 1.1, obs) == pytest.approx(10.0)
    assert K.rmse(obs + 0.5, obs) == pytest.approx(0.5) and K.mean_error(obs - 0.25, obs) == pytest.approx(-0.25)
    assert K.log_nse(obs, obs) == 1.0


def test_tables(tmp_path):
    p = tmp_path / "parameters.csv"
    p.write_text("name,model,target,value,lower,upper,transform\nk,modflow6,npf.k:all,5,0.5,50,log\n")
    assert K.read_parameters(p)[0].transform == "log"
    p.write_text("name,model,target,value,lower,upper,transform\nk,modflow6,npf.k:all,500,0.5,50,log\n")
    with pytest.raises(ValueError, match="outside"):
        K.read_parameters(p)
    o = tmp_path / "observations.csv"
    o.write_text("name,kind,object,date,value,weight,source\nq1,flow,cha01,2010-02-01,0.002,1,test\n"
                 "h1,head,pw1,2010-02-01,88.5,2,test\n")
    obs = K.read_observations(o)
    assert [x.kind for x in obs] == ["flow", "head"] and obs[1].weight == 2.0


def test_simulated_and_statistics_by_site():
    obs = [K.Observation("a", "flow", "cha01", "2010-01-02", 1.0), K.Observation("b", "flow", "cha01", "2010-01-03", 3.0),
           K.Observation("c", "head", "pw1", "2010-01-02", 88.0), K.Observation("d", "head", "pw1", "2010-01-03", 89.0)]
    series = {"flow:cha01": [0.5, 1.0, 3.0], "head:pw1": [87.0, 88.5, 89.5]}
    sim = K.simulated(series, ["2010-01-01", "2010-01-02", "2010-01-03"], obs)
    assert sim.tolist() == [1.0, 3.0, 88.5, 89.5]
    st = K.statistics(obs, sim)
    assert st["flow:cha01"]["nse"] == 1.0 and st["head:pw1"]["mean_error"] == pytest.approx(0.5)
    with pytest.raises(KeyError, match="no 'flow:cha09'"):
        K.simulated(series, ["2010-01-01"], [K.Observation("x", "flow", "cha09", "2010-01-01", 1.0)])


def test_forward_model(tmp_path, swatplus):
    """The forward model runs the example and writes simulated values at each observation."""
    import shutil
    d = tmp_path / "synthetic"
    shutil.copytree(ROOT / "examples" / "synthetic", d, ignore=shutil.ignore_patterns("__marimo__"))
    text = (d / "basin.toml").read_text().replace('folder = "../../runs/synthetic"', f'folder = "{RUNS / "fwd"}"')
    text = text.replace('[run]\nfolder', '[run]\ndays = 60\nfolder')
    text += '\n[outputs]\nheads = [{name = "pw1", cell = [0, 20, 28]}]\n'
    (d / "basin.toml").write_text(text)
    o = tmp_path / "obs.csv"
    o.write_text("name,kind,object,date,value,weight,source\nq1,flow,cha01,2010-02-01,0.002,1,synthetic\n"
                 "h1,head,pw1,2010-02-15,88.6,1,synthetic\n")
    out = K.run_forward_model(d / "basin.toml", None, o, tmp_path / "sim.csv")
    rows = list(csv.DictReader(open(out)))
    assert [r["name"] for r in rows] == ["q1", "h1"]
    assert 80 < float(rows[1]["simulated"]) < 100 and float(rows[0]["simulated"]) >= 0
    p = tmp_path / "pars.csv"
    p.write_text("name,model,target,value,lower,upper,transform\nk,modflow6,npf.k:all,5,0.5,50,log\n")
    out2 = K.run_forward_model(d / "basin.toml", p, o, tmp_path / "sim2.csv")
    rows2 = list(csv.DictReader(open(out2)))
    assert float(rows2[1]["simulated"]) != float(rows[1]["simulated"])  # K five times the head moves


def test_baseflow_separation():
    # a flow that only recedes is all base flow once the filter settles; a spike adds quick flow only
    q = 10 * 0.95 ** np.arange(200)
    alpha = float(np.quantile(K.recession_ratios(q), 0.9))
    assert alpha == pytest.approx(0.95)
    b = K.eckhardt(q, alpha)
    assert (b <= q + 1e-12).all() and b[-1] / q[-1] > 0.95
    spike = q.copy()
    spike[100] += 50
    bs = K.eckhardt(spike, alpha)
    assert bs[100] < spike[100] / 2 and bs.sum() < spike.sum()


# ---------- parameters ----------

def par(name, model, target, value, lower=-1e9, upper=1e9, transform="none"):
    return K.Parameter(name, model, target, value, lower, upper, transform)


def test_swatplus_targets():
    from swatplus_bmi.calibration import Change
    assert K.swatplus_change(par("cn", "swatplus", "cn2:pctchg", -10)) == Change("cn2", -10.0, "pctchg")
    assert K.swatplus_change(par("a", "swatplus", "awc:absval:layers=1-2", 0.15)) == \
        Change("awc", 0.15, "absval", layers=(1, 2))
    assert K.swatplus_change(par("c", "swatplus", "cn2:pctchg:objects=1,3:hsg=D", 5)) == \
        Change("cn2", 5.0, "pctchg", objects=(1, 3), conditions=(("hsg", "D"),))
    with pytest.raises(ValueError, match="not 'cn2'"):
        K.swatplus_change(par("x", "swatplus", "cn2", 1))
    with pytest.raises(ValueError, match="'zone' in"):
        K.swatplus_change(par("x", "swatplus", "cn2:pctchg:zone=2", 1))


def test_swatplus_parameters_go_into_calibration_cal(tmp_path, swatplus):
    import shutil
    d = tmp_path / "p"
    shutil.copytree(swatplus / "Osu_1hru", d)
    own = int((d / "calibration.cal").read_text().splitlines()[1].split()[0])  # the project's own rows come first
    K.apply_parameters("swatplus", d, [par("cn", "swatplus", "cn2:pctchg", -10)])
    cio = [line.split() for line in (d / "file.cio").read_text().splitlines()]
    assert next(f for f in cio if f[:1] == ["chg"])[2] == "calibration.cal"
    rows = (d / "calibration.cal").read_text().splitlines()
    assert int(rows[1].split()[0]) == own + 1 and rows[-1].split()[:2] == ["cn2", "pctchg"]
    with pytest.raises(ValueError, match="wouldn't make as asked"):  # swatplus-bmi's evidence: read, never applied
        K.apply_parameters("swatplus", d, [par("e", "swatplus", "evlai:absval", 3)])


def mf6_model(tmp_path, name, sfr=False):
    import datetime as dt
    import importlib.util
    spec = importlib.util.spec_from_file_location("synthetic_model", ROOT / "examples" / "synthetic" / "model.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    ws = tmp_path / name
    (m.build_small_sfr if sfr else m.build_small)(ws, 10, dt.date(2010, 1, 1))
    return ws


def load(ws):
    import flopy
    return flopy.mf6.MFSimulation.load(sim_ws=str(ws), verbosity_level=0).get_model("gwf")


def test_modflow6_multipliers_by_where(tmp_path):
    ws = mf6_model(tmp_path, "mult")
    (tmp_path / "zones.csv").write_text("layer,row,col,zone\n" + "".join(f"0,{r},0,1\n" for r in range(5)))
    before = load(ws)
    k0, sy0 = before.npf.k.array.copy(), before.sto.sy.array.copy()
    cond0 = before.riv.stress_period_data.get_data(0)["cond"].copy()
    K.apply_parameters("modflow6", ws, [par("k", "modflow6", "npf.k:all", 2.0),
                                        par("sy", "modflow6", "sto.sy:layer=0", 0.5),
                                        par("c", "modflow6", "riv.cond:zone=1", 3.0)], zones=tmp_path / "zones.csv")
    after = load(ws)
    np.testing.assert_allclose(after.npf.k.array, 2.0 * k0, rtol=1e-15)
    np.testing.assert_allclose(after.sto.sy.array, 0.5 * sy0, rtol=1e-15)
    cond = after.riv.stress_period_data.get_data(0)["cond"]
    np.testing.assert_allclose(cond, np.r_[3.0 * cond0[:5], cond0[5:]], rtol=1e-15)


def test_sfr_bed_conductivity(tmp_path):
    ws = mf6_model(tmp_path, "sfr", sfr=True)
    rhk0 = load(ws).sfr.packagedata.get_data()["rhk"].copy()
    K.apply_parameters("modflow6", ws, [par("bed", "modflow6", "sfr.rhk:all", 0.1)])
    np.testing.assert_allclose(load(ws).sfr.packagedata.get_data()["rhk"], 0.1 * rhk0, rtol=1e-15)


def test_a_multiplier_of_one_changes_nothing(tmp_path, libmf6):
    """The changed package written again, at full precision: MODFLOW 6's heads are the same bit for bit."""
    from basinledger.coupling.modflow6 import Modflow6
    heads = []
    for name, params in (("plain", []), ("one", [par("k", "modflow6", "npf.k:all", 1.0)])):
        ws = mf6_model(tmp_path, name)
        K.apply_parameters("modflow6", ws, params)
        mf = Modflow6(ws)
        try:
            for _ in range(10):
                mf.prepare_day()
                mf.set_well("wel", 0, -50.0)
                mf.solve_day()
            heads.append(mf.head())
        finally:
            mf.finalize()
    assert np.array_equal(heads[0], heads[1])


@pytest.mark.parametrize("target, refused", [
    ("npf.icelltype:all", "npf.icelltype isn't a value Basin Ledger multiplies \\(npf: k, k22, k33\\)"),
    ("npf.k:zone=1", "zone=1 needs \\[calibration\\] zones"),
    ("npf.k:layer=1", "layer 1 isn't in the grid's 1"),
    ("evt.rate:all", "has no package 'evt'"),
    ("npf.k:row=2", "where is all, layer=<n> or zone=<id>"),
])
def test_modflow6_targets_refused(tmp_path, target, refused):
    ws = mf6_model(tmp_path, "refused")
    with pytest.raises(ValueError, match=refused):
        K.apply_parameters("modflow6", ws, [par("x", "modflow6", target, 2.0)])


def test_the_examples_set_and_the_forward_models_values():
    own = [par("a", "swatplus", "cn2:pctchg", 1), par("b", "modflow6", "npf.k:all", 2)]
    assert [(p.name, p.value) for p in K.merged(own, [par("b", "modflow6", "npf.k:all", 3)])] == [("a", 1), ("b", 3)]


def test_parameters_for_a_model_not_run_are_refused():
    from conftest import synthetic_config
    from basinledger.driver import Basin
    cfg = synthetic_config()
    cfg.mf6_builder = cfg.mf6_workspace = None
    with pytest.raises(ValueError, match="parameters for modflow6, which this example doesn't run"):
        Basin(cfg, run="no-mf6", days=1, parameters=[par("k", "modflow6", "npf.k:all", 2.0)])
    with pytest.raises(NotImplementedError, match="set them in the user's config"):
        K.apply_parameters("user", RUNS, [par("c", "user", "cycles", 4.0)])


def test_a_coupled_run_with_parameters(swatplus):
    """Through the driver: a zero change to SWAT+'s cn2 and a K multiplier of one give the plain run bit for bit,
    and a real change to each moves the flow and the head."""
    from basinledger import config as C
    from basinledger.driver import Basin

    def run(name, params):
        cfg = C.load(ROOT / "examples" / "synthetic" / "two_way.toml")
        cfg.run_folder = RUNS
        b = Basin(cfg, run=name, days=120, parameters=params)
        b.run()
        out = {k: np.array(v) for k, v in b.series.items()}
        b.finalize()
        return out

    plain = run("par-plain", [])
    same = run("par-same", [par("cn", "swatplus", "cn2:pctchg", 0.0), par("k", "modflow6", "npf.k:all", 1.0)])
    moved = run("par-moved", [par("cn", "swatplus", "cn2:pctchg", -20.0), par("k", "modflow6", "npf.k:all", 0.5)])
    for k in plain:
        assert np.array_equal(plain[k], same[k]), k
    assert not np.array_equal(plain["flow:cha01"], moved["flow:cha01"])
    assert not np.array_equal(plain["head:pw1"], moved["head:pw1"])


# ---------- PEST++ ----------

def test_pest_files(tmp_path):
    p = tmp_path / "parameters.csv"
    p.write_text("name,model,target,value,lower,upper,transform\n"
                 "cn2,swatplus,cn2:pctchg,0,-25,25,none\n"
                 "riv_cond,modflow6,riv.cond:all,1,0.1,10,log\n")
    o = tmp_path / "observations.csv"
    o.write_text("name,kind,object,date,value,weight,source\n"
                 "q1,flow,cha01,2010-02-01,0.002,1,test\nq2,flow,cha01,2010-02-02,0.003,1,test\n"
                 "h1,head,pw1,2010-02-01,88.5,2,test\n")
    pst = K.write_pest_files(ROOT / "examples" / "synthetic" / "two_way.toml", p, o, tmp_path / "pest",
                             reals=20, python="python")
    d = pst.parent
    # the template, filled in as PEST++ would, is a parameter table
    tpl = (d / "parameters.csv.tpl").read_text().splitlines()
    assert tpl[0] == "ptf ~"
    filled = "\n".join(tpl[1:]).replace(f"~{'cn2':^24}~", f"{-12.5:<26}").replace(f"~{'riv_cond':^24}~", f"{3.0:<26}")
    (d / "parameters.csv").write_text(filled + "\n")
    assert [(q.name, q.value) for q in K.read_parameters(d / "parameters.csv")] == [("cn2", -12.5), ("riv_cond", 3.0)]
    # the instructions find each observation's value in the forward model's output
    ins = (d / "simulated.csv.ins").read_text().splitlines()
    assert ins == ["pif ~", "l1", "l1 ~,~ !q1!", "l1 ~,~ !q2!", "l1 ~,~ !h1!"]
    # the control file's counts and sections agree
    lines = pst.read_text().splitlines()
    assert lines[0] == "pcf" and lines[3].split() == ["2", "3", "2", "0", "2"]
    sec = {s: i for i, s in enumerate(lines) if s.startswith("*")}
    pdata = lines[sec["* parameter data"] + 1:sec["* observation groups"]]
    assert [r.split()[:6] for r in pdata] == [["cn2", "none", "relative", "0.0", "-25.0", "25.0"],
                                              ["riv_cond", "log", "factor", "1.0", "0.1", "10.0"]]
    odata = lines[sec["* observation data"] + 1:sec["* model command line"]]
    assert [r.split() for r in odata] == [["q1", "0.002", "1.0", "flow_cha01"], ["q2", "0.003", "1.0", "flow_cha01"],
                                          ["h1", "88.5", "2.0", "head_pw1"]]
    cmd = lines[sec["* model command line"] + 1]
    assert cmd.startswith("python -m basinledger forward /") and cmd.endswith("simulated.csv --run-folder runs")
    assert lines[-1] == "++ ies_num_reals(20)" and (d / "observations.csv").exists()


@pytest.mark.parametrize("rows, refused", [
    ("a b,swatplus,cn2:pctchg,0,-1,1,none\n", "PEST\\+\\+ names are 1 to 200 characters with no spaces"),
    ("K,swatplus,cn2:pctchg,0,-1,1,none\nk,modflow6,npf.k:all,1,0.5,2,none\n", "repeat"),
    ("k,modflow6,npf.k:all,1,0,2,log\n", "lower bound must be above 0"),
])
def test_pest_names_and_bounds_refused(tmp_path, rows, refused):
    p = tmp_path / "parameters.csv"
    p.write_text("name,model,target,value,lower,upper,transform\n" + rows)
    o = tmp_path / "observations.csv"
    o.write_text("name,kind,object,date,value,weight,source\nq1,flow,cha01,2010-02-01,0.002,1,test\n")
    with pytest.raises(ValueError, match=refused):
        K.write_pest_files(ROOT / "examples" / "synthetic" / "two_way.toml", p, o, tmp_path / "pest")


def test_the_forward_command(tmp_path, swatplus):
    """`basinledger forward`, as PEST++'s workers run it, with the run's files in the worker's own folder."""
    import shutil

    from basinledger.__main__ import main
    d = tmp_path / "synthetic"
    shutil.copytree(ROOT / "examples" / "synthetic", d, ignore=shutil.ignore_patterns("__marimo__", "maps"))
    cfg = d / "two_way.toml"
    cfg.write_text(cfg.read_text().replace("[run]\n", "[run]\ndays = 40\n"))
    p = tmp_path / "parameters.csv"
    p.write_text("name,model,target,value,lower,upper,transform\nk,modflow6,npf.k:all,2,0.5,5,log\n")
    o = tmp_path / "observations.csv"
    o.write_text("name,kind,object,date,value,weight,source\nh1,head,pw1,2010-02-01,88.5,1,test\n")
    assert main(["forward", str(cfg), str(p), str(o), str(tmp_path / "sim.csv"), "--run-folder",
                 str(tmp_path / "runs")]) == 0
    rows = list(csv.DictReader(open(tmp_path / "sim.csv")))
    assert rows[0]["name"] == "h1" and 80 < float(rows[0]["simulated"]) < 100
    assert (tmp_path / "runs" / "calibration").is_dir()


def test_a_restart_with_parameters_matches_the_unbroken_run(swatplus):
    """A restart writes the parameters into the new run's files again; SWAT+'s saved state already has them."""
    from basinledger import config as C
    from basinledger.driver import Basin
    params = [par("cn", "swatplus", "cn2:pctchg", -20.0), par("k", "modflow6", "npf.k:all", 0.5)]

    def cfg():
        c = C.load(ROOT / "examples" / "synthetic" / "two_way.toml")
        c.run_folder = RUNS
        return c

    a = Basin(cfg(), run="par-unbroken", days=120, parameters=params)
    a.run()
    ta = a.ledger.table()
    a.finalize()
    w = Basin(cfg(), run="par-warm", days=120, parameters=params)
    w.run(60)
    w.save(RUNS / "par-state")
    w.finalize()
    b = Basin(cfg(), run="par-restarted", days=120, restart=RUNS / "par-state", parameters=params)
    b.run()
    tb = b.ledger.table()
    b.finalize()
    for k in ta:
        assert np.array_equal(ta[k], tb[k], equal_nan=True), k


def test_kriging_weights():
    """Ordinary kriging, exponential covariance: the weights add to 1 everywhere and pick a point at the point; a
    nugget spreads them between points."""
    px, py = np.array([0.0, 100.0, 0.0]), np.array([0.0, 0.0, 100.0])
    x, y = np.array([0.0, 100.0, 50.0, 1e5]), np.array([0.0, 0.0, 50.0, 1e5])
    w = K.kriging_weights(px, py, x, y, range_=300.0)
    np.testing.assert_allclose(w.sum(axis=1), 1.0, atol=1e-12)
    np.testing.assert_allclose(w[0], [1.0, 0.0, 0.0], atol=1e-12)
    np.testing.assert_allclose(w[1], [0.0, 1.0, 0.0], atol=1e-12)
    # far from all, the points' declustered mean: K^-1 1 normalized (the corner, near both others, weighs least)
    d = np.hypot(px[:, None] - px[None, :], py[:, None] - py[None, :])
    g = np.linalg.solve(np.exp(-3.0 * d / 300.0), np.ones(3))
    np.testing.assert_allclose(w[3], g / g.sum(), atol=1e-9)
    assert w[3][0] < w[3][1] == pytest.approx(w[3][2])
    near = K.kriging_weights(px, py, np.array([10.0]), np.array([0.0]), range_=300.0)[0]
    spread = K.kriging_weights(px, py, np.array([10.0]), np.array([0.0]), range_=300.0, nugget=0.5)[0]
    assert spread[0] < near[0]  # a nugget spreads the weights off the nearest point (it still honours it at it)
    with pytest.raises(ValueError, match="range must be above 0"):
        K.kriging_weights(px, py, x, y, range_=0.0)


def test_modflow6_by_pilot_point(tmp_path):
    """Three log-transformed points on the small model's layer: K multiplied by their kriged field, in log10, each
    point's own value at its cell; a set's points on a layer, an array changed once by each set."""
    import flopy
    ws = mf6_model(tmp_path, "pp")
    grid = load(ws).modelgrid
    xc, yc = np.asarray(grid.xcellcenters), np.asarray(grid.ycellcenters)
    cells = [(0, 0), (5, 2), (9, 4)]
    with open(tmp_path / "pp.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "set", "layer", "x", "y"])
        for i, (r, c) in enumerate(cells):
            w.writerow([f"k{i}", "kpp", 0, xc[r, c], yc[r, c]])
    values = [2.0, 0.5, 4.0]
    params = [K.Parameter(f"k{i}", "modflow6", "npf.k:pp=kpp", v, 0.1, 10.0, "log") for i, v in enumerate(values)]
    k0 = load(ws).npf.k.array.copy()
    pp = {"file": tmp_path / "pp.csv", "range": 200.0, "nugget": 0.0}
    K.apply_parameters("modflow6", ws, params, pilot_points=pp)
    k = flopy.mf6.MFSimulation.load(sim_ws=str(ws), verbosity_level=0).get_model("gwf").npf.k.array
    mult = (k / k0)[0]
    for (r, c), v in zip(cells, values):
        assert mult[r, c] == pytest.approx(v, rel=1e-12)
    wts = K.kriging_weights(np.array([xc[r, c] for r, c in cells]), np.array([yc[r, c] for r, c in cells]),
                            xc.ravel(), yc.ravel(), 200.0)
    np.testing.assert_allclose(mult.ravel(), 10.0 ** (wts @ np.log10(values)), rtol=1e-12)
    with pytest.raises(ValueError, match="pp=<set> needs \\[calibration\\] pilot_points"):
        K.apply_parameters("modflow6", ws, params)
    with pytest.raises(ValueError, match="parameter kx: not a pilot point"):
        K.apply_parameters("modflow6", ws, [K.Parameter("kx", "modflow6", "npf.k:pp=kpp", 1.0, 0.1, 10.0, "log")],
                           pilot_points=pp)
    with pytest.raises(ValueError, match="all log-transformed or none"):
        K.apply_parameters("modflow6", ws, [params[0], K.Parameter("k1", "modflow6", "npf.k:pp=kpp", 1.0, 0.1,
                                                                   10.0)], pilot_points=pp)
    with pytest.raises(ValueError, match="pilot points multiply a griddata array"):
        K.apply_parameters("modflow6", ws, [K.Parameter("k0", "modflow6", "riv.cond:pp=kpp", 1.0, 0.1, 10.0, "log")],
                           pilot_points=pp)


def test_pilot_points_settings_read(tmp_path):
    from basinledger import config as C
    base = (ROOT / "examples" / "synthetic" / "two_way.toml").read_text().replace('"maps/', f'"{ROOT}/examples/synthetic/maps/')
    (tmp_path / "x.toml").write_text(base + '\n[calibration]\npilot_points = {file = "pp.csv", range = 500.0}\n')
    cfg = C.load(tmp_path / "x.toml")
    assert cfg.calibration_pilot_points == {"file": tmp_path / "pp.csv", "range": 500.0, "nugget": 0.0}
    (tmp_path / "y.toml").write_text(base + '\n[calibration]\npilot_points = "pp.csv"\n')
    with pytest.raises(ValueError, match="pilot_points is \\{file, range, nugget\\}"):
        C.load(tmp_path / "y.toml")
