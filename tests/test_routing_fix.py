"""swatplus-bmi's routing fix passed through (``[swatplus] routing_fix``; docs/design.md, The routing fix passed
through to swatplus-bmi):
- the setting, read and checked; unset by default;
- swatplus-bmi's start file: the project and each option given, in order; with none, no file, so a run without it
  starts as before; the water-table link's file as it was;
- on SWAT+ (swatplus-bmi 0.18.0), the synthetic example with Muskingum routing (rte_cha 1): the key reaches the
  start file and SWAT+ routes differently; unset, no start file;
- refused on a swatplus-bmi before 0.18.0;
- a restart: swatplus-bmi keeps the setting a state was saved with, so a restart whose config says otherwise is
  refused, and one that agrees matches the unbroken run to the bit."""

from __future__ import annotations

import json
import warnings

import numpy as np
import pytest
from conftest import RUNS, SYNTHETIC, synthetic_config

from basinledger import config as C
from basinledger.coupling import swatplus as S
from basinledger.driver import Basin


def _load(tmp_path, line: str):
    text = SYNTHETIC.read_text().replace('[swatplus]\n', f'[swatplus]\n{line}\n', 1)
    (tmp_path / "basin.toml").write_text(text)
    return C.load(tmp_path / "basin.toml")


# ---------- the setting ----------

def test_unset_by_default():
    assert C.load(SYNTHETIC).routing_fix is None


@pytest.mark.parametrize("fix", ["muskingum_daily", "none"])
def test_the_setting_is_read(tmp_path, fix):
    assert _load(tmp_path, f'routing_fix = "{fix}"').routing_fix == fix


@pytest.mark.parametrize("value", ['"all"', '"muskingum"', "true"])
def test_another_value_is_refused(tmp_path, value):
    with pytest.raises(ValueError, match=r"\[swatplus\] routing_fix is one of"):
        _load(tmp_path, f"routing_fix = {value}")


# ---------- swatplus-bmi's start file ----------

def test_the_start_file_holds_the_project_and_each_option_in_order(tmp_path):
    p = (tmp_path / "proj").resolve()
    f = S.write_start_file(tmp_path / "a.toml", p, held_lateral="none", routing_fix="muskingum_daily")
    assert f.read_text() == f'project = "{p.as_posix()}"\nheld_lateral = "none"\nrouting_fix = "muskingum_daily"\n'
    f = S.write_start_file(tmp_path / "b.toml", p, routing_fix="muskingum_daily")
    assert f.read_text() == f'project = "{p.as_posix()}"\nrouting_fix = "muskingum_daily"\n'


def test_the_water_table_links_start_file_is_as_it_was(tmp_path):
    """Without routing_fix the water-table link's file is byte for byte what 0.20.0's branch wrote before."""
    p = (tmp_path / "proj").resolve()
    f = S.write_start_file(tmp_path / "c.toml", p, held_lateral="swat")
    assert f.read_text() == f'project = "{p.as_posix()}"\nheld_lateral = "swat"\n'


# ---------- on SWAT+: the synthetic example with Muskingum routing ----------

DAYS = 365


@pytest.fixture(scope="module")
def fix(swatplus):
    if not S.has_routing_fix():
        pytest.skip("the routing fix needs swatplus-bmi 0.18.0")
    return swatplus


def muskingum(model, folder):
    """Osu_1hru's channel routed by Muskingum (rte_cha 1 in codes.bsn; the example has 0, variable storage)."""
    if model != "swatplus":
        return
    p = folder / "codes.bsn"
    lines = p.read_text().splitlines()
    head, vals = lines[1].split(), lines[2].split()
    vals[head.index("rte_cha")] = "1"
    lines[2] = "  ".join(vals)
    p.write_text("\n".join(lines) + "\n")


def _basin(name, routing_fix=None, days=DAYS, restart=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return Basin(synthetic_config(routing_fix=routing_fix), run=name, days=days, prepare=muskingum,
                     restart=restart)


def _outflow(b, days):
    """cha01's outflow, m3 s-1, day by day, from SWAT+."""
    v = np.zeros(b.swat.get_grid_size(b.swat.get_var_grid(S.QOUT)))
    out = []
    for _ in range(days):
        b.step()
        b.swat.get_value(S.QOUT, v)
        out.append(v[0])
    return np.array(out)


def _run(name, routing_fix=None):
    b = _basin(name, routing_fix)
    try:
        q = _outflow(b, DAYS)
        b.write_outputs()
        return q, b.folder
    finally:
        b.finalize()


@pytest.mark.slow
def test_the_fix_reaches_the_start_file_and_swatplus_and_unset_there_is_none(fix):
    stock, off = _run("rf-unset")
    fixed, on = _run("rf-fixed", "muskingum_daily")
    assert not (off / "swatplus-bmi.toml").exists()          # unset: SWAT+ started from the project, as before
    assert "routing_fix" not in json.loads((off / "outputs" / "summary.json").read_text())
    swat = (on / "swatplus").resolve().as_posix()
    assert (on / "swatplus-bmi.toml").read_text() == f'project = "{swat}"\nrouting_fix = "muskingum_daily"\n'
    assert json.loads((on / "outputs" / "summary.json").read_text())["routing_fix"] == "muskingum_daily"
    assert np.abs(fixed - stock).max() > 1e-3                 # SWAT+ routed it: Muskingum keeps its substeps
    t = np.arange(DAYS)
    assert (t * fixed).sum() / fixed.sum() != (t * stock).sum() / stock.sum()


def test_a_swatplus_bmi_before_0_18_is_refused(fix, monkeypatch):
    monkeypatch.setattr(S, "_swatplus_bmi_version", lambda: (0, 17))
    with pytest.raises(ValueError, match="routing_fix needs swatplus-bmi 0.18.0"):
        _basin("rf-old", "muskingum_daily", days=10)


@pytest.mark.slow
def test_a_restart_keeps_the_saved_setting_and_refuses_another(fix):
    n, k = 120, 50
    a = _basin("rf-unbroken", "muskingum_daily", days=n)
    try:
        unbroken = _outflow(a, n)
    finally:
        a.finalize()
    for fix_, state in (("muskingum_daily", "rf-state-fixed"), (None, "rf-state-stock")):
        w = _basin(f"{state}-warm", fix_, days=n)
        try:
            w.run(k)
            w.save(RUNS / state)
        finally:
            w.finalize()
    # swatplus-bmi's state carries the setting and a load takes it, so a config that says otherwise is refused
    with pytest.raises(ValueError, match="routing_fix is 'none' but the saved state"):
        _basin("rf-restart-unset", None, days=n, restart=RUNS / "rf-state-fixed")
    with pytest.raises(ValueError, match="routing_fix is 'muskingum_daily' but the saved state"):
        _basin("rf-restart-fixed-on-stock", "muskingum_daily", days=n, restart=RUNS / "rf-state-stock")
    b = _basin("rf-restarted", "muskingum_daily", days=n, restart=RUNS / "rf-state-fixed")
    try:
        assert b.day == k
        restarted = _outflow(b, n - k)
    finally:
        b.finalize()
    assert restarted.tobytes() == unbroken[k:].tobytes()

