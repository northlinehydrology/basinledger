"""`basinledger check`: every problem an example's start-up finds, without running a day (Basin.check).

A good example has none; one with several mistakes has each listed, the rest of its start-up going on; a run with
the same mistakes still stops at the first (tier 1).
"""

from __future__ import annotations

import shutil

import numpy as np
import pytest

from basinledger import config as C
from basinledger.__main__ import main
from basinledger.driver import Basin
from conftest import ROOT, RUNS

SYNTHETIC = ROOT / "examples" / "synthetic"


def broken(tmp_path, name, edits):
    d = tmp_path / "synthetic"
    if not d.exists():
        shutil.copytree(SYNTHETIC, d, ignore=shutil.ignore_patterns("__marimo__"))
    text = (d / name).read_text()
    for old, new in edits:
        assert old in text, old
        text = text.replace(old, new, 1)
    (d / f"broken_{name}").write_text(text)
    cfg = C.load(d / f"broken_{name}")
    cfg.run_folder = RUNS
    return cfg


def test_a_good_example_has_no_problems(swatplus):
    cfg = C.load(SYNTHETIC / "sfr.toml")
    cfg.run_folder = RUNS
    problems, notes = Basin.check(cfg)
    assert problems == []
    assert any("empty soils" in n for n in notes)  # the start-up's warnings come back as notes


def test_every_problem_is_listed(tmp_path, swatplus):
    cfg = broken(tmp_path, "two_way.toml", [
        ("cell = [0, 5, 6]", "cell = [0, 5, 16]"),                                   # the well: outside the grid
        ('heads = [{name = "pw1", cell = [0, 5, 6]}]', 'heads = [{name = "pw1", cell = [0, 56]}]'),  # DISV form
        ('map = "all:cha01"\n\n#', 'map = "all:cha09"\n\n#'),                        # no such channel
    ])
    cfg.rules = [{"name": "nobody", "kind": "curtailment", "user": "ghost", "limit": 1.0}]
    problems, _ = Basin.check(cfg)
    text = "\n".join(problems)
    assert "well pw1: cell (0, 5, 16) is outside the grid" in text
    assert "head output pw1: cell (0, 56): this model's grid is DIS" in text
    assert "rule nobody: no user 'ghost'" in text
    assert "cha09" in text
    assert len(problems) >= 4
    # a run with the same mistakes stops at the first
    with pytest.raises((ValueError, KeyError)):
        Basin(cfg, run="broken-run", days=1)


def test_an_sfr_package_with_seconds_and_a_bad_well(tmp_path, swatplus):
    cfg = broken(tmp_path, "sfr.toml", [("cell = [0, 5, 6]", "cell = [0, 50, 6]")])
    cfg.mf6_builder_args = {"unit_conversion": None}
    problems, _ = Basin.check(cfg)
    text = "\n".join(problems)
    assert "well pw1: cell (0, 50, 6) is outside the grid" in text
    assert "[exchange]: the SFR package can't be coupled" in text and "Manning's constant is 1" in text


def test_the_command(tmp_path, swatplus, capsys):
    assert main(["check", str(SYNTHETIC / "two_way.toml")]) == 0
    assert "synthetic-two-way: no problems found" in capsys.readouterr().out
    cfg = broken(tmp_path, "two_way.toml", [("cell = [0, 5, 6]", "cell = [0, 5, 16]")])
    assert main(["check", str(cfg.path)]) == 1
    out = capsys.readouterr().out
    assert "1 problem" in out and "- well pw1: cell (0, 5, 16) is outside the grid" in out


def test_the_runs_start_up_hooks_are_checked(tmp_path, swatplus):
    """A study's runs pass prepare, parameters, heads or wrap_surface; check takes the same, so it checks what the run
    runs. The project's stored calibration.cal turns flo_min back on; the study's prepare hook writes its own."""
    project = tmp_path / "Osu_1hru"
    shutil.copytree(swatplus / "Osu_1hru", project)
    cio = project / "file.cio"
    cio.write_text(cio.read_text().replace("cal_parms.cal     null", "cal_parms.cal     calibration.cal"))
    (project / "calibration.cal").write_text("calibration.cal\n 1\n NAME CHG_TYP VAL\n"
                                             "flo_min absval 5.0 0 0 0 0 0 0 0 0\n")
    cfg = C.load(SYNTHETIC / "two_way.toml")
    cfg.run_folder = RUNS
    cfg.swat_project = project
    problems, _ = Basin.check(cfg)
    assert any("flo_min" in p for p in problems)

    def prepare(model, folder):
        if model == "swatplus":
            (folder / "calibration.cal").write_text("calibration.cal\n 0\n NAME CHG_TYP VAL\n")

    assert Basin.check(cfg, prepare=prepare)[0] == []

    from basinledger.calibration import Parameter
    evt = Parameter("x", "modflow6", "evt.rate:all", 2.0, 0.0, 4.0, "none")  # the synthetic model has no EVT
    problems, _ = Basin.check(cfg, prepare=prepare, parameters=[evt])
    assert len(problems) == 1 and "has no package 'evt'" in problems[0]

    wrapped = []
    problems, _ = Basin.check(cfg, prepare=prepare, wrap_surface=lambda m: wrapped.append(m) or m)
    assert problems == [] and len(wrapped) == 1
    problems, _ = Basin.check(cfg, prepare=prepare, heads=np.zeros(3))
    assert problems == ["the start-up stopped: saved heads don't match this model's nodes"]
