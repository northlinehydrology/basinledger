"""`basinledger check`: every problem an example's start-up finds, without running a day (Basin.check).

A good example has none; one with several mistakes has each listed, the rest of its start-up going on; a run with
the same mistakes still stops at the first (tier 1).
"""

from __future__ import annotations

import shutil

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
