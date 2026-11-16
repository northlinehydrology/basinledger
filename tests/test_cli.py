"""The command line: run, save, restart, status."""

from __future__ import annotations

import csv
import json
import shutil

from basinledger.__main__ import main
from conftest import ROOT, RUNS


def copy_example(tmp_path):
    d = tmp_path / "synthetic"
    shutil.copytree(ROOT / "examples" / "synthetic", d, ignore=shutil.ignore_patterns("__marimo__", "*.pyc"))
    text = (d / "basin.toml").read_text().replace('folder = "../../runs/synthetic"', f'folder = "{(RUNS / "cli").as_posix()}"')
    (d / "basin.toml").write_text(text)
    return d / "basin.toml"


def test_run_writes_outputs(tmp_path, swatplus, capsys):
    cfg = copy_example(tmp_path)
    assert main(["run", str(cfg), "--days", "30", "--out", str(tmp_path / "out")]) == 0
    s = json.loads((tmp_path / "out" / "summary.json").read_text())
    assert s["status"] == "screening" and s["days"] == 30 and "screening run" in s["status_note"]
    rows = list(csv.DictReader(open(tmp_path / "out" / "ledger.csv")))
    assert len(rows) == 30 and rows[0]["date"] == "2010-01-01" and "recharge:sent" in rows[0]
    assert "screening" in capsys.readouterr().out


def test_save_and_restart(tmp_path, swatplus):
    cfg = copy_example(tmp_path)
    assert main(["run", str(cfg), "--days", "40", "--save-at", "20", "--save", str(tmp_path / "state"),
                 "--out", str(tmp_path / "first")]) == 0
    assert main(["run", str(cfg), "--days", "40", "--restart", str(tmp_path / "state"), "--run", "second",
                 "--out", str(tmp_path / "second")]) == 0
    assert main(["run", str(cfg), "--days", "40", "--run", "whole", "--out", str(tmp_path / "whole")]) == 0
    a = list(csv.DictReader(open(tmp_path / "second" / "ledger.csv")))
    b = list(csv.DictReader(open(tmp_path / "whole" / "ledger.csv")))
    assert a == b


def test_status(capsys):
    assert main(["status", str(ROOT / "examples" / "synthetic" / "basin.toml")]) == 0
    assert capsys.readouterr().out.strip() == "synthetic: screening"
