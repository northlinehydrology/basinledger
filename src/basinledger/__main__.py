"""The command line.

    python -m basinledger run basin.toml [--run NAME] [--days N] [--baseline] [--out DIR]
                                         [--save-at DAY --save DIR] [--restart DIR]
    python -m basinledger status basin.toml
    python -m basinledger forward basin.toml parameters.csv observations.csv simulated.csv [--run-folder DIR]

``run`` steps the example's models through its period (or ``--days``), checking the ledger every day, and writes
ledger.csv, series.csv and summary.json to the run's outputs folder. ``--baseline`` runs it without its water users.
``--save-at`` stops after that day and saves the state to ``--save``; ``--restart`` carries on from a saved state.
``status`` prints whether the example is calibrated or a screening run. ``forward`` is the forward model PEST++
runs (calibration.py): the parameter set applied, the example run, and its value at each observation written;
``--run-folder`` puts the run's files there instead of the example's run folder.
"""

from __future__ import annotations

import argparse
import json
import sys

from . import config as C
from .driver import Basin


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="basinledger", description="Basin Ledger: coupled basin water accounting.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run an example")
    r.add_argument("config")
    r.add_argument("--run", default=None, help="run name (default: base, or baseline with --baseline)")
    r.add_argument("--days", type=int, default=None)
    r.add_argument("--baseline", action="store_true", help="without the example's water users")
    r.add_argument("--out", default=None, help="outputs folder (default: <run folder>/outputs)")
    r.add_argument("--save-at", type=int, default=None, help="stop after this many days and save the state")
    r.add_argument("--save", default=None, help="folder for the saved state")
    r.add_argument("--restart", default=None, help="carry on from a saved state")
    s = sub.add_parser("status", help="calibrated or screening run")
    s.add_argument("config")
    f = sub.add_parser("forward", help="the forward model for PEST++")
    for name in ("config", "parameters", "observations", "out"):
        f.add_argument(name)
    f.add_argument("--run-folder", default=None, help="where the run's files go (default: the example's)")
    a = ap.parse_args(argv)

    if a.cmd == "forward":
        from .calibration import run_forward_model
        run_forward_model(a.config, a.parameters, a.observations, a.out, run_folder=a.run_folder)
        return 0

    if a.cmd == "status":
        cfg = C.load(a.config)
        print(f"{cfg.name}: {cfg.calibration_status}")
        return 0

    if (a.save_at is None) != (a.save is None):
        ap.error("--save-at and --save go together")
    run = a.run or ("baseline" if a.baseline else "base")
    b = Basin.from_config(a.config, run=run, users=not a.baseline, days=a.days, restart=a.restart)
    try:
        if a.save_at is not None:
            b.run(a.save_at - b.day)
            b.save(a.save)
        else:
            b.run()
        out = b.write_outputs(a.out)
    finally:
        b.finalize()
    summary = json.loads((out / "summary.json").read_text())
    print(f"{summary['example']} ({summary['status']}): {summary['days']} days, ledger closed every day")
    print(f"outputs: {out}")
    if summary["status_note"]:
        print("note: " + summary["status_note"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
