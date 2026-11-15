"""Pleasant Lake's lake, as the lake link reads it (coupling/lakes.py): 1,382 vertical connections in the LGR inset, no
table, no outlets. On a copy in daily steps, the steady first day, then the first transient day: the change in the
lake's volume at its stage, as the link computes it from the connections, is MODFLOW 6's own lake storage (QSTO), and
the stage at a volume inverts it. Skipped when the model isn't fetched (examples/central-sands/fetch.py).
"""

from __future__ import annotations

import importlib.util
import shutil

import pytest
from conftest import ROOT, RUNS

from basinledger.coupling import lakes as LK
from basinledger.coupling.modflow6 import Modflow6
from basinledger.coupling.timesteps import daily_steps

spec = importlib.util.spec_from_file_location("central_sands_fetch", ROOT / "examples" / "central-sands" / "fetch.py")
F = importlib.util.module_from_spec(spec)
spec.loader.exec_module(F)
WS = F.workspace()
pytestmark = pytest.mark.skipif(not (WS / "mfsim.nam").exists(),
                                reason="Pleasant Lake not fetched: uv run python examples/central-sands/fetch.py")


def test_pleasant_lakes_volume_is_modflow_6s(libmf6):
    ws = RUNS / "pleasant-lake"
    if ws.exists():
        shutil.rmtree(ws)
    shutil.copytree(WS, ws)
    daily_steps(ws)
    mf = Modflow6(ws, "plsnt_lgr_inset")
    try:
        assert LK.problems(mf, [LK.LakeConfig("lake1", "lak_0")], ["lake1"], None) == []
        mf.prepare_day()
        mf.solve_day()  # the steady first day
        mf.prepare_day()
        c = LK.curve(mf, "lak_0", 1)
        assert len(c.sarea) == 1382
        mf.solve_day()
        old = float(mf.get("XOLDPAK", mf.model, "LAK_0")[0])
        new = float(mf.get("XNEWPAK", mf.model, "LAK_0")[0])
        qsto = float(mf.get_model_units("QSTO", mf.model, "LAK_0")[0])
        assert old != new
        assert c.volume(new) - c.volume(old) == pytest.approx(-qsto, rel=1e-8)
        assert c.stage(c.volume(new)) == pytest.approx(new, abs=1e-8)
        assert 1e6 < c.volume(new) < 1e8  # Pleasant Lake: about 0.5 km2, several meters deep
    finally:
        mf.finalize()
