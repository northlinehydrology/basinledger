"""MODFLOW 6 through its own BMI/XMI (libmf6), driven from Python with USGS modflowapi.

One MODFLOW 6 time step is one day. Each day the driver calls, in order:

    prepare_day()          MODFLOW 6 reads the day's stress period data (rp)
    set_recharge(...)      overrides what it read, for this day
    set_well(...)
    solve_day()            outer iterations to convergence, then the budget

Values must be set between prepare and solve. MODFLOW 6.5 and later keep a package's input in named arrays
(RCH/RECHARGE, WEL/Q) and copy them into BOUND each step, so it's the named arrays that are set here; a value written to
BOUND, as in older API examples, is overwritten before the solve.

After a solve, ``simvals(package)`` gives each entry's simulated flow (m3 d-1, positive into the aquifer), which is what
the ledger checks against what was sent.

The model must use meters and days, with one-day time steps. libmf6 holds one model per process; runs can follow one
another in a process (finalize, then a new Modflow6), but not overlap.

**Restarts.** MODFLOW 6's state, for the packages used so far, is its heads. A scenario restarts from a saved day with
a simulation covering only the days left (``remaining_simulation``) and the saved heads put straight into memory
(``set_heads``) after initialize and before the first step. MODFLOW 6 takes the start of a step's heads as its old
heads, so the restarted run continues exactly. Writing the heads to the IC file would round them. Advanced packages
with their own state (SFR stages, MAW, UZF) will need their state saved too.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
from modflowapi import ModflowApi

LIBNAME = {"darwin": "libmf6.dylib", "win32": "libmf6.dll"}.get(sys.platform, "libmf6.so")


def find_library(path: str | None = None) -> Path:
    """libmf6: the given path, $LIBMF6, or where flopy's get-modflow puts it (``get-modflow :flopy --repo modflow6``)."""
    for cand in (path, os.environ.get("LIBMF6"), Path.home() / ".local" / "share" / "flopy" / "bin" / LIBNAME):
        if cand and Path(cand).exists():
            return Path(cand)
    raise FileNotFoundError(f"{LIBNAME} not found: set LIBMF6, or run `uv run get-modflow :flopy --repo modflow6 "
                            f"--subset mf6,{LIBNAME}`")


class Modflow6:
    def __init__(self, workspace: str | Path, model: str = "gwf", library: str | None = None):
        self.ws = Path(workspace).resolve()
        self.model = model.upper()
        self.api = ModflowApi(str(find_library(library)), working_directory=str(self.ws))
        self.api.initialize()
        self.day = 0
        if int(self.get("ITMUNI", "", "TDIS")[0]) != 4:
            raise ValueError("MODFLOW 6 model must use days (TDIS TIME_UNITS days)")
        if int(self.get("LENUNI", self.model, "DIS")[0]) != 2:
            raise ValueError("MODFLOW 6 model must use meters (DIS LENGTH_UNITS meters)")
        self.cell_area = self.get("AREA", self.model, "DIS").astype(float)
        self.nodes = len(self.cell_area)
        shape = self.get("MSHAPE", self.model, "DIS")
        self.shape = tuple(int(s) for s in shape)  # (nlay, nrow, ncol) for DIS
        reduced = self._try("NODEREDUCED", self.model, "DIS")
        # NODEREDUCED is only allocated (to the user node count) when IDOMAIN removes cells
        nuser = int(np.prod(self.shape))
        self._reduced = reduced.astype(int) if reduced is not None and len(reduced) == nuser else None
        self._rows: dict[str, np.ndarray] = {}

    # ---------- memory access ----------
    def address(self, var: str, model: str, package: str) -> str:
        return self.api.get_var_address(var, model, package) if model else self.api.get_var_address(var, package)

    def get(self, var: str, model: str, package: str) -> np.ndarray:
        return self.api.get_value(self.address(var, model, package)).copy()

    def _try(self, var, model, package):
        try:
            return self.get(var, model, package)
        except Exception:
            return None

    def set(self, var: str, package: str, values: np.ndarray):
        self.api.set_value(self.address(var, self.model, package.upper()), values)

    # ---------- cells ----------
    def node(self, layer: int, row: int, col: int) -> int:
        """0-based (layer, row, col) of a DIS grid -> 0-based reduced node."""
        nlay, nrow, ncol = self.shape
        user = layer * nrow * ncol + row * ncol + col
        if self._reduced is None:
            return user
        n = int(self._reduced[user])
        if n <= 0:
            raise ValueError(f"cell {(layer, row, col)} is inactive")
        return n - 1

    def package_nodes(self, package: str) -> np.ndarray:
        """0-based reduced node of each entry in a boundary package, in its order."""
        return self.get("NODELIST", self.model, package.upper()).astype(int)[: self.nbound(package)] - 1

    def nbound(self, package: str) -> int:
        return int(self.get("NBOUND", self.model, package.upper())[0])

    def entry(self, package: str, node: int) -> int:
        idx = np.flatnonzero(self.package_nodes(package) == node)
        if len(idx) != 1:
            raise ValueError(f"{package}: {len(idx)} entries at node {node}")
        return int(idx[0])

    # ---------- the day ----------
    def prepare_day(self):
        # MODFLOW 6 takes its step length from TDIS (known only once the step is prepared), so it is checked after
        self.api.prepare_time_step(1.0)
        dt = self.api.get_time_step()
        if abs(dt - 1.0) > 1e-12:
            raise ValueError(f"MODFLOW 6 time step is {dt} days; Basin Ledger steps one day at a time")

    def set_recharge(self, package: str, rate_per_node: np.ndarray):
        """Recharge as a depth rate per node (m d-1); entries of the package take their node's rate."""
        nodes = self.package_nodes(package)
        self.set("RECHARGE", package, rate_per_node[nodes].astype(float))

    def set_well(self, package: str, entry: int, q: float):
        """A well's rate, m3 d-1 (negative pumps out)."""
        qq = self.get("Q", self.model, package.upper())
        qq[entry] = q
        self.set("Q", package, qq)

    def set_stage(self, package: str, stage: np.ndarray):
        """River stages, m, one per entry of a river package."""
        self.set("STAGE", package, np.asarray(stage, dtype=float))

    def solve_day(self, callback=None, max_rounds: int = 30):
        """Solve the day prepared. ``callback()``, if given, runs after each outer iteration and returns whether it
        changed a boundary; the day is done once the solve has converged and the callback changes nothing (or it
        has changed something ``max_rounds`` times after convergence). How a coupling sets a boundary from the heads
        within the day: moving the boundary with the iterations, rather than after each convergence, takes about
        half as many iterations.

        A simulation with transport models (``coupling/transport.py``) has a solution for each, solved after the
        flow's: ``do_time_step`` solves them all, and with a callback they are solved here once the flow is done."""
        if callback is None:
            self.api.do_time_step()
        else:
            self.api.prepare_solve(1)
            iterations, rounds = 0, 0
            while True:
                converged = self.api.solve(1)
                iterations += 1
                changed = callback()
                if converged and not changed:
                    break
                if converged:
                    rounds += 1
                    if rounds >= max_rounds:
                        break
                if iterations >= 2000:
                    raise RuntimeError(f"MODFLOW 6 didn't converge on day {self.day}")
            self.api.finalize_solve(1)
            for sol in range(2, self.api.get_subcomponent_count() + 1):
                self._solve(sol)
        self.api.finalize_time_step()
        self.day += 1

    def _solve(self, sol: int, max_iterations: int = 2000):
        """One solution's outer iterations to convergence (as modflowapi's ``do_time_step`` does for each)."""
        self.api.prepare_solve(sol)
        for _ in range(max_iterations):
            if self.api.solve(sol):
                break
        else:
            raise RuntimeError(f"MODFLOW 6 solution {sol} didn't converge on day {self.day}")
        self.api.finalize_solve(sol)

    def simvals(self, package: str) -> np.ndarray:
        """Each entry's simulated flow for the last day, m3 d-1, positive into the aquifer."""
        return self.get("SIMVALS", self.model, package.upper())[: self.nbound(package)].astype(float)

    def set_heads(self, heads: np.ndarray):
        """Heads, m, per reduced node: call after initialize and before the first prepare_day."""
        addr = self.api.get_var_address("X", self.model)
        if self.api.get_value(addr).shape != np.shape(heads):
            raise ValueError("saved heads don't match this model's nodes")
        self.api.set_value(addr, np.asarray(heads, dtype=np.float64))

    def head(self) -> np.ndarray:
        """Heads, m, per reduced node."""
        return self.api.get_value(self.api.get_var_address("X", self.model)).copy()

    @property
    def end(self) -> float:
        return self.api.get_end_time()

    @property
    def time(self) -> float:
        return self.api.get_current_time()

    def finalize(self):
        self.api.finalize()


def remaining_simulation(src: Path, dst: Path, day: int):
    """A copy of a MODFLOW 6 simulation covering only the days from ``day`` on, for a restart. Supports one stress
    period of daily steps (what the driver's builders write); models with several stress periods need their period
    data re-keyed, which isn't built yet."""
    import shutil

    import flopy

    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    sim = flopy.mf6.MFSimulation.load(sim_ws=str(dst), verbosity_level=0)
    pd = sim.tdis.perioddata.get_data()
    if len(pd) != 1:
        raise NotImplementedError("restarts of MODFLOW 6 models with several stress periods aren't built yet")
    perlen, nstp, tsmult = float(pd[0][0]), int(pd[0][1]), float(pd[0][2])
    if abs(perlen - nstp) > 1e-9 or tsmult != 1.0 or not 0 < day < nstp:
        raise ValueError(f"restart at day {day}: needs one period of daily steps longer than that")
    sim.tdis.perioddata = [(perlen - day, nstp - day, 1.0)]
    sim.write_simulation(silent=True)
    return dst
