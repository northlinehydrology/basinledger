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

**Units.** The model stays in its own units; everything this class gives or takes is in the ledger's: meters, cubic
meters and days (kilograms for a transport model's mass). The model's units are DIS's LENGTH_UNITS (feet, meters or
centimeters) and TDIS's TIME_UNITS (seconds, minutes, hours or days). Where either is undefined, as in a model
converted by mf5to6, the study states it (``[modflow6] length_units``, ``time_units``); a stated unit that differs
from the model's own is refused, and so are years (a day isn't a whole step of them). ``get`` and ``set`` convert
each variable by its dimension (``DIMENSION``): a head or stage by meters per unit, a conductance by its square over
days per unit, a flow by its cube over days per unit, a concentration by its inverse cube. In a model not in meters
and days, a real-valued variable whose dimension isn't listed is refused rather than passed through unconverted.
Integer variables (counts, node numbers, status) are never converted. In meters and days nothing is converted.

**Cells.** A structured grid (DIS) names a cell by ``(layer, row, col)``, a vertex grid (DISV, such as a quadtree
refined around a site) by ``(layer, icell2d)``, all 0-based as flopy numbers them; ``node`` turns either into the
reduced node, and a cell in the other grid's form is refused. An unstructured grid (DISU) has no layers in its cell
numbers, so maps on it give reduced nodes (the ``node`` column) and its wells and heads aren't supported yet.

Time steps must be one day long: the driver rewrites a model's time discretization to one-day steps in the run's
copy first (``timesteps.daily_steps``), each stress period kept. libmf6 holds one model per process; runs can follow one
another in a process (finalize, then a new Modflow6), but not overlap.

**Restarts.** MODFLOW 6's state, for the packages used so far, is its heads. A scenario restarts from a saved day with
a simulation covering only the days left (``remaining_simulation``: the periods over left out, the rest renumbered)
and the saved heads put straight into memory (``set_heads``) after initialize and before the first step. MODFLOW 6
takes the start of a step's heads as its old heads, so the restarted run continues exactly. Writing the heads to the
IC file would round them. SFR without its storage option keeps no water from day to day; advanced packages with their
own state (MAW, UZF, LAK, SFR with storage) will need their state saved too.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
from modflowapi import ModflowApi

LIBNAME = {"darwin": "libmf6.dylib", "win32": "libmf6.dll"}.get(sys.platform, "libmf6.so")
GRIDS = {3: "dis", 2: "disv", 1: "disu"}  # by the length of the discretization's shape (MSHAPE)
CELL_FORM = {"dis": "layer, row, col", "disv": "layer, icell2d", "disu": "reduced nodes (give node)"}

LENUNI = {1: "feet", 2: "meters", 3: "centimeters"}                 # DIS LENGTH_UNITS, 0 undefined
ITMUNI = {1: "seconds", 2: "minutes", 3: "hours", 4: "days", 5: "years"}  # TDIS TIME_UNITS, 0 undefined
METERS = {"feet": 0.3048, "meters": 1.0, "centimeters": 0.01}       # per unit (the international foot)
DAYS = {"seconds": 1.0 / 86400.0, "minutes": 1.0 / 1440.0, "hours": 1.0 / 24.0, "days": 1.0}  # per unit

L, A, V, Q, LT, C, MT, ONE = (1, 0), (2, 0), (3, 0), (3, -1), (1, -1), (-3, 0), (0, -1), (0, 0)
# (length, time) exponents of each variable Basin Ledger reads or sets, by name; (package type, name) where a name
# means different things in different packages, and ("gwt", name) in a transport model (mass taken as kg)
DIMENSION = {
    "X": L, "TOP": L, "BOT": L, "AREA": A, "K": LT, "K22": LT, "K33": LT, "SY": ONE, "SS": (-1, 0),
    "STRGSS": Q, "STRGSY": Q, "SIMVALS": Q, "Q": Q, "COND": (2, -1), "RECHARGE": LT,
    "STAGE": L, "RBOT": L, "ELEV": L, "BHEAD": L, "HEAD": L, "SURFACE": L, "DEPTH": L,
    "INFLOW": Q, "RUNOFF": Q, "QOUTFLOW": Q, "QEXTOUTFLOW": Q, "QEXTINFLOW": Q, "QFROMMVR": Q, "QTOMVR": Q,
    "RAIN": LT, "EVAP": LT, "STRTOP": L, "SSTAGE": L, "LENGTH": L, "WIDTH": L, "USTRF": ONE, "SLOPE": ONE,
    "RAINFALL": LT, "EVAPORATION": LT, "WITHDRAWAL": Q, "XNEWPAK": L, "XOLDPAK": L,
    "RATESIM": Q, ("maw", "RATE"): Q, ("evt", "RATE"): LT,
    ("lak", "TELEV"): L, ("lak", "BELEV"): L, ("lak", "SAREA"): A, ("lak", "LAKEBOT"): L, ("lak", "LAKETOP"): L,
    ("lak", "STRT"): L,
    ("gwt", "X"): C, ("gwt", "SIMVALS"): MT, ("gwt", "FLOW"): Q, ("gwt", "AUXVAR"): C,
    ("gwt", "RATESTO"): MT, ("gwt", "RATESRB"): MT,
}


def units_factor(dimension: tuple[int, int], meters: float, days: float) -> float:
    """What a value of this (length, time) dimension in the model's units is multiplied by to be in meters and days."""
    return meters ** dimension[0] * days ** dimension[1]


def find_library(path: str | None = None) -> Path:
    """libmf6: the given path, $LIBMF6, or where flopy's get-modflow puts it (``get-modflow :flopy --repo modflow6``)."""
    for cand in (path, os.environ.get("LIBMF6"), Path.home() / ".local" / "share" / "flopy" / "bin" / LIBNAME):
        if cand and Path(cand).exists():
            return Path(cand)
    raise FileNotFoundError(f"{LIBNAME} not found: set LIBMF6, or run `uv run get-modflow :flopy --repo modflow6 "
                            f"--subset mf6,{LIBNAME}`")


class Modflow6:
    def __init__(self, workspace: str | Path, model: str = "gwf", library: str | None = None,
                 length_units: str | None = None, time_units: str | None = None):
        """``length_units``, ``time_units``: the model's, where its DIS or TDIS leaves them undefined (``[modflow6]``
        in the config); a unit stated that differs from the model's own is refused."""
        self.ws = Path(workspace).resolve()
        self.model = model.upper()
        self.api = ModflowApi(str(find_library(library)), working_directory=str(self.ws))
        self.api.initialize()
        self.day = 0
        self.meters = self.days = 1.0  # nothing converted while the units are read
        self._types: dict[tuple[str, str], str] = {}
        try:
            self.length_units, self.length_stated = self._units(
                "length", LENUNI, int(self.get("LENUNI", self.model, "DIS")[0]), length_units, METERS,
                "DIS LENGTH_UNITS", "length_units")
            self.time_units, self.time_stated = self._units(
                "time", ITMUNI, int(self.get("ITMUNI", "", "TDIS")[0]), time_units, DAYS,
                "TDIS TIME_UNITS", "time_units")
        except ValueError:
            self.api.finalize()
            raise
        self.meters, self.days = METERS[self.length_units], DAYS[self.time_units]
        self._read_grid()

    def _read_grid(self):
        """The model's grid: cell areas, shape, and the user-to-reduced node numbering."""
        self.cell_area = self.get("AREA", self.model, "DIS").astype(float)
        self.nodes = len(self.cell_area)
        shape = self.get("MSHAPE", self.model, "DIS")
        self.shape = tuple(int(s) for s in shape)  # (nlay, nrow, ncol) for DIS, (nlay, ncpl) for DISV, (nodes,) DISU
        self.grid = GRIDS[len(self.shape)]
        reduced = self._try("NODEREDUCED", self.model, "DIS")
        # NODEREDUCED is only allocated (to the user node count) when IDOMAIN removes cells
        nuser = int(np.prod(self.shape))
        self._reduced = reduced.astype(int) if reduced is not None and len(reduced) == nuser else None
        self._rows: dict[str, np.ndarray] = {}

    def view(self, model: str) -> Modflow6:
        """Another groundwater model of the same simulation, through the same libmf6 (one simulation per process):
        its own grid and packages, the simulation's time and units. A view doesn't step the simulation. Models in
        other length units than this one are refused (the conversion is one per run)."""
        v = object.__new__(Modflow6)
        v.__dict__.update({k: x for k, x in self.__dict__.items() if k not in ("_types",)})
        v.model = model.upper()
        v._types = {}
        try:
            code = int(v.get_model_units("LENUNI", v.model, "DIS")[0])
        except Exception as e:  # noqa: BLE001  (no such model)
            raise ValueError(f"[modflow6] models: {model.lower()!r} isn't a groundwater model of the simulation") from e
        own = LENUNI.get(code)
        if own is not None and own != self.length_units:
            raise ValueError(f"[modflow6] models: {model} is in {own} and {self.model.lower()} in "
                             f"{self.length_units}; coupled models must share their length units")
        v._read_grid()
        return v

    # ---------- units ----------
    @staticmethod
    def _units(what, codes, code, stated, known, where, key):
        own = codes.get(code)
        if stated is not None and stated not in known:
            raise ValueError(f"[modflow6] {key} = {stated!r}: Basin Ledger converts {', '.join(known)}")
        if own is None:
            if stated is None:
                raise ValueError(f"MODFLOW 6 model's {what} units are undefined ({where}): state them in the config, "
                                 f"[modflow6] {key} = one of {', '.join(repr(k) for k in known)}")
            return stated, True
        if own not in known:
            raise ValueError(f"MODFLOW 6 model in {own} ({where}): a one-day step isn't a whole number of them; "
                             f"Basin Ledger converts {', '.join(known)}")
        if stated is not None and stated != own:
            raise ValueError(f"[modflow6] {key} = {stated!r}, but the model's {where} is {own}")
        return own, False

    @property
    def si(self) -> bool:
        """Whether the model is in meters and days, the ledger's units (nothing converted)."""
        return self.meters == 1.0 and self.days == 1.0

    def factor(self, dimension: tuple[int, int]) -> float:
        """Model units to meters and days, for a (length, time) dimension."""
        return units_factor(dimension, self.meters, self.days)

    def units_summary(self) -> str:
        """The model's units and the factors applied at every link, for the start-up check."""
        stated = [w for w, s in (("length", self.length_stated), ("time", self.time_stated)) if s]
        where = f" ({' and '.join(stated)} stated in [modflow6], undefined in the model)" if stated else ""
        head = f"{self.length_units} and {self.time_units}{where}"
        if self.si:
            return f"{head}, the ledger's units: nothing converted"
        f = lambda d: f"{self.factor(d):.6g}"  # noqa: E731
        return (f"{head}; converted at every link to meters and days: lengths and heads x {f(L)}, areas x {f(A)}, "
                f"flows x {f(Q)} (to m3 d-1), depth rates x {f(LT)} (to m d-1), conductances x {f((2, -1))} "
                f"(to m2 d-1), concentrations x {f(C)} (to kg m-3)")

    def _ptype(self, model: str, package: str) -> str:
        key = (model, package)
        if key not in self._types:
            if model != self.model and self._has(model, "MST"):
                self._types[key] = "gwt"
            elif package:
                try:
                    self._types[key] = str(self.api.get_value(self.api.get_var_address(
                        "PACKAGE_TYPE", model, package))[0]).strip().lower()
                except Exception:
                    self._types[key] = ""
            else:
                self._types[key] = ""
        return self._types[key]

    def _has(self, model: str, package: str) -> bool:
        prefix = f"{model}/{package}/"
        return any(n.startswith(prefix) for n in self.api.get_input_var_names()) or \
            any(n.startswith(prefix) for n in self.api.get_output_var_names())

    def _dimension(self, var: str, model: str, package: str) -> tuple[int, int]:
        t = self._ptype(model, package) if model else ""
        d = DIMENSION.get((t, var)) or (DIMENSION.get(var) if t != "gwt" else None)
        if d is None:
            raise ValueError(f"{model}/{package}/{var}: Basin Ledger doesn't know its dimension, so it can't convert "
                             f"it from {self.length_units} and {self.time_units} (coupling/modflow6.py, DIMENSION)")
        return d

    def _scale(self, values: np.ndarray, var: str, model: str, package: str, inverse: bool = False):
        if self.si or not np.issubdtype(np.asarray(values).dtype, np.floating):
            return values
        f = self.factor(self._dimension(var, model, package))
        return values / f if inverse else values * f

    # ---------- memory access ----------
    def address(self, var: str, model: str, package: str) -> str:
        if not model:
            return self.api.get_var_address(var, package)
        return self.api.get_var_address(var, model, package) if package else self.api.get_var_address(var, model)

    def get(self, var: str, model: str, package: str) -> np.ndarray:
        """A variable from MODFLOW 6's memory, in meters and days (``package`` "" for one of the model's own)."""
        v = self.api.get_value(self.address(var, model, package)).copy()
        return self._scale(v, var, model.upper(), package.upper())

    def get_model_units(self, var: str, model: str, package: str) -> np.ndarray:
        """A variable as MODFLOW 6 holds it, in the model's units."""
        return self.api.get_value(self.address(var, model, package)).copy()

    def _try(self, var, model, package):
        try:
            return self.get(var, model, package)
        except Exception:
            return None

    def set(self, var: str, package: str, values: np.ndarray, model: str | None = None):
        """A variable of the coupled model's package (or ``model``'s), given in meters and days."""
        m = (model or self.model).upper()
        values = np.asarray(values)
        self.api.set_value(self.address(var, m, package.upper()), self._scale(values, var, m, package.upper(), True))

    # ---------- cells ----------
    def node(self, *cell: int, model: str | None = None, node: int | None = None) -> int:
        """0-based ``(layer, row, col)`` of a DIS grid or ``(layer, icell2d)`` of a DISV grid -> 0-based reduced
        node. A cell in the other grid's form, outside the grid or inactive is refused. ``model``, if given, must be
        this one (``Models`` takes the others); ``node`` instead of a cell is a reduced node, checked."""
        if model is not None and model.upper() != self.model:
            raise ValueError(f"a cell of model {model}, but the model coupled is {self.model.lower()} (several models "
                             "are coupled with [modflow6] models)")
        if node is not None:
            if not 0 <= int(node) < self.nodes:
                raise ValueError(f"node {node} is outside the model's {self.nodes} nodes")
            return int(node)
        cell = tuple(int(c) for c in cell)
        if self.grid == "disu" or len(cell) != len(self.shape):
            raise ValueError(f"cell {cell}: this model's grid is {self.grid.upper()}, whose cells are "
                             f"{CELL_FORM[self.grid]}")
        if any(c < 0 or c >= n for c, n in zip(cell, self.shape)):
            raise ValueError(f"cell {cell} is outside the grid ({self.grid.upper()} {self.shape})")
        user = int(np.ravel_multi_index(cell, self.shape))
        if self._reduced is None:
            return user
        n = int(self._reduced[user])
        if n <= 0:
            raise ValueError(f"cell {cell} is inactive")
        return n - 1

    def package_nodes(self, package: str) -> np.ndarray:
        """0-based reduced node of each entry in a boundary package, in its order."""
        return self.get("NODELIST", self.model, package.upper()).astype(int)[: self.nbound(package)] - 1

    def nbound(self, package: str) -> int:
        return int(self.get("NBOUND", self.model, package.upper())[0])

    def has_package(self, package: str) -> bool:
        return self._has(self.model, package.upper())

    def aux_multiplier(self, package: str) -> bool:
        """Whether a boundary package multiplies its values by an auxiliary variable."""
        return int(self.get("IAUXMULTCOL", self.model, package.upper())[0]) > 0

    def entry(self, package: str, node: int) -> int:
        idx = np.flatnonzero(self.package_nodes(package) == node)
        if len(idx) != 1:
            raise ValueError(f"{package}: {len(idx)} entries at node {node}")
        return int(idx[0])

    # ---------- the day ----------
    def prepare_day(self):
        # MODFLOW 6 takes its step length from TDIS (known only once the step is prepared), so it is checked after
        self.api.prepare_time_step(1.0)
        dt = self.api.get_time_step() * self.days
        if abs(dt - 1.0) > 1e-12:
            raise ValueError(f"MODFLOW 6 time step is {dt} days; Basin Ledger steps one day at a time (a run "
                             f"rewrites the model's time discretization to daily steps: timesteps.daily_steps)")

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
        if self.api.get_value(self.address("X", self.model, "")).shape != np.shape(heads):
            raise ValueError("saved heads don't match this model's nodes")
        self.set("X", "", np.asarray(heads, dtype=np.float64))

    def head(self) -> np.ndarray:
        """Heads, m, per reduced node."""
        return self.get("X", self.model, "")

    @property
    def end(self) -> float:
        """The simulation's end, days."""
        return self.api.get_end_time() * self.days

    def exchanges(self) -> list[tuple[str, str, int]]:
        """The GWF-GWF exchanges joining the coupled model to another groundwater model of the simulation, from
        mfsim.nam: (libmf6's name for it, such as "GWF-GWF_1"; the other model, lower case; +1 where the coupled model
        is its first model, -1 where it is its second). libmf6 names an exchange by its type and its place in the
        exchanges block."""
        from .mover import _blocks
        out = []
        for name, _, lines in _blocks(self.ws / "mfsim.nam", self.ws):
            if name == "EXCHANGES":
                for n, words in enumerate(lines, start=1):
                    kind = words[0].upper().replace("6", "")
                    a, b = words[2].upper(), words[3].upper()
                    if kind == "GWF-GWF" and self.model in (a, b):
                        out.append((f"{kind}_{n}", (b if a == self.model else a).lower(), 1 if a == self.model else -1))
        return out

    def exchange_flow(self, name: str, sign: int) -> np.ndarray:
        """A GWF-GWF exchange's flow for the last day across each of its connections, m3 d-1, positive into the
        coupled model (libmf6's SIMVALS are positive into the exchange's first model)."""
        return sign * self.api.get_value(f"{name}/SIMVALS").astype(float) * self.factor(Q)

    @property
    def period(self) -> int:
        """The stress period of the day prepared (1-based)."""
        return int(self.api.get_value("TDIS/KPER")[0])

    @property
    def time(self) -> float:
        """The simulation's time, days."""
        return self.api.get_current_time() * self.days

    def finalize(self):
        self.api.finalize()


class Models:
    """Several groundwater models of one simulation coupled together (``[modflow6] models``; docs/design.md, "Several
    groundwater models coupled together"): a view of each, and their cells numbered in one sequence, each model's
    after the one before, so a map, the recharge and the wells take a cell of any of them. A package named in more
    than one model is the same name in each (the recharge package): its entries are each model's in turn. What runs
    the simulation (prepare, solve, time) is the first model's, the one libmf6 was started for.

    Links that are one model's (the stream exchange, a lake, rivers) take that model's view (``view_of``)."""

    def __init__(self, first: Modflow6, models: list[str]):
        names = [m.upper() for m in models]
        if names[0] != first.model:
            raise ValueError(f"[modflow6] models: the first is {models[0]}, the model started is {first.model.lower()}")
        if len(set(names)) != len(names):
            raise ValueError(f"[modflow6] models: {models} names a model twice")
        self.first = first
        self.views = [first] + [first.view(m) for m in names[1:]]
        self.offset = np.concatenate([[0], np.cumsum([v.nodes for v in self.views])]).astype(int)
        self.cell_area = np.concatenate([v.cell_area for v in self.views])
        self.nodes = int(self.offset[-1])
        self.model = first.model

    def __getattr__(self, name):  # the simulation's own: api, ws, units, prepare_day, solve_day, end, ...
        if name == "first":
            raise AttributeError(name)
        return getattr(self.first, name)

    @property
    def names(self) -> list[str]:
        return [v.model for v in self.views]

    def view_of(self, model: str | None) -> Modflow6:
        if model is None:
            return self.first
        for v in self.views:
            if v.model == model.upper():
                return v
        raise ValueError(f"model {model!r} isn't one of the coupled models {[n.lower() for n in self.names]}")

    def model_of(self, node: int) -> tuple[Modflow6, int]:
        """The model a node of the sequence is in, and its node there."""
        i = int(np.searchsorted(self.offset, node, side="right")) - 1
        return self.views[i], int(node - self.offset[i])

    def node(self, *cell: int, model: str | None = None, node: int | None = None) -> int:
        v = self.view_of(model)
        return int(self.offset[self.views.index(v)]) + v.node(*cell, node=node)

    def _parts(self, package: str):
        return [(v, int(self.offset[i])) for i, v in enumerate(self.views) if v.has_package(package)]

    def package_nodes(self, package: str) -> np.ndarray:
        parts = self._parts(package)
        if not parts:
            raise ValueError(f"no coupled model has a package {package!r}")
        return np.concatenate([np.where(n >= 0, n + off, n) for v, off in parts
                               for n in [v.package_nodes(package)]])

    def nbound(self, package: str) -> int:
        return sum(v.nbound(package) for v, _ in self._parts(package))

    def has_package(self, package: str) -> bool:
        return bool(self._parts(package))

    def aux_multiplier(self, package: str) -> bool:
        return any(v.aux_multiplier(package) for v, _ in self._parts(package))

    def entry(self, package: str, node: int) -> int:
        idx = np.flatnonzero(self.package_nodes(package) == node)
        if len(idx) != 1:
            raise ValueError(f"{package}: {len(idx)} entries at node {node}")
        return int(idx[0])

    def _entry_of(self, package: str, entry: int) -> tuple[Modflow6, int]:
        for v, _ in self._parts(package):
            n = v.nbound(package)
            if entry < n:
                return v, entry
            entry -= n
        raise IndexError(f"{package}: no entry {entry}")

    def set_recharge(self, package: str, rate_per_node: np.ndarray):
        for v, off in self._parts(package):
            v.set_recharge(package, rate_per_node[off:off + v.nodes])

    def set_well(self, package: str, entry: int, q: float):
        v, i = self._entry_of(package, entry)
        v.set_well(package, i, q)

    def simvals(self, package: str) -> np.ndarray:
        return np.concatenate([v.simvals(package) for v, _ in self._parts(package)])

    def head(self) -> np.ndarray:
        return np.concatenate([v.head() for v in self.views])

    def set_heads(self, heads: np.ndarray):
        heads = np.asarray(heads, float)
        if heads.shape != (self.nodes,):
            raise ValueError("saved heads don't match the coupled models' nodes")
        for i, v in enumerate(self.views):
            v.set_heads(heads[self.offset[i]:self.offset[i + 1]])


def remaining_simulation(src: Path, dst: Path, day: int, time_units: str | None = None):
    """A copy of a MODFLOW 6 simulation covering only the days from ``day`` on, for a restart (``src`` and ``dst``
    the same folder: changed in place). Its steps must be one day each (what the driver's builders write, or any
    model once ``timesteps.daily_steps`` has rewritten it). The periods already over are left out, the one ``day``
    falls in is shortened, every input file's period blocks are renumbered to match and time series moved back
    (``timesteps.remaining_days``); the files are changed as text, so nothing else in them is rewritten.
    ``time_units``: the model's, where its TDIS leaves them undefined."""
    import shutil

    from .timesteps import remaining_days

    src, dst = Path(src), Path(dst)
    if src.resolve() != dst.resolve():  # or in place, in a run's copy
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(src, dst)
    remaining_days(dst, day, time_units)
    return dst
