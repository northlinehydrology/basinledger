"""Phased growth: a water user whose size follows a dated schedule (docs/design.md, "Scenarios of water users").

    phases = [{from = "2027-01-01", mw = 60}, {from = "2030-01-01", mw = 300, config = {design = "hybrid"}}]

Each phase starts on its date and lasts until the next; before the first the user is off. A phase may give ``mw``
(CoolingWaterBmi's ``it_capacity_mw``), ``config`` (any of the component's config keys for the phase) and
``multiplier`` (a factor on everything the user takes, returns and uses; default 1). A phase without ``mw`` or
``config`` keeps the config before it.

``PhasedUser`` is a BMI component holding one component per distinct config, each built from the user's config file
with the phase's keys put in. All are stepped together from the run's start with the same inputs, so each phase is
exactly the component at that size, whatever state or calendar it keeps; each is dropped once its last phase has
passed. The active phase gets the withdrawal limit (divided by its multiplier) and dry mode; the others run without
them. Outputs in a volume, volume flux, mass flux or power are the active component's times its multiplier;
a running total of time (units h or s) counts the active phase's increase; anything else (a temperature) is the
active component's own.
"""

from __future__ import annotations

import copy
import datetime as dt
import json
import math
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np
from bmipy import Bmi

from . import units as U
from .users import CONTROLS

SIZE_KEY = "it_capacity_mw"  # CoolingWaterBmi's IT capacity, what a phase's ``mw`` sets
SCALED = {**U.VOLUME_FLUX, **U.VOLUME, **U.MASS_FLUX, **U.POWER}  # units the multiplier applies to
RUNNING = U.TIME  # units of a running total of time


@dataclass(frozen=True)
class Phase:
    start: dt.date | None          # None: from the run's start (a scale's "before" phase)
    config: tuple = ()             # sorted (key, value) overrides of the component's config
    multiplier: float = 1.0
    off: bool = False              # the user does nothing (before a schedule's first phase)

    @property
    def overrides(self) -> dict:
        return {k: json.loads(v) for k, v in self.config}


def parse(where: str, phases) -> list[Phase]:
    """``phases`` as configured: a list of {from, mw, config, multiplier}, dates increasing."""
    if not isinstance(phases, list) or not phases:
        raise ValueError(f"{where}: phases is a list of {{from = \"YYYY-MM-DD\", mw, config, multiplier}}")
    out: list[Phase] = []
    current: dict = {}
    for p in phases:
        if not isinstance(p, dict) or "from" not in p or set(p) - {"from", "mw", "config", "multiplier"}:
            raise ValueError(f"{where}: each phase has from (a date) and may have mw, config and multiplier; got {p!r}")
        try:
            start = dt.date.fromisoformat(str(p["from"]))
        except ValueError:
            raise ValueError(f"{where}: phase from = {p['from']!r} is not a date (YYYY-MM-DD)") from None
        if out and start <= out[-1].start:
            raise ValueError(f"{where}: phases' dates must increase ({out[-1].start}, then {start})")
        cfg = p.get("config", None)
        if cfg is not None and not isinstance(cfg, dict):
            raise ValueError(f"{where}: a phase's config is a table of the component's config keys; got {cfg!r}")
        if "mw" in p or cfg is not None:
            current = dict(cfg or {})
            if "mw" in p:
                if not _number(p["mw"]) or p["mw"] < 0:
                    raise ValueError(f"{where}: a phase's mw is a capacity, MW >= 0; got {p['mw']!r}")
                if SIZE_KEY in current:
                    raise ValueError(f"{where}: a phase gives mw or config.{SIZE_KEY}, not both")
                current[SIZE_KEY] = float(p["mw"])
        m = p.get("multiplier", 1.0)
        if not _number(m) or m < 0:
            raise ValueError(f"{where}: a phase's multiplier is a number >= 0; got {m!r}")
        out.append(Phase(start, _frozen(current), float(m)))
    return out


def scale_schedule(where: str, multiplier=None, start: dt.date | None = None, phases=None) -> list[Phase]:
    """A ``scale`` (or ``remove``) on an existing user: 1 before it, then its multiplier from ``start``, or its own
    multiplier-only phases."""
    if phases is not None:
        if multiplier is not None:
            raise ValueError(f"{where}: a scale takes multiplier (with from) or phases, not both")
        ps = parse(where, phases)
        if any(p.config for p in ps):
            raise ValueError(f"{where}: a scale's phases take multipliers only; a change of size is the user's own "
                             "phases")
    else:
        if not _number(multiplier) or multiplier < 0:
            raise ValueError(f"{where}: a scale needs multiplier, a number >= 0 (or phases); got {multiplier!r}")
        ps = [Phase(start, (), float(multiplier))]
    return [Phase(None)] + ps


def combine(own: list[Phase] | None, scale: list[Phase]) -> list[Phase]:
    """A user's own phases (None: always on, as configured) with a scale's schedule: at each date of either, the own
    phase's config and the product of the two multipliers."""
    own = own or [Phase(None)]
    if own[0].start is not None:  # off before its first phase
        own = [Phase(None, (), 0.0, off=True)] + own
    dates = sorted({p.start for p in own + scale if p.start is not None})
    points = [None] + dates

    def at(sched, d):
        cur = sched[0]
        for p in sched:
            if p.start is None or (d is not None and p.start <= d):
                cur = p
        return cur

    out = []
    for d in points:
        a, b = at(own, d), at(scale, d)
        out.append(Phase(d, a.config, a.multiplier * b.multiplier, a.off))
    return out


def _number(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def _frozen(d: dict) -> tuple:
    return tuple(sorted((k, json.dumps(v, sort_keys=True)) for k, v in d.items()))


def as_config(phases: list[Phase]) -> list[dict]:
    """The phases back as config tables (for fingerprints and the run's record)."""
    return [{"from": None if p.start is None else p.start.isoformat(), "config": p.overrides,
             "multiplier": p.multiplier, "off": p.off} for p in phases]


# ---------- config files for a phase ----------

def read_config(path: Path) -> dict:
    """A component's config file as a table: TOML, or YAML/JSON (pyyaml reads both)."""
    if path.suffix == ".toml":
        with open(path, "rb") as f:
            return tomllib.load(f)
    if path.suffix in (".json",):
        return json.loads(path.read_text())
    if path.suffix in (".yaml", ".yml"):
        import yaml
        return yaml.safe_load(path.read_text()) or {}
    raise ValueError(f"{path}: a phased user's config is read to put a phase's keys in: TOML, YAML or JSON")


def write_config(path: Path, cfg: dict, like: Path) -> Path:
    """A phase's config, in the format of the file it came from (JSON for YAML, which reads it), relative paths in it
    made absolute against that file's folder by the caller."""
    if like.suffix == ".toml":
        path = path.with_suffix(".toml")
        path.write_text("".join(f"{_toml_key(k)} = {_toml(v)}\n" for k, v in cfg.items()))
    else:
        path = path.with_suffix(like.suffix)
        path.write_text(json.dumps(cfg, indent=1))
    return path


def _toml_key(k: str) -> str:
    return k if k and all(c.isalnum() or c in "_-" for c in k) else json.dumps(k)


def _toml(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return "inf" if v == math.inf else "-inf" if v == -math.inf else "nan" if v != v else repr(v)
    if isinstance(v, str):
        return json.dumps(v)
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(_toml(x) for x in v) + "]"
    if isinstance(v, dict):
        return "{" + ", ".join(f"{_toml_key(k)} = {_toml(x)}" for k, x in v.items()) + "}"
    if isinstance(v, (dt.date, dt.datetime)):
        return v.isoformat()
    raise ValueError(f"a config value TOML can't hold: {v!r}")


# ---------- the component ----------

@dataclass
class _Part:
    bmi: object
    last: float                     # the time (component units) its last phase ends; inf: never
    outputs: dict = field(default_factory=dict)


class PhasedUser(Bmi):
    """A water user stepped through its phases (module docstring). ``make()`` returns a new, uninitialized component;
    ``start`` is the run's first day, which the phases' dates count from; ``folder`` takes the phases' config files;
    ``limit`` is the component's withdrawal limit input (``variables.limit``, or the common name)."""

    def __init__(self, make: Callable[[], object], phases: list[Phase], start: dt.date, folder: Path,
                 limit: str | None = None, dry_mode: str | None = None):
        self._make = make
        self._phases = list(phases)
        self._start = start
        self._folder = Path(folder)
        self._limit = limit or CONTROLS["limit"]
        self._dry = dry_mode or CONTROLS["dry_mode"]

    # ---------- control ----------
    def initialize(self, config_file=None):
        base_path = Path(config_file) if config_file else None
        base = read_config(base_path) if base_path and any(p.config for p in self._phases) else None
        by_config: dict[tuple, int] = {}
        self._parts: list[_Part] = []
        self._schedule = []  # (start time in component units, part index or None when off, multiplier)
        files = 0
        for p in self._phases:
            if p.config not in by_config:
                comp = self._make()
                if p.config:
                    if base_path is None:
                        raise ValueError("a phase that changes the component's config needs the user's config file")
                    cfg = copy.deepcopy(base)
                    cfg.update(p.overrides)
                    _absolute(cfg, base_path.parent)
                    self._folder.mkdir(parents=True, exist_ok=True)
                    files += 1
                    path = write_config(self._folder / f"phase{files}", cfg, base_path)
                    comp.initialize(str(path))
                else:
                    comp.initialize(str(base_path) if base_path else None)
                by_config[p.config] = len(self._parts)
                self._parts.append(_Part(comp, math.inf))
        first = self._parts[0].bmi
        self._seconds = U.time_factor(first.get_time_units())
        self._dt = first.get_time_step()
        for part in self._parts[1:]:
            if part.bmi.get_time_step() != self._dt or part.bmi.get_time_units() != first.get_time_units():
                raise ValueError("a phased user's components step differently (each phase's config keeps the time "
                                 "step)")
        for p in self._phases:
            t = 0.0 if p.start is None else (p.start - self._start).days * U.DAY / self._seconds
            self._schedule.append((t, None if p.off else by_config[p.config], p.multiplier))
        for i, part in enumerate(self._parts):  # when each part's last phase ends
            later = [k for k, (_, j, _) in enumerate(self._schedule) if j == i]
            k = later[-1] + 1
            part.last = self._schedule[k][0] if k < len(self._schedule) else math.inf
        self._t = first.get_current_time()
        self._t0 = self._t
        self._out = {v: np.zeros(self._size(first, v), dtype=np.float64) for v in first.get_output_var_names()}
        self._kind = {}
        for v in self._out:
            u = first.get_var_units(v).strip()
            self._kind[v] = "scaled" if u in SCALED else "running" if u in RUNNING else "own"
        self._inputs: dict[str, np.ndarray] = {}
        self._ins = set(first.get_input_var_names())

    @staticmethod
    def _size(bmi, var):
        return bmi.get_grid_size(bmi.get_var_grid(var))

    def _active(self, t: float):
        """(part index or None, multiplier) of the phase in effect at component time ``t``."""
        cur = (None, 0.0)
        eps = 1e-9 * max(1.0, abs(t))
        for start, part, m in self._schedule:
            if start <= t + eps:
                cur = (part, m)
        return cur

    def update(self):
        t = self._t
        active, m = self._active(t)
        live = [i for i, p in enumerate(self._parts) if p.bmi is not None]
        before = {}
        for i in live:
            b = self._parts[i].bmi
            for v, x in self._inputs.items():
                if v == self._limit:
                    lim = float(x.reshape(-1)[0])
                    x = np.full_like(x, np.inf if i != active or m == 0 else lim / m)
                elif v == self._dry and i != active:
                    x = np.zeros_like(x)
                b.set_value(v, x)
            if i == active:
                before = {v: self._get(b, v) for v, k in self._kind.items() if k == "running"}
            b.update()
        times = {self._parts[i].bmi.get_current_time() for i in live}
        if len(times) != 1:
            raise RuntimeError("a phased user's components fell out of step")
        self._t = times.pop()
        for v, k in self._kind.items():
            if active is None:
                if k != "running":
                    self._out[v][:] = 0.0
                continue
            x = self._get(self._parts[active].bmi, v)
            if k == "scaled":
                self._out[v][:] = x * m
            elif k == "running":
                self._out[v] += x - before[v]
            else:
                self._out[v][:] = x
        for i in live:  # a part whose last phase has passed is done
            p = self._parts[i]
            if p.last <= self._t + 1e-9 * max(1.0, abs(self._t)) and len(live) > 1:
                p.bmi.finalize()
                p.bmi = None

    def _get(self, bmi, var) -> np.ndarray:
        buf = np.zeros(self._out[var].shape, dtype=np.float64)
        bmi.get_value(var, buf)
        return buf

    def update_until(self, time):
        while self._t + self._dt <= time + 1e-9 * max(1.0, abs(time)):
            self.update()

    def finalize(self):
        for p in self._parts:
            if p.bmi is not None:
                p.bmi.finalize()
                p.bmi = None

    # ---------- information: the first component's ----------
    def _any(self):
        for p in self._parts:
            if p.bmi is not None:
                return p.bmi
        raise RuntimeError("a phased user after finalize")

    def get_component_name(self):
        return f"{self._any().get_component_name()} (in phases)"

    def get_input_item_count(self):
        return len(self._ins)

    def get_output_item_count(self):
        return len(self._out)

    def get_input_var_names(self):
        return self._any().get_input_var_names()

    def get_output_var_names(self):
        return tuple(self._out)

    def get_var_grid(self, name):
        return self._any().get_var_grid(name)

    def get_var_type(self, name):
        return "float64" if name in self._out else self._any().get_var_type(name)

    def get_var_units(self, name):
        return self._any().get_var_units(name)

    def get_var_itemsize(self, name):
        return 8 if name in self._out else self._any().get_var_itemsize(name)

    def get_var_nbytes(self, name):
        return self._out[name].nbytes if name in self._out else self._any().get_var_nbytes(name)

    def get_var_location(self, name):
        return self._any().get_var_location(name)

    def get_current_time(self):
        return self._t

    def get_start_time(self):
        return self._t0

    def get_end_time(self):
        return self._any().get_end_time()

    def get_time_units(self):
        return self._any().get_time_units()

    def get_time_step(self):
        return self._dt

    def get_value(self, name, dest):
        if name in self._out:
            dest[:] = self._out[name]
            return dest
        return self._any().get_value(name, dest)

    def get_value_ptr(self, name):
        return self._out[name]

    def get_value_at_indices(self, name, dest, inds):
        dest[:] = self._out[name].take(inds)
        return dest

    def set_value(self, name, src):
        if name not in self._ins:
            raise KeyError(f"{name} is not an input variable")
        self._inputs[name] = np.array(src, dtype=np.float64, copy=True)

    def set_value_at_indices(self, name, inds, src):
        cur = self._inputs.get(name)
        if cur is None:
            cur = np.zeros(self._size(self._any(), name))
        cur = cur.copy()
        cur[np.asarray(inds)] = src
        self.set_value(name, cur)

    def get_grid_rank(self, grid):
        return self._any().get_grid_rank(grid)

    def get_grid_size(self, grid):
        return self._any().get_grid_size(grid)

    def get_grid_type(self, grid):
        return self._any().get_grid_type(grid)

    def get_grid_node_count(self, grid):
        return self._any().get_grid_node_count(grid)

    def get_grid_edge_count(self, grid):
        return self._any().get_grid_edge_count(grid)

    def get_grid_face_count(self, grid):
        return self._any().get_grid_face_count(grid)

    def _passed(name):
        def f(self, *a):
            return getattr(self._any(), name)(*a)
        return f

    get_grid_shape = _passed("get_grid_shape")
    get_grid_spacing = _passed("get_grid_spacing")
    get_grid_origin = _passed("get_grid_origin")
    get_grid_x = _passed("get_grid_x")
    get_grid_y = _passed("get_grid_y")
    get_grid_z = _passed("get_grid_z")
    get_grid_edge_nodes = _passed("get_grid_edge_nodes")
    get_grid_face_edges = _passed("get_grid_face_edges")
    get_grid_face_nodes = _passed("get_grid_face_nodes")
    get_grid_nodes_per_face = _passed("get_grid_nodes_per_face")
    del _passed


def _absolute(cfg: dict, base: Path):
    """Relative file paths in a phase's config (a ReportedUser's series, a forcing file) made absolute against the
    original config's folder, since the phase's file is written elsewhere."""
    for k, v in cfg.items():
        if isinstance(v, str) and ("/" in v or "." in v) and not Path(v).is_absolute() and (base / v).exists():
            cfg[k] = str((base / v).resolve())
