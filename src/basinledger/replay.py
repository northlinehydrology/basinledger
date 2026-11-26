"""Replay: record a coupled SWAT+ / MODFLOW 6 run's flow side once, then run the salt again from the record, without
SWAT+ or MODFLOW 6's flow model (docs/design.md, Replay). Built Nov 22, 2026 for the salt calibration of long, heavy
runs: flow never depends on transport, and salt sends nothing back to SWAT+, so one recorded flow run can be replayed
many times, cheaply and in parallel, with only salt's settings changed.

**Recording** (``[replay] record = "<folder>"`` on a normal run; ``Recorder``): after each day, what transport and the
salt balances read from the flow side:

- MODFLOW 6's flow, as GWT reads it through the GWF-GWT exchange: the face flows (FLOWJA, the upper triangle; the
  lower is its negative, checked every day), each cell's saturation, specific discharge and storage terms, the cells
  MODFLOW 6 has made dry, and every flow package's flow by entry (the coupled rivers, drains and lakes, recharge, the
  soil link, wells, and the model's own), with each package's entries checked to stay where they were;
- SWAT+'s day as the salt balances read it (the ``hru_wb`` columns, percolation, channel and reservoir volumes,
  flows, rain and seepage, point sources, withdrawals, fertilizer), the soil link's plan by map row (gwflow's rule),
  the users' ledger lines and controls, ``recharge:sent``, and the flow, base flow and head series;
- once: the run's identity (``identity``), the starting state, the recharge and soil maps, the exchange's entries,
  the users, and the SWAT+ files the salt reads at start-up (routing, HRU data, the stations' weather).

MODFLOW 6's arrays are float32 by default (``[replay] precision``), SWAT+'s are float32 as SWAT+ keeps them (exact).
Each chunk of days is stored per array: XOR with the day before, bytes shuffled, then zlib.

**Replaying** (``Basin(cfg, replay="<folder>")``, ``basinledger run --replay``, ``basinledger forward --replay``;
``ReplayBasin``): the same MODFLOW 6 simulation is built (its flow model's input checked against the record's), with
the flow model left unsolved: each day, after MODFLOW 6 prepares it, the flow model's memory is set from the record
and only the transport solutions are solved. The salt balances read SWAT+'s day from the record. Users who take in
the aquifer's water are run again (their water checked against the record); the others' lines come from the record.
Anything that would change the flow is refused (``identity``).
"""

from __future__ import annotations

import contextlib
import dataclasses
import datetime as dt
import hashlib
import json
import shutil
import time
import warnings
import zlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from . import config as C

FORMAT = 1
DHDRY = -1e30  # MODFLOW 6's head in a dry cell
CHUNK_BYTES = 256 * 2 ** 20  # raw bytes of a chunk of days, about (its first day stored whole, the rest XORed)
LEVEL = 6  # zlib


class ReplayRefused(ValueError):
    """A replay that would change the flow, or a record that can't serve it."""


# ---------- identity ----------

FLOW_FIELDS_LEFT_OUT = {  # ExampleConfig's fields that don't change the flow (or that the flow model's files carry)
    "path", "name", "description", "run_folder", "days", "transport", "salt", "salt_conductance", "outputs",
    "calibration_record", "calibration_parameters", "replay_record", "replay_precision", "replay_mantissa_bits",
    "replay_specific_discharge", "mf6_builder_args", "mf6_workspace",
}
SALT_PARAMETER_MODELS = ("salt",)


def _hash_bytes(*parts: bytes) -> str:
    h = hashlib.sha256()
    for p in parts:
        h.update(p)
    return h.hexdigest()


def _hash_path(p: Path) -> str:
    """A file's bytes, or a folder's files (names and bytes) in order."""
    h = hashlib.sha256()
    files = [p] if p.is_file() else sorted(x for x in p.rglob("*") if x.is_file())
    for f in files:
        h.update(str(f.relative_to(p) if f != p else f.name).encode())
        with open(f, "rb") as fh:
            for block in iter(lambda: fh.read(1 << 20), b""):
                h.update(block)
    return h.hexdigest()


def _canonical(x, folder: Path, files: dict):
    """A config value as JSON, every path in it (or a string naming a file beside the config) by its content."""
    if dataclasses.is_dataclass(x):
        return {f.name: _canonical(getattr(x, f.name), folder, files) for f in dataclasses.fields(x)}
    if isinstance(x, Path):
        if x.exists() and not (folder == x or folder.is_relative_to(x)):
            files.setdefault(str(x), _hash_path(x))
        return str(x)
    if isinstance(x, dict):
        return {str(k): _canonical(v, folder, files) for k, v in sorted(x.items(), key=lambda kv: str(kv[0]))}
    if isinstance(x, (list, tuple)):
        return [_canonical(v, folder, files) for v in x]
    if isinstance(x, str) and x and not x.startswith("all:"):
        with contextlib.suppress(OSError, ValueError):
            q = (folder / x).resolve()
            if q.is_file():
                files.setdefault(str(q), _hash_path(q))
    if isinstance(x, float):
        return repr(x)
    return x


def identity(cfg: C.ExampleConfig, parameters: list, users: bool, heads) -> dict:
    """What the flow of a run depends on outside MODFLOW 6's flow model's own files (``gwf_inputs``): the config's
    flow side (every field but the salt's, the outputs and the run's place and length), the content of every file it
    names (the SWAT+ project, maps, users' configs, weather, the builder's module), the parameters of the models
    (not salt's), whether the users run, and the starting heads given. Each part hashed."""
    files: dict[str, str] = {}
    flow = {f.name: _canonical(getattr(cfg, f.name), cfg.folder, files)
            for f in dataclasses.fields(cfg) if f.name not in FLOW_FIELDS_LEFT_OUT}
    if cfg.mf6_builder:  # the builder's code: its arguments are checked through the files it writes
        mod = cfg.mf6_builder.partition(":")[0]
        if (cfg.folder / f"{mod}.py").is_file():
            files[str(cfg.folder / f"{mod}.py")] = _hash_path(cfg.folder / f"{mod}.py")
    params = sorted((p.name, p.model, p.target, repr(float(p.value))) for p in parameters
                    if p.model not in SALT_PARAMETER_MODELS)
    return {
        "config": _hash_bytes(json.dumps(flow, sort_keys=True, default=str).encode()),
        "files": files,
        "parameters": _hash_bytes(json.dumps(params).encode()),
        "users": bool(users),
        "heads": None if heads is None else _hash_bytes(np.asarray(heads, np.float64).tobytes()),
    }


def _file_lines(path: Path, ws: Path, h, seen: set):
    """A MODFLOW 6 input file into the hash: its words, comments left out, and the files it reads (OPEN/CLOSE, FILEIN)
    in place (a binary array by its bytes)."""
    if path in seen or not path.is_file():
        return
    seen.add(path)
    h.update(f"<{path.name}>".encode())
    for line in path.read_text(errors="replace").splitlines():
        words = line.split("#", 1)[0].split("!", 1)[0].split()
        if not words:
            continue
        up = [w.upper() for w in words]
        if up[0] == "OPEN/CLOSE":
            f = ws / words[1].strip("'\"")
            if "(BINARY)" in up or "BINARY" in up:
                h.update(_hash_path(f).encode()) if f.is_file() else None
            else:
                _file_lines(f, ws, h, seen)
            h.update(" ".join(words[2:]).encode())
            continue
        if "FILEIN" in up:
            i = up.index("FILEIN")
            if i + 1 < len(words):
                _file_lines(ws / words[i + 1].strip("'\""), ws, h, seen)
        h.update(" ".join(words).encode() + b"\n")


def gwf_inputs(ws: Path, model: str) -> str:
    """A hash of everything MODFLOW 6's flow model ``model`` reads: the simulation's options and time discretization,
    the flow model's name file and every package file (with the files they read), its solution's settings, and any
    GWF-GWF exchange. The transport models, their solutions and the GWF-GWT exchanges are left out, so a replay may
    change them."""
    from .coupling.mover import _blocks
    ws = Path(ws)
    h = hashlib.sha256()
    seen: set = set()
    model = model.upper()
    blocks = _blocks(ws / "mfsim.nam", ws)
    gwt = {w[2].upper() for b, _, lines in blocks if b == "MODELS" for w in lines if w[0].upper().startswith("GWT")}
    for name, extra, lines in blocks:
        if name == "OPTIONS":
            h.update(repr(lines).encode())
        elif name == "TIMING":
            for w in lines:
                _file_lines(ws / w[1].strip("'\""), ws, h, seen)
        elif name == "MODELS":
            for w in lines:
                if w[2].upper() == model:
                    _file_lines(ws / w[1].strip("'\""), ws, h, seen)
                    for b, _, plines in _blocks(ws / w[1].strip("'\""), ws):
                        if b == "PACKAGES":
                            for p in plines:
                                if p[0].upper() not in ("OC6", "OBS6"):  # outputs only
                                    _file_lines(ws / p[1].strip("'\""), ws, h, seen)
        elif name == "EXCHANGES":
            for w in lines:
                if w[0].upper().startswith("GWF6-GWF6"):
                    h.update(" ".join(w).encode())
                    _file_lines(ws / w[1].strip("'\""), ws, h, seen)
        elif name == "SOLUTIONGROUP":
            for w in lines:
                if model in (x.upper() for x in w[2:]) and not {x.upper() for x in w[2:]} & gwt:
                    _file_lines(ws / w[1].strip("'\""), ws, h, seen)
    return h.hexdigest()


def build_modflow6(cfg: C.ExampleConfig, ws: Path, days: int, start: dt.date, prepare=None):
    """The run's copy of MODFLOW 6's files, as ``Basin`` makes it for a fresh run (driver.py): built or copied, every
    step one day, users' added wells, storage, an added recharge package, and ``prepare("modflow6", ws)`` (the
    parameters). Notes are left out (the recorded run gave them)."""
    from .coupling import existing as EX
    from .coupling.timesteps import daily_steps
    from .driver import Basin, _example_function
    ws = Path(ws)
    if ws.exists():
        shutil.rmtree(ws)
    if cfg.mf6_builder:
        _example_function(cfg, cfg.mf6_builder, "build")(ws, days, start, **cfg.mf6_builder_args)
    else:
        shutil.copytree(cfg.mf6_workspace, ws)
    daily_steps(ws, cfg.mf6_time_units)
    by_model: dict[str, list] = {}
    for w in cfg.wells:
        if w.add:
            by_model.setdefault((w.model or cfg.mf6_model).lower(), []).append((w.name, w.package, w.cell))
    for model, added in by_model.items():
        EX.add_wells(ws, model, added)
    if cfg.mf6_storage is not None:
        st = cfg.mf6_storage
        EX.add_storage(ws, cfg.mf6_model, st.get("sy"), st.get("ss"), int(st.get("iconvert", 1)), cfg.mf6_length_units)
    if cfg.recharge_add:
        for model, cells in Basin._map_cells(cfg).items():
            EX.add_recharge(ws, model, cfg.recharge_package, cells)
    if prepare is not None:
        prepare("modflow6", ws)
    return ws


def _strip_outputs(ws: Path, model: str):
    """The flow model's output control and observations left out of a replay's copy: it isn't solved, so what it
    would write is meaningless."""
    from .coupling.mover import _blocks
    for _, _, lines in _blocks(ws / "mfsim.nam", ws):
        for w in lines:
            if len(w) >= 3 and w[0].upper().startswith("GWF") and w[2].upper() == model.upper():
                nam = ws / w[1].strip("'\"")
                keep = [ln for ln in nam.read_text().splitlines()
                        if not (ln.split()[:1] and ln.split()[0].upper() in ("OC6", "OBS6"))]
                nam.write_text("\n".join(keep) + "\n")
                return


# ---------- storage: chunks of days ----------

def _encode(a: np.ndarray, xor: bool) -> bytes:
    a = np.ascontiguousarray(a)
    w = a.dtype.itemsize
    if xor and a.ndim == 2 and a.shape[0] > 1 and a.dtype.kind in "fiu":
        u = a.view({1: np.uint8, 2: np.uint16, 4: np.uint32, 8: np.uint64}[w])
        d = u.copy()
        d[1:] ^= u[:-1]
        a = d
    b = np.ascontiguousarray(a.view(np.uint8).reshape(-1, w).T) if w > 1 else a.view(np.uint8)
    return zlib.compress(b.tobytes(), LEVEL)


def _decode(raw: bytes, dtype: str, shape: list, xor: bool) -> np.ndarray:
    dt_ = np.dtype(dtype)
    w = dt_.itemsize
    b = np.frombuffer(zlib.decompress(raw), np.uint8)
    if w > 1:
        b = np.ascontiguousarray(b.reshape(w, -1).T)
    a = b.view(dt_).reshape(shape) if shape else b.view(dt_)
    if xor and a.ndim == 2 and a.shape[0] > 1:
        u = a.view({1: np.uint8, 2: np.uint16, 4: np.uint32, 8: np.uint64}[w])
        a = np.bitwise_xor.accumulate(u, axis=0).view(dt_)
    return np.array(a)


def write_chunk(path: Path, days: list[dict]) -> int:
    """Days (each {name: 1-D array} and ``"json"``: small values) into one file; returns its bytes. An array the same
    length every day is stacked (and XORed day to day); one that changes length is stored end to end with its
    lengths."""
    names = sorted({k for d in days for k in d if k != "json"})
    out, index = {}, {"days": len(days), "json": [d.get("json", {}) for d in days], "arrays": {}}
    for k in names:
        parts = [np.asarray(d[k]) for d in days]
        if len({p.shape for p in parts}) == 1:
            a = np.stack(parts)
            out[k] = np.frombuffer(_encode(a, True), np.uint8)
            index["arrays"][k] = {"dtype": a.dtype.str, "shape": list(a.shape), "xor": True}
        else:
            a = np.concatenate([p.ravel() for p in parts])
            out[k] = np.frombuffer(_encode(a, False), np.uint8)
            index["arrays"][k] = {"dtype": a.dtype.str, "shape": list(a.shape), "xor": False,
                                  "lengths": [int(p.size) for p in parts]}
    out["__index__"] = np.frombuffer(json.dumps(index).encode(), np.uint8)
    tmp = path.with_suffix(".tmp.npz")
    np.savez(tmp, **out)
    tmp.replace(path)
    return path.stat().st_size


def read_chunk(path: Path) -> list[dict]:
    with np.load(path) as z:
        index = json.loads(bytes(z["__index__"]).decode())
        arrays = {k: _decode(bytes(z[k]), v["dtype"], v["shape"], v["xor"]) for k, v in index["arrays"].items()}
    days = [{"json": j} for j in index["json"]]
    for k, v in index["arrays"].items():
        a = arrays[k]
        if "lengths" in v:
            ends = np.cumsum(v["lengths"])
            for d, (s, e) in zip(days, zip(np.concatenate([[0], ends[:-1]]), ends)):
                d[k] = a[s:e]
        else:
            for i, d in enumerate(days):
                d[k] = a[i]
    return days


# ---------- what is read from the flow side ----------

def surface_variables(sw) -> tuple[list[str], list[str]]:
    """SWAT+'s variables the salt balances read each day, and those read once at start-up, that this run has."""
    from . import salt as SA
    from .coupling import swatplus as S
    daily = [*SA.HRU_READS.values(), SA.PERCOLATION, SA.CH_PRECIP, SA.CH_VOLUME, S.QOUT, SA.RES_VOLUME, SA.RES_OUT,
             SA.RES_SEEP, SA.RES_PRECIP, SA.POINT, SA.WITHDRAWN, *SA.FERTILIZER]
    start = [SA.SOIL_WATER, SA.WILTING, SA.CH_VOLUME, SA.RES_VOLUME]
    return [v for v in dict.fromkeys(daily) if sw.has(v)], [v for v in dict.fromkeys(start) if sw.has(v)]


SWATPLUS_FILES = ("file.cio", "hru-data.hru", "snow.sno", "weather-sta.cli", "rout_unit.def", "rout_unit.ele")


def copy_swatplus_files(src: Path, dst: Path):
    """The SWAT+ files the salt reads at start-up: the routing (``*.con``, routing units), HRU data, snow, and the
    stations' precipitation and temperature records (the annual road salt's days)."""
    dst.mkdir(parents=True, exist_ok=True)
    names = set(SWATPLUS_FILES) | {p.name for p in src.glob("*.con")}
    cli = src / "weather-sta.cli"
    if cli.exists():
        lines = cli.read_text().splitlines()
        head = lines[1].split()
        for ln in lines[2:]:
            row = dict(zip(head, ln.split()))
            for k in ("pcp", "tmp"):
                if row.get(k) not in (None, "sim", "null"):
                    names.add(row[k])
    for n in sorted(names):
        if (src / n).is_file():
            shutil.copy2(src / n, dst / n)


class Precision:
    """How MODFLOW 6's flows are stored: float64 (bit for bit), float32 (the default), or float32 with its mantissa
    rounded to ``bits`` (of 23) [A]: 16 bits keep each value within 2^-17 (7.6e-6) of the run's and store the
    Rappahannock grid's flows in 30% less. SWAT+'s values are always float32, as SWAT+ holds them (exact)."""

    def __init__(self, precision: str = "float32", bits: int | None = None):
        if precision not in ("float32", "float64"):
            raise ValueError(f"[replay] precision is float32 or float64, not {precision!r}")
        if bits is not None and (precision != "float32" or not 8 <= int(bits) <= 23):
            raise ValueError("[replay] mantissa_bits (8 to 23) rounds float32 values")
        self.precision, self.bits = precision, None if bits in (None, 23) else int(bits)
        self.exact = precision == "float64"

    def __call__(self, a) -> np.ndarray:
        a = np.asarray(a, np.float64)
        if self.exact:
            return a.copy()
        a = a.astype(np.float32)
        if self.bits is None:
            return a
        drop = 23 - self.bits
        u = a.view(np.uint32).astype(np.uint64)
        u = ((u + (1 << (drop - 1))) >> drop) << drop  # to nearest, ties up
        return u.astype(np.uint32).view(np.float32)


class _Gwf:
    """MODFLOW 6's flow model's arrays GWT reads, by memory address, in the model's own units (nothing converted: they
    go back where they came from)."""

    def __init__(self, mf6):
        from .coupling.existing import flow_packages
        self.api, self.m = mf6.api, mf6.model
        names = set(self.api.get_input_var_names()) | set(self.api.get_output_var_names())
        self.has = lambda v: f"{self.m}/{v}" in names
        ia = self.api.get_value(f"{self.m}/CON/IA").astype(np.int64) - 1
        ja = self.api.get_value(f"{self.m}/CON/JA").astype(np.int64) - 1
        isym = self.api.get_value(f"{self.m}/CON/ISYM").astype(np.int64) - 1
        self.nodes, self.nja = len(ia) - 1, len(ja)
        row = np.repeat(np.arange(self.nodes), np.diff(ia))
        self.upper = np.flatnonzero(ja > row)
        self.lower = isym[self.upper]
        self.diag = ia[:-1]
        self.structure = _hash_bytes(ia.tobytes(), ja.tobytes())
        self.spdis = self.has("NPF/SPDIS") and int(self.api.get_value(f"{self.m}/NPF/ICALCSPDIS")[0]) != 0
        self.storage = [v for v in ("STO/STRGSS", "STO/STRGSY") if self.has(v)]
        self.packages = flow_packages(mf6)
        self.movers = [p for p in self.packages if self.has(f"{p}/IMOVER")
                       and int(self.api.get_value(f"{self.m}/{p}/IMOVER")[0]) != 0]
        # a transport model correcting for the flow's imbalance reads FLOWJA's diagonal
        self.flow_error = any(n.endswith("/FMI/IFLOWERR") and int(self.api.get_value(n)[0]) != 0 for n in names)

    def get(self, v):
        return self.api.get_value(f"{self.m}/{v}")

    def set(self, v, x):
        self.api.set_value(f"{self.m}/{v}", np.asarray(x, np.float64))

    def nbound(self, p) -> int:
        return int(self.get(f"{p}/NBOUND")[0])

    def nodelist_crc(self, p) -> int:
        return zlib.crc32(self.get(f"{p}/NODELIST")[: self.nbound(p)].astype(np.int32).tobytes())

    def state(self, cast: Precision) -> dict:
        """The day's flow, as stored. Transport reads only the sum of the two storage terms (MST and FMI's old
        saturation: ``vold = vnew + (strgss + strgsy) dt``), so it's stored as one, except in float64."""
        f = self.get("FLOWJA")
        up = f[self.upper]
        if not np.array_equal(f[self.lower], -up):
            bad = float(np.max(np.abs(f[self.lower] + up)))
            raise ValueError(f"replay: MODFLOW 6's face flows aren't antisymmetric (by {bad:g}); the record keeps one "
                             "side of each face")
        out = {"flowja": cast(up), "sat": cast(self.get("NPF/SAT")),
               "dry": np.flatnonzero(self.get("X") == DHDRY).astype(np.int32)}
        if self.flow_error:
            out["flowja_diag"] = cast(f[self.diag])
        if self.spdis:
            out["spdis"] = cast(self.get("NPF/SPDIS").ravel())
        if cast.exact:
            for v in self.storage:
                out[v.split("/")[1].lower()] = cast(self.get(v))
        elif self.storage:
            out["storage"] = cast(sum(self.get(v) for v in self.storage))
        for p in self.packages:
            n = self.nbound(p)
            out[f"q:{p}"] = cast(self.get(f"{p}/SIMVALS")[:n])
            if p in self.movers:
                out[f"mvr:{p}"] = cast(self.get(f"{p}/SIMTOMVR")[:n])
        return out

    def structure_today(self) -> dict:
        return {p: [self.nbound(p), self.nodelist_crc(p)] for p in self.packages}

    def put(self, day: dict, x0: np.ndarray, before: bool = False):
        """Set the flow model's memory from a recorded day: all of it after MODFLOW 6 prepares the day, or, before a
        replay's first day (``before``), what GWT reads as it prepares (the dry cells, saturation, specific discharge
        and face flows of the day before)."""
        f = np.zeros(self.nja)
        up = np.asarray(day["flowja"], np.float64)
        f[self.upper] = up
        f[self.lower] = -up
        if "flowja_diag" in day:
            f[self.diag] = day["flowja_diag"]
        self.set("FLOWJA", f)
        self.set("NPF/SAT", day["sat"])
        x = np.array(x0, np.float64)
        x[np.asarray(day["dry"], int)] = DHDRY
        self.set("X", x)
        if self.spdis and "spdis" in day:
            self.set("NPF/SPDIS", np.asarray(day["spdis"], np.float64).reshape(self.nodes, -1))
        if before:
            return
        if "storage" in day:  # the two terms' sum in the first, zero in the other: the same to transport
            for i, v in enumerate(self.storage):
                self.set(v, day["storage"] if i == 0 else np.zeros(self.nodes))
        else:
            for v in self.storage:
                self.set(v, day[v.split("/")[1].lower()])
        for p in self.packages:
            for key, var in ((f"q:{p}", "SIMVALS"), (f"mvr:{p}", "SIMTOMVR")):
                if key in day:
                    cur = self.get(f"{p}/{var}")
                    cur[: len(day[key])] = day[key]
                    self.set(f"{p}/{var}", cur)


# ---------- recording ----------

class Recorder:
    """The record of a running ``Basin`` (``[replay] record``): created at the end of its start-up, ``day`` before the
    ledger closes each day, ``close`` at finalize. A restarted run carries on the record it was saved from, its days
    after the restart written again."""

    def __init__(self, basin, restart=None, users: bool = True, heads=None, prepare=None):
        """``users``, ``heads``, ``prepare``: the run's (``Basin``'s arguments, ``prepare`` with the parameters)."""
        cfg = basin.cfg
        from .coupling.modflow6 import Models
        if cfg.surface != "swatplus" or basin.mf6 is None or isinstance(basin.mf6, Models):
            raise ValueError("[replay] records a run of SWAT+ with one MODFLOW 6 groundwater model")
        self.folder = Path(cfg.replay_record)
        self.cast = Precision(cfg.replay_precision, cfg.replay_mantissa_bits)
        self.gwf = _Gwf(basin.mf6)
        self.gwf.spdis = self.gwf.spdis and cfg.replay_specific_discharge  # only dispersion reads it
        self.ident = identity(cfg, basin.parameters, users, heads)
        self.daily, self.start_vars = surface_variables(basin.sw)
        meta_path = self.folder / "record.json"
        if restart is not None:
            if not meta_path.exists():
                raise ValueError(f"[replay] record {self.folder}: a restarted run carries on a record made by the run it "
                                 "was saved from, and there is none")
            self.meta = json.loads(meta_path.read_text())
            refused = compare(self.meta, self.ident)
            ws = self.folder.parent / f".{self.folder.name}-check"
            try:
                build_modflow6(cfg, ws, self.meta["days_total"], dt.date.fromisoformat(self.meta["start"]),
                               prepare)
                if gwf_inputs(ws, cfg.mf6_model) != self.meta["modflow6"]:
                    refused.append("MODFLOW 6's flow model's input files")
            finally:
                shutil.rmtree(ws, ignore_errors=True)
            if refused:
                raise ValueError(f"[replay] record {self.folder} is of another run: {', '.join(refused)} differ")
            if self.meta["days"] < basin.day:
                raise ValueError(f"[replay] record {self.folder} holds {self.meta['days']} days; the run restarts on day "
                                 f"{basin.day}")
            self._keep(basin.day)
        else:
            if self.folder.exists():
                if any(self.folder.iterdir()) and not meta_path.exists():
                    raise ValueError(f"[replay] record {self.folder} is a folder with other files in it; give a new or "
                                     "an old record's folder")
                shutil.rmtree(self.folder)
            (self.folder / "days").mkdir(parents=True)
            self.meta = {"format": FORMAT, "start": basin.start.isoformat(), "days_total": basin.days, "days": 0,
                         "identity": self.ident, "modflow6": gwf_inputs(basin.mf6.ws, cfg.mf6_model),
                         "precision": cfg.replay_precision, "mantissa_bits": self.cast.bits, "chunk_days": None,
                         "bytes": 0, "seconds": 0.0,
                         "nodes": self.gwf.nodes, "nja": self.gwf.nja, "structure": self.gwf.structure}
            self._initial = self._start_state(basin)
            copy_swatplus_files(basin.folder / "swatplus", self.folder / "swatplus")
            self.buffer, self.first = [], 0
        self.static_done = restart is not None

    def _start_state(self, basin) -> dict:
        """What the replay needs before its first day: SWAT+'s start-up reads, and the flow model's memory that GWT
        reads as it prepares the first day."""
        g = self.gwf
        out = {f"start:{v}": basin.sw.get(v).astype(np.float32) for v in self.start_vars}
        out["x0"] = g.get("X").astype(np.float64)
        f = g.get("FLOWJA")
        out["pre:flowja"] = self.cast(f[g.upper])
        out["pre:sat"] = self.cast(g.get("NPF/SAT"))
        if g.spdis:
            out["pre:spdis"] = self.cast(g.get("NPF/SPDIS").ravel())
        return out

    def _keep(self, day: int):
        """A restart's record: the days before ``day`` kept, the chunk it falls in read back for the rest of it."""
        k = self.meta["chunk_days"] or 1
        self.first = (day // k) * k
        path = self.folder / "days" / f"{self.first:06d}.npz"
        self.buffer = read_chunk(path)[: day - self.first] if path.exists() and day > self.first else []
        for p in (self.folder / "days").glob("*.npz"):
            if int(p.stem) > self.first or (int(p.stem) == self.first and not self.buffer):
                p.unlink()
        self.meta["days"] = day

    def _static(self, basin):
        """What doesn't change through the run, written on its first day (the links are built then)."""
        sw, ex, gs = basin.sw, basin.exchange, basin.gw_soil
        info = {"surface": {v: sw.names(v) for v in dict.fromkeys(self.daily + self.start_vars)},
                "daily": self.daily, "start_vars": self.start_vars, "direct_exchange": basin.direct_exchange,
                "packages": self.gwf.packages, "movers": self.gwf.movers, "flow_error": self.gwf.flow_error,
                "storage": [v.split("/")[1].lower() for v in self.gwf.storage], "spdis": self.gwf.spdis,
                "series": [k for k in basin.series if k.split(":")[0] in ("flow", "baseflow", "head")],
                "users": [], "exchange": None, "gw_soil": None}

        def link(pkg, cmap):
            return {"package": pkg, "channels": list(cmap.channels),
                    "entry_channel": np.asarray(cmap.entry_channel).tolist()}
        if ex is not None:
            if ex.kind != "riv" or getattr(ex, "parts", None):
                raise ValueError("[replay] records a stream exchange through one river (RIV) package")
            info["exchange"] = {"kind": ex.kind, **link(ex.package, ex.map),
                                "drains": link(*ex.drains) if ex.drains else None,
                                "lakes": link(*ex.lakes) if ex.lakes else None}
        arrays = dict(self._initial)
        rm = basin.recharge.map if basin.recharge is not None else None
        if rm is not None:
            info["recharge_objects"] = list(rm.objects)
            arrays.update({"map:obj": rm.obj, "map:node": rm.node, "map:area": rm.area})
        if gs is not None:
            info["gw_soil"] = {"package": gs.package, "rule": gs.rule, "cap": gs.cap, "objects": list(gs.map.objects)}
            arrays.update({"soil:obj": gs.map.obj, "soil:node": gs.map.node, "soil:area": gs.map.area,
                           "soil:nodes": gs.nodes})
        for u, wu in basin.users:
            ins = [v for v in wu.bmi.get_input_var_names() if v.startswith("water_user_intake_water_")]
            info["users"].append({"name": u.name, "constituents": list(wu.constituents), "intake": ins})
        (self.folder / "static.json").write_text(json.dumps(info))
        np.savez_compressed(self.folder / "static.npz", **arrays)
        self.static_done = True

    def day(self, basin):
        """The day just run, before the ledger closes."""
        t0 = time.perf_counter()
        if not self.static_done:
            self._static(basin)
        d = self.gwf.state(self.cast)
        for v in self.daily:
            d[f"sw:{v}"] = basin.sw.get(v).astype(np.float32)  # SWAT+'s own 4-byte values: exact
        gs = basin.gw_soil
        if gs is not None and gs.rule == "gwflow":
            d["soil:row_planned"] = self.cast(gs.row_planned)
        L = basin.ledger
        lines = {k: L.today(k) for k in L.rows if k.startswith("user:") and len(L.rows[k]) > L.days}
        lines["recharge:sent"] = L.today("recharge:sent")
        d["json"] = {"ledger": lines, "controls": {u.name: dict(wu.controls) for u, wu in basin.users},
                     "series": {k: basin.series[k][-1] for k in basin.series
                                if k.split(":")[0] in ("flow", "baseflow", "head")},
                     "structure": self.gwf.structure_today()}
        if self.meta["chunk_days"] is None:
            raw = sum(np.asarray(v).nbytes for k, v in d.items() if k != "json")
            self.meta["chunk_days"] = int(max(1, min(64, CHUNK_BYTES // max(raw, 1))))
        self.buffer.append(d)
        if len(self.buffer) == self.meta["chunk_days"]:
            self._flush()
        self.meta["seconds"] += time.perf_counter() - t0

    def _flush(self):
        if not self.buffer:
            return
        path = self.folder / "days" / f"{self.first:06d}.npz"
        old = path.stat().st_size if path.exists() else 0
        self.meta["bytes"] += write_chunk(path, self.buffer) - old
        self.meta["days"] = self.first + len(self.buffer)
        if len(self.buffer) == self.meta["chunk_days"]:
            self.first += len(self.buffer)
            self.buffer = []
        (self.folder / "record.json").write_text(json.dumps(self.meta, indent=1))

    def close(self):
        self._flush()
        (self.folder / "record.json").write_text(json.dumps(self.meta, indent=1))


def compare(meta: dict, ident: dict) -> list[str]:
    """What differs between a record's run and this one, by part."""
    rec = meta["identity"]
    out = [k for k in ("config", "parameters", "users", "heads") if rec[k] != ident[k]]
    files = sorted(set(rec["files"]) | set(ident["files"]))
    out += [f"file {f}" for f in files if rec["files"].get(f) != ident["files"].get(f)]
    return out


# ---------- reading ----------

class Record:
    """A record, read day by day (a chunk at a time)."""

    def __init__(self, folder):
        self.folder = Path(folder)
        try:
            self.meta = json.loads((self.folder / "record.json").read_text())
            self.info = json.loads((self.folder / "static.json").read_text())
        except FileNotFoundError:
            raise ReplayRefused(f"{self.folder} isn't a replay record (record.json, static.json)") from None
        if self.meta.get("format") != FORMAT:
            raise ReplayRefused(f"{self.folder}: record format {self.meta.get('format')}, this Basin Ledger reads "
                                f"{FORMAT}")
        with np.load(self.folder / "static.npz") as z:
            self.static = {k: z[k] for k in z.files}
        self.start = dt.date.fromisoformat(self.meta["start"])
        self.days = int(self.meta["days"])
        self._chunk, self._first = None, None

    def day(self, d: int) -> dict:
        if not 0 <= d < self.days:
            raise ReplayRefused(f"{self.folder}: the record holds days 0 to {self.days - 1}; day {d} asked")
        k = self.meta["chunk_days"]
        first = (d // k) * k
        if first != self._first:
            self._chunk, self._first = read_chunk(self.folder / "days" / f"{first:06d}.npz"), first
        return self._chunk[d - first]

    def bytes_per_cell_day(self) -> float:
        return self.meta["bytes"] / max(self.days, 1) / max(self.meta["nodes"], 1)


class RecordedSurface:
    """SWAT+'s side as the salt balances read it (``coupling.swatplus.Swat``'s ``has``, ``names``, ``get``): the
    start-up values until the first day, then each day's."""

    def __init__(self, rec: Record):
        self.rec = rec
        self.names_ = rec.info["surface"]
        self.values = {v: rec.static[f"start:{v}"].astype(np.float64) for v in rec.info["start_vars"]}

    def has(self, var):
        return var in self.names_

    def names(self, var):
        try:
            return list(self.names_[var])
        except KeyError:
            raise KeyError(f"the replay record holds no {var}") from None

    def get(self, var) -> np.ndarray:
        try:
            return self.values[var].copy()
        except KeyError:
            raise KeyError(f"the replay record holds no {var}") from None

    def set_day(self, day: dict):
        for v in self.rec.info["daily"]:
            self.values[v] = np.asarray(day[f"sw:{v}"], np.float64)

    def set(self, var, values):
        raise ValueError("a replay sets nothing in SWAT+")


class _NoSurface:
    """In place of SWAT+ in a replay: nothing to step, save or finalize."""

    def finalize(self):
        pass

    def save_state(self, folder):
        pass


@dataclasses.dataclass
class _RecordedUser:
    """A user that isn't run again in a replay: its constituents (for the salt's point sources); its lines come from
    the record."""
    constituents: dict

    def finalize(self):
        pass

    def fast_forward(self, day):
        pass


# ---------- replaying ----------

from .driver import Basin  # (driver imports this module lazily, so no cycle)


class ReplayBasin(Basin):
    """A run of the salt (transport, the land and the streams) from a record (module docstring). ``Basin(cfg,
    replay=folder)`` makes one; the arguments are ``Basin``'s, with ``users`` the record's unless given, and
    ``restart`` a state saved by a replay or by a full run. ``heads`` and ``wrap_surface`` are refused (the flow is
    the record's)."""

    def __init__(self, cfg: C.ExampleConfig, run: str = "replay", users: bool | None = None, days: int | None = None,
                 restart: str | Path | None = None, prepare=None, heads=None, wrap_surface=None,
                 parameters: list | None = None, _problems: list | None = None, replay=None):
        try:
            self._replay_start(cfg, run, users, days, restart, prepare, heads, wrap_surface, parameters, replay)
        except BaseException:
            self._finalize_started()
            raise

    def _replay_start(self, cfg, run, users, days, restart, prepare, heads, wrap_surface, parameters, replay):
        from . import driver as D
        from .calibration import apply_parameters, merged, read_parameters
        from .coupling.modflow6 import Modflow6, remaining_simulation
        from .ledger import Ledger
        from .mapping import AreaMap
        from .users import WaterUser, load_component
        self.cfg, self._problems, self.recorder = cfg, None, None
        self.status = cfg.calibration_status
        if heads is not None or wrap_surface is not None:
            raise ReplayRefused("a replay takes the flow from its record: heads and wrap_surface aren't taken")
        rec = self.record = Record(replay)
        info = rec.info
        own = read_parameters(cfg.calibration_parameters) if cfg.calibration_parameters else []
        self.parameters = merged(own, parameters or [])
        users = rec.meta["identity"]["users"] if users is None else users
        ident = identity(cfg, self.parameters, users, None)
        ident["heads"] = rec.meta["identity"]["heads"]  # a replay starts from the recorded run's state, heads and all
        refused = compare(rec.meta, ident)
        if refused:
            raise ReplayRefused(f"replay from {rec.folder}: {', '.join(refused)} differ from the recorded run's, which "
                                "would change the flow; only the salt's settings and parameters may change")

        def prep(model, folder):
            apply_parameters(model, folder, self.parameters, cfg.calibration_zones, cfg.mf6_model,
                             cfg.calibration_pilot_points)
            if prepare:
                prepare(model, folder)
        self.folder = cfg.run_folder / run
        if self.folder.exists():
            shutil.rmtree(self.folder)
        self.folder.mkdir(parents=True)
        saved = json.loads((Path(restart) / "basin.json").read_text()) if restart else None
        self.start, self.day = rec.start, (saved["day"] if saved else 0)
        self.days = days or rec.meta["days_total"]
        if self.days > rec.days:
            raise ReplayRefused(f"replay from {rec.folder}: {self.days} days asked, the record holds {rec.days}")
        self.sw, self.swat, self.weather = RecordedSurface(rec), _NoSurface(), {}
        self.direct_exchange = info["direct_exchange"]

        # MODFLOW 6: the recorded run's simulation, its flow model checked, then left unsolved
        ws = self.folder / "mf6"
        build_modflow6(cfg, ws, rec.meta["days_total"], rec.start, prep)
        if gwf_inputs(ws, cfg.mf6_model) != rec.meta["modflow6"]:
            raise ReplayRefused(f"replay from {rec.folder}: MODFLOW 6's flow model's input files differ from the "
                                "recorded run's, which would change the flow")
        if self.day:
            remaining_simulation(ws, ws, self.day, cfg.mf6_time_units)
        _strip_outputs(ws, cfg.mf6_model)
        self.mf6 = Modflow6(ws, cfg.mf6_model, length_units=cfg.mf6_length_units, time_units=cfg.mf6_time_units)
        self.gwf = _Gwf(self.mf6)
        if (self.gwf.nodes, self.gwf.nja, self.gwf.structure) != (rec.meta["nodes"], rec.meta["nja"],
                                                                  rec.meta["structure"]):
            raise ReplayRefused(f"replay from {rec.folder}: MODFLOW 6's grid differs from the recorded run's")
        if self.gwf.packages != info["packages"]:
            raise ReplayRefused(f"replay from {rec.folder}: the flow model's packages {self.gwf.packages} differ from "
                                f"the record's {info['packages']}")
        if self.gwf.flow_error and not info["flow_error"]:
            raise ReplayRefused("a transport model corrects for the flow's imbalance (FMI FLOW_IMBALANCE_CORRECTION), "
                                "which reads FLOWJA's diagonal, and the record has none: record with it on")
        self.transport = self._transport_models(cfg)
        dispersing = sorted(t.model.lower() for t in self.transport.values()
                            if int(self.mf6.api.get_value(f"{t.model}/INDSP")[0]) != 0)  # a DSP package's unit
        if dispersing and not info["spdis"]:
            raise ReplayRefused(f"replay: {', '.join(dispersing)} disperse, which reads the specific discharge, and the "
                                "record was made without it ([replay] specific_discharge = false)")
        for sp, t in self.transport.items():
            if saved:
                t.set_concentration(np.load(Path(restart) / f"concentration_{sp}.npy"))
        self.x0 = rec.static["x0"]
        pre = (rec.day(self.day - 1) if self.day else
               {k[4:]: v for k, v in rec.static.items() if k.startswith("pre:")} | {"dry": np.flatnonzero(
                   self.x0 == DHDRY)})
        self.gwf.put(pre, self.x0, before=True)

        # the links, as recorded
        self.ledger = L = Ledger()
        self.map = self.ground_map = None
        self.recharge = None
        if "recharge_objects" in info:
            self.map = AreaMap(info["recharge_objects"], rec.static["map:obj"], rec.static["map:node"],
                               rec.static["map:area"])
            self.recharge = SimpleNamespace(map=self.map)
        self.gw_soil = None
        if info["gw_soil"] is not None:
            from .coupling.gwsoil import cell_shares
            g = info["gw_soil"]
            gmap = AreaMap(g["objects"], rec.static["soil:obj"], rec.static["soil:node"], rec.static["soil:area"])
            self.gw_soil = SimpleNamespace(package=g["package"], rule=g["rule"], cap=g["cap"], map=gmap,
                                           nodes=rec.static["soil:nodes"], shares=cell_shares(gmap, self.mf6.nodes),
                                           planned=None, row_planned=None, due_volume=lambda: None)

        def link(x):
            return None if x is None else (x["package"], SimpleNamespace(
                channels=x["channels"], entry_channel=np.asarray(x["entry_channel"], int)))
        ex = info["exchange"]
        self.exchange = None if ex is None else SimpleNamespace(
            kind=ex["kind"], package=ex["package"], map=link(ex)[1], drains=link(ex["drains"]), lakes=link(ex["lakes"]))
        self.due, self.lakes, self.irrigation, self.rivers, self.rules = {}, None, None, None, []
        # the stage-2 species' state (coupling/surface_mass.py): none in a replay, which carries salt only
        self.surface, self._track_transport, self.swat_aquifer_ha = {}, False, {}
        self.rule_memory = saved["rule_memory"] if saved else {}
        self._wells_left_out = set()
        if cfg.ground or any(w.cells for w in cfg.wells):
            raise ReplayRefused("replay: returns to the ground and wells over several cells aren't built in a replay "
                                "(it re-runs only the users who take in the aquifer's water)")
        wells = {w.name: (w.package, self._well_node(w)) for w in cfg.wells}
        # the routes as the transport's links read them (driver.Basin's router): wells only
        self.router = SimpleNamespace(wells=wells, grounds={}, groups={}, used=set(), req={})
        self._ground_kg, self._ground_cfg, self._well_c = {}, {}, {}

        # users: those taking in the aquifer's water run again; the others' lines are the record's
        recorded = {u["name"]: u for u in info["users"]}
        self.users, self.rerun, self.intakes = [], set(), []
        for u in (cfg.users if users else []):
            r = recorded.get(u.name)
            if r is None:
                raise ReplayRefused(f"replay: user {u.name} isn't in the record")
            takes = [sp for sp in self.transport if f"water_user_intake_water_{sp}__mass_concentration" in r["intake"]]
            k, name = D.kind(u.withdrawal)
            if takes and k == "well" and name in wells:
                if u.forcing:
                    raise ReplayRefused(f"replay: user {u.name} takes in the aquifer's water and has forcing; running it "
                                        "again from a record isn't built")
                comp = load_component(u.component)
                wu = WaterUser(u.name, comp, str(u.config) if u.config else None, variables=u.variables or None,
                               constituents=u.constituents, report=u.report or None)
                wu.inputs = self._forcing(u, comp)
                self.users.append((u, wu))
                self.rerun.add(u.name)
                L.pair(f"user:{u.name}:withdrawal", [f"user:{u.name}:{x}" for x in ("consumption", "return",
                                                                                    "storage_less_gain")])
                for sp in takes:
                    self.intakes.append((u.name, wu, f"water_user_intake_water_{sp}__mass_concentration", sp, name))
            else:
                self.users.append((u, _RecordedUser({c: c for c in r["constituents"]})))
        self._set_intakes()

        self.salt = {}
        for sp in cfg.salt:
            self._salt(cfg, sp, rec.folder / "swatplus", restart)

        # outputs: the recorded series, and the salt's
        channels = rec.info["surface"].get("channel_water__volume_outflow_rate", [])
        self.out_channels = cfg.outputs.get("channels", channels)
        self._ch_index, self._nexus, self.out_heads = {}, {}, []
        self.series = {k: [] for k in info["series"]}
        for sp in self.salt:
            self.series.update({f"{sp}_mg_l:{c}": [] for c in self.out_channels})
        if self.salt and cfg.salt_conductance is not None:
            self.series.update({f"conductance_us_cm:{c}": [] for c in self.out_channels})
        if saved:
            L.load(saved["ledger"])
            self.series = {k: v for k, v in saved["series"].items() if k in self.series} | {
                k: [] for k in self.series if k not in saved["series"]}
            for _, wu in self.users:
                wu.fast_forward(self.day)
        self.times = {"record": 0.0, "salt": 0.0, "modflow6": 0.0}

    # ---------- one day ----------
    def step(self):
        d, L, rec = self.day, self.ledger, self.record
        t0 = time.perf_counter()
        day = rec.day(d)
        j = day["json"]
        self.sw.set_day(day)
        t1 = time.perf_counter()
        for u, wu in self.users:
            k = f"user:{u.name}:"
            if u.name not in self.rerun:
                for a, v in j["ledger"].items():
                    if a.startswith(k):
                        L.enter(a, v)
                continue
            wu.set_controls(**j["controls"].get(u.name, {}))
            use = wu.advance_day(d)
            for a, v in (("withdrawal", use.withdrawal), ("gain", use.gain), ("storage_less_gain",
                         use.storage_change - use.gain), ("consumption", use.consumption), ("return", use.returned)):
                was = j["ledger"].get(k + a, 0.0)
                if abs(v - was) > 1e-12 * max(abs(was), 1.0):
                    raise ReplayRefused(f"replay day {d}: user {u.name}'s {a} is {v!r}, the recorded run's {was!r}: its "
                                        "water depends on what it takes in from the aquifer, so changing the salt "
                                        "would change the flow")
                L.enter(k + a, v)
            for c, kg in use.constituents.items():
                L.enter(k + f"{c}:kg", kg)
            for sp in self.transport:  # the intake from the aquifer, as Basin.step books it
                if (u.name, sp) in self.intake:
                    L.enter(k + f"intake:{sp}", self.intake[(u.name, sp)])
                    L.enter(k + f"intake:{sp}:kg", use.withdrawal * self.intake[(u.name, sp)])
            if use.energy is not None:
                L.enter(k + "MWh", use.energy)
            for name, v in use.reported.items():
                L.enter(k + name, v)
        reads = {}
        for s in self.salt.values():
            s.after_surface(L, self.start + dt.timedelta(days=d), reads)
        t2 = time.perf_counter()

        self.mf6.prepare_day()
        today = self.gwf.structure_today()
        if today != j["structure"]:
            bad = sorted(p for p in today if today[p] != j["structure"].get(p))
            raise ReplayRefused(f"replay day {d}: the flow model's entries in {bad} aren't where the recorded run had them")
        # the recharge's water as the recorded run sent it, before the transport's inflows read it (driver's order)
        L.enter("recharge:sent", j["ledger"].get("recharge:sent", 0.0))
        if self.transport:
            self._transport_day(L)
        self.gwf.put(day, self.x0)
        gs = self.gw_soil
        if gs is not None and gs.rule == "gwflow":
            gs.row_planned = np.asarray(day["soil:row_planned"], np.float64)
            gs.planned = np.bincount(gs.map.node, weights=gs.row_planned, minlength=self.mf6.nodes)[gs.nodes]
        api = self.mf6.api
        for sol in range(2, api.get_subcomponent_count() + 1):  # the transport's solutions; the flow's (1) isn't solved
            self.mf6._solve(sol)
        api.finalize_time_step()
        self.mf6.day += 1
        if self.transport:
            self._transport_accounts(L)
        t3 = time.perf_counter()

        for k in self.series:
            if k in j["series"]:
                self.series[k].append(j["series"][k])
        if self.salt:
            from .salt import conductance
            sc = self.cfg.salt_conductance
            for c in self.out_channels:
                mg_l = {sp: s.concentration(c) * 1000.0 for sp, s in self.salt.items()}
                for sp, v in mg_l.items():
                    self.series[f"{sp}_mg_l:{c}"].append(v)
                if sc is not None:
                    self.series[f"conductance_us_cm:{c}"].append(conductance(sc["base"], sc["slopes"], mg_l))
        L.close_day()
        self.day += 1
        self.times["record"] += t1 - t0
        self.times["salt"] += t2 - t1
        self.times["modflow6"] += t3 - t2

    # ---------- restarts and outputs ----------
    def save(self, folder: str | Path):
        """The replay's state after the last day run: the concentrations, the salt, the ledger and the series, for a
        replay's ``restart`` (not a full run's: there is no SWAT+ or flow state in it)."""
        folder = Path(folder)
        if folder.exists():
            shutil.rmtree(folder)
        folder.mkdir(parents=True)
        for sp, t in self.transport.items():
            np.save(folder / f"concentration_{sp}.npy", t.concentration())
        for sp, s in self.salt.items():
            np.savez(folder / f"salt_{sp}.npz", **s.state())
        (folder / "basin.json").write_text(json.dumps({
            "day": self.day, "ledger": self.ledger.state(), "series": self.series, "rule_memory": self.rule_memory,
            "replay": str(self.record.folder)}))

    def write_outputs(self, folder: str | Path | None = None) -> Path:
        out = super().write_outputs(folder)
        s = json.loads((out / "summary.json").read_text())
        s["replay"] = {"record": str(self.record.folder), "note": "salt replayed from a recorded flow run; the water "
                       "accounts are the recorded run's (its own outputs)"}
        (out / "summary.json").write_text(json.dumps(s, indent=2))
        return out

    def finalize(self):
        for _, wu in self.users:
            wu.finalize()
        if self.mf6 is not None:
            self.mf6.finalize()


def replay_basin(cfg, replay, **kw) -> ReplayBasin:
    """``Basin(cfg, replay=replay, **kw)``."""
    with warnings.catch_warnings():
        return ReplayBasin(cfg, replay=replay, **kw)
