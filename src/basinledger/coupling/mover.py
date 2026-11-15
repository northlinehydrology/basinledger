"""MODFLOW 6's mover (MVR) read from its files, so the SFR line-up rule can treat its moves as connections
(coupling/exchange.py, ``sfr_problems``; docs/design.md, "Extending an existing MODFLOW 6 model").

A mover takes water from one package's entry (a provider: an SFR reach's outflow, a drain's, a lake's outlet...) and
gives it to another's (a receiver), by a factor or an amount each stress period. A model's own mover (``MVR6`` in its
name file) moves water within the model; a GWF-GWF exchange's (``MVR6`` in the exchange's options) moves it from one
model to the other, as Pleasant Lake's carries the inset's streams into the parent's.

libmf6 keeps each mover's types and values in memory (``<model or exchange>/MVR/IMVRTYPE`` and ``VALUE``, one per
move in the stress period's block, in its order) but not which entries they join, so those are read here from the
mover's file. A move can be cut for a day, after the day's period data are read, by setting its value so it moves
nothing: 0 for a factor or an upper limit (FACTOR, UPTO), more than any flow for an excess or a threshold (EXCESS,
THRESHOLD). The provider keeps the water then: an SFR reach's outflow goes on down its own connections, or leaves SFR
as its external outflow where it has none, as a cut upstream fraction does.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

KINDS = ("FACTOR", "EXCESS", "THRESHOLD", "UPTO")
NOTHING = {"FACTOR": 0.0, "UPTO": 0.0, "EXCESS": 1e30, "THRESHOLD": 1e30}  # a value that moves nothing


@dataclass(frozen=True)
class Entry:
    """A package's entry: model and package names lower case, the entry 0-based."""
    model: str
    package: str
    id: int

    def __str__(self):
        return f"{self.model}'s {self.package} entry {self.id}"


@dataclass(frozen=True)
class Move:
    """One row of a mover's period block."""
    mover: str     # the mover's memory path: "GWF/MVR", or "GWF-GWF_1/MVR" for an exchange's
    period: int    # the period block it is in (1-based); it holds until the next block
    index: int     # its place in that block, and in libmf6's IMVRTYPE and VALUE
    provider: Entry
    receiver: Entry
    kind: str      # FACTOR, EXCESS, THRESHOLD or UPTO
    value: float


@dataclass
class Mover:
    path: str                       # memory path
    file: Path
    blocks: dict[int, list[Move]]   # first period -> its moves

    def active(self, period: int) -> list[Move]:
        """The moves in effect in a stress period (1-based): the last block at or before it."""
        k = [p for p in self.blocks if p <= period]
        return self.blocks[max(k)] if k else []

    @property
    def moves(self) -> list[Move]:
        return [m for b in self.blocks.values() for m in b]


def _lines(path: Path, ws: Path) -> list[list[str]]:
    """The file's lines as words, comments and blank lines left out, OPEN/CLOSE files read in place."""
    out = []
    for line in path.read_text().splitlines():
        words = line.split("#", 1)[0].split("!", 1)[0].split()
        if not words:
            continue
        if words[0].upper() == "OPEN/CLOSE":
            out.extend(_lines(ws / words[1].strip("'\""), ws))
        else:
            out.append(words)
    return out


def _blocks(path: Path, ws: Path) -> list[tuple[str, list[str], list[list[str]]]]:
    """(block name, the BEGIN line's other words, the block's lines) in file order."""
    out, cur = [], None
    for words in _lines(path, ws):
        key = words[0].upper()
        if key == "BEGIN":
            cur = (words[1].upper(), words[2:], [])
            out.append(cur)
        elif key == "END":
            cur = None
        elif cur is not None:
            cur[2].append(words)
    return out


def movers(ws: str | Path) -> list[Mover]:
    """Every mover in the simulation in ``ws``: each model's own, then each exchange's."""
    ws = Path(ws)
    sim = _blocks(ws / "mfsim.nam", ws)
    found: list[tuple[str, Path, str | None]] = []  # (memory path, file, the model a model's own mover is in)
    for name, _, lines in sim:
        if name == "MODELS":
            for words in lines:
                for b, _, plines in _blocks(ws / words[1].strip("'\""), ws):
                    if b == "PACKAGES":
                        for p in plines:
                            if p[0].upper() == "MVR6":
                                found.append((f"{words[2].upper()}/MVR", ws / p[1].strip("'\""), words[2].lower()))
    n = 0
    for name, _, lines in sim:
        if name == "EXCHANGES":
            for words in lines:
                n += 1  # libmf6 names an exchange by its type and its place in the block: GWF-GWF_1
                kind = words[0].upper().replace("6", "")
                for b, _, olines in _blocks(ws / words[1].strip("'\""), ws):
                    if b == "OPTIONS":
                        for o in olines:
                            if o[0].upper() == "MVR6":
                                f = o[[w.upper() for w in o].index("FILEIN") + 1] if len(o) > 2 else o[1]
                                found.append((f"{kind}_{n}/MVR", ws / f.strip("'\""), None))
    return [_read(path, f, model, ws) for path, f, model in found]


def _read(path: str, f: Path, model: str | None, ws: Path) -> Mover:
    blocks: dict[int, list[Move]] = {}
    modelnames = False
    for name, head, lines in _blocks(f, ws):
        if name == "OPTIONS":
            modelnames = any(w[0].upper() == "MODELNAMES" for w in lines)
        elif name == "PERIOD":
            k = int(head[0])
            rows = []
            for i, w in enumerate(lines):
                if modelnames:
                    p = Entry(w[0].lower(), w[1].lower(), int(w[2]) - 1)
                    r = Entry(w[3].lower(), w[4].lower(), int(w[5]) - 1)
                    kind, value = w[6].upper(), float(w[7])
                else:
                    p = Entry(model, w[0].lower(), int(w[1]) - 1)
                    r = Entry(model, w[2].lower(), int(w[3]) - 1)
                    kind, value = w[4].upper(), float(w[5])
                if kind not in KINDS:
                    raise ValueError(f"{f.name}: period {k}, move {i + 1}: type {kind}, not one of {KINDS}")
                rows.append(Move(path, k, i, p, r, kind, value))
            blocks[k] = rows
    return Mover(path, f, blocks)


def cut(mf6, movers: list[Mover], moves: set[Move] | list[Move]):
    """Set each of ``moves``, where in effect this stress period, to move nothing today: call after the day's period
    data are read (``Modflow6.prepare_day``) and before the solve. A value of 0 or 1e30 means the same in any units."""
    moves = set(moves)
    period = mf6.period
    for mv in movers:
        active = mv.active(period)
        todo = [m for m in active if m in moves]
        if not todo:
            continue
        value = mf6.api.get_value(f"{mv.path}/VALUE").copy()
        n = int(mf6.api.get_value(f"{mv.path}/NMVR")[0])
        if n != len(active):
            raise ValueError(f"{mv.file.name}: libmf6 has {n} moves in period {period}, the file {len(active)}")
        for m in todo:
            value[m.index] = NOTHING[m.kind]
        mf6.api.set_value(f"{mv.path}/VALUE", np.asarray(value, dtype=np.float64))
