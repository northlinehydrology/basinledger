"""What Basin Ledger takes over in a MODFLOW 6 model, and what it leaves as the model's (docs/design.md, "Extending an
existing MODFLOW 6 model").

A model someone else built comes with its own boundaries: recharge from its own estimate, rivers, drains, general
and fixed heads, its wells, evapotranspiration. Basin Ledger sets only what it couples:

- recharge, in the package ``[recharge] package`` names (every entry, from SWAT+);
- the stream exchange, on the entries of ``[exchange] package`` (and its drains) that the map gives a SWAT+ channel;
- groundwater up into the soil, in ``[groundwater_to_soil] package``;
- users' wells, the entries ``[[wells]]`` names.

Everything else stays the model's, and the ledger books its flow each day as water from or to outside the coupled
models (``groundwater:boundary:<package>``), as it always has. The start-up check lists each package with what
it was taken over for and how many of its entries, so a reviewer sees what the coupling changed. It also refuses
what would count water twice:

- another recharge package left as the model's (RCH): the model's own recharge would be added to SWAT+'s;
- a UZF package: it takes infiltration of its own, which is recharge again (coupling UZF is not built);
- an EVT package with groundwater up into the soil: SWAT+'s plants then take groundwater from the soil, and EVT
  would take it a second time.

And it notes, without refusing, an EVT package without that link (groundwater ET stays the model's), and an SFR or
LAK package left as the model's while SWAT+'s channels are coupled: their water is the model's, so if they stand for
the same streams or lakes as SWAT+'s channels, the river is carried twice.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


class TakeoverNote(UserWarning):
    """Something about the model's own packages a reviewer should know: a start-up warning, listed by a check."""


@dataclass
class PackageRole:
    package: str         # the package's name in the model, lower case
    type: str            # its type, lower case: rch, riv, drn, ghb, chd, wel, evt, sfr, lak, maw, uzf, ...
    entries: int         # its entries this period (reaches for SFR, connections for LAK and MAW)
    taken: int           # how many Basin Ledger sets
    what: str            # what Basin Ledger takes it over for, or "" (left as the model's)

    @property
    def role(self) -> str:
        if self.taken == 0:
            return "the model's"
        if self.taken == self.entries:
            return "Basin Ledger's"
        return "partly Basin Ledger's"

    def describe(self) -> str:
        left = f"booked as groundwater:boundary:{self.package}"
        if self.taken == 0:
            return f"the model's: water from or to outside the coupled models, {left}"
        if self.taken == self.entries:
            return f"Basin Ledger's: {self.what}"
        return (f"{self.taken} of {self.entries} entries Basin Ledger's ({self.what}); the other "
                f"{self.entries - self.taken} the model's, {left}")


def flow_packages(mf6, model: str | None = None) -> list[str]:
    """A groundwater model's packages with flows (SIMVALS), upper case, sorted: ``model``, or the one Basin Ledger
    couples."""
    prefix = (model or mf6.model).upper() + "/"
    names = set(mf6.api.get_input_var_names()) | set(mf6.api.get_output_var_names())
    return sorted({n.split("/")[1] for n in names if n.startswith(prefix) and n.endswith("/SIMVALS")})


def inventory(mf6, model: str | None = None) -> list[PackageRole]:
    """Each package with flows in a model as given, before anything is coupled: what a first look at someone's model
    lists (call after ``prepare_day``, when the packages have their entries). ``model``: any groundwater model of the
    simulation, the one ``mf6`` couples if not given."""
    from .exchange import package_type
    m = (model or mf6.model).upper()
    out = []
    for p in flow_packages(mf6, m):
        n = int(mf6.get("NBOUND", m, p)[0])
        t = str(mf6.get("PACKAGE_TYPE", m, p)[0]).strip().lower()
        out.append(PackageRole(p.lower(), t, n, 0, ""))
    return out


def roles(mf6, coupled: dict[str, tuple[np.ndarray, list[str]]]) -> list[PackageRole]:
    """Each package with flows (SIMVALS), in the model's order of names: ``coupled`` maps a package's upper-case name to
    (a mask of its entries Basin Ledger sets, the links that set them)."""
    from .exchange import package_type
    out = []
    for p, (mask, links) in coupled.items():
        out.append(PackageRole(p.lower(), package_type(mf6, p), int(len(mask)), int(mask.sum()),
                               ", ".join(dict.fromkeys(links))))
    return out


OVER_CELL = 1e-3  # SWAT+'s share of a cell plus a package's own, beyond the whole cell, that is counted twice


def recharge_cover(share: np.ndarray, nodes: np.ndarray, mult: np.ndarray,
                   cell_area: np.ndarray) -> tuple[int, int, int, float]:
    """Where a recharge package left as the model's lies against SWAT+'s recharge. ``share``: each cell's (by node)
    share SWAT+'s recharge map covers; ``nodes`` and ``mult``: the package's entries' nodes and multipliers (its
    auxiliary multiplier, the share of the cell it stands for, or 1). On a cell where the two add up to more than the
    whole cell (by more than ``OVER_CELL``) the rest is recharged twice. Returns (entries, cells over, cells shared
    within the cell, the area recharged twice in m2)."""
    ok = nodes >= 0
    n, m = nodes[ok], np.asarray(mult, float)[ok]
    c = share[n]
    over = c + m - 1.0
    twice = float(np.sum(np.clip(over, 0.0, None) * cell_area[n]))
    return len(nodes), int((over > OVER_CELL).sum()), int(((c > 0) & (over <= OVER_CELL)).sum()), twice


def problems(rows: list[PackageRole], groundwater_to_soil: bool,
             cover: dict[str, tuple[int, int, int]] | None = None) -> tuple[list[str], list[str]]:
    """What the model's own packages would count twice (problems), and what a reviewer should know (notes).
    ``cover``: for each recharge package left as the model's, ``recharge_cover``'s (entries, cells over, cells
    shared, area recharged twice). On ground SWAT+'s recharge doesn't cover (a catchment outside SWAT+'s project) the
    package's recharge is the model's, booked as water from outside, and on a cell SWAT+ covers in part it may stand
    for the rest (its auxiliary multiplier the rest's share); where SWAT+'s share and its own add up to more than the
    cell, that ground is recharged twice, which is refused. Without ``cover`` (no recharge link) a recharge package
    left as the model's is refused, as SWAT+'s recharge can't be told apart."""
    bad, notes = [], []
    cover = cover or {}
    for r in rows:
        if r.taken:
            continue
        if r.type == "rch" and r.package in cover and cover[r.package][1] == 0:
            n, _, shared, _ = cover[r.package]
            part = (f"; {shared} of them shared with SWAT+'s recharge, the package standing for the rest of the "
                    "cell") if shared else ""
            notes.append(f"{r.package}: recharge left as the model's on {n} cells off SWAT+'s recharge, booked as "
                         f"water from outside{part}")
        elif r.type == "rch" and r.package in cover:
            _, over, _, twice = cover[r.package]
            bad.append(f"{r.package}: a recharge package left as the model's recharges {over} cells SWAT+'s recharge "
                       f"also covers, beyond the whole cell ({twice:.0f} m2 recharged twice). Its multiplier "
                       "(AUXMULTNAME) on a cell SWAT+ covers in part is the share of the cell left to it; name it in "
                       "[recharge] package if it is the one SWAT+'s recharge replaces")
        elif r.type == "rch":
            bad.append(f"{r.package}: a recharge package left as the model's; its recharge would be added to SWAT+'s. "
                       "Name it in [recharge] package if it is the one SWAT+'s recharge replaces, or leave it out of "
                       "the run's copy of the model")
        elif r.type == "uzf":
            bad.append(f"{r.package}: UZF takes infiltration of its own, which is recharge again; coupling UZF to "
                       "SWAT+'s percolation isn't built")
        elif r.type == "evt" and groundwater_to_soil:
            bad.append(f"{r.package}: EVT takes groundwater ET, and with [groundwater_to_soil] SWAT+'s plants take it "
                       "from the soil too; leave EVT out of the run's copy of the model")
        elif r.type == "evt":
            notes.append(f"{r.package}: groundwater ET stays the model's (EVT), apart from SWAT+'s ET from the soil")
        elif r.type in ("sfr", "lak"):
            notes.append(f"{r.package}: the model's {r.type.upper()} carries water SWAT+ doesn't see; if it stands for "
                         "the same streams or lakes as SWAT+'s channels, that water is carried twice")
    return bad, notes


def table(rows: list[PackageRole]) -> str:
    """The start-up check's list, one package per line."""
    w = max([len(r.package) for r in rows] + [7])
    lines = [f"{'package':<{w}}  type  {'entries':>7}  taken over, or left as the model's"]
    for r in rows:
        lines.append(f"{r.package:<{w}}  {r.type:<4}  {r.entries:>7}  {r.describe()}")
    return "\n".join(lines)



def _model_nam(ws, model: str):
    """The model's name file, the package names it has (lower case) and its grid type ("dis", "disv", "disu")."""
    from .mover import _blocks
    nam = None
    for name, _, lines in _blocks(ws / "mfsim.nam", ws):
        if name == "MODELS":
            for words in lines:
                if words[2].lower() == model.lower():
                    nam = ws / words[1].strip("'\"")
    if nam is None:
        raise ValueError(f"[modflow6] model {model!r} isn't in mfsim.nam")
    have, grid = set(), None
    for name, _, lines in _blocks(nam, ws):
        if name == "PACKAGES":
            for words in lines:
                if len(words) > 2:
                    have.add(words[2].lower())
                if words[0].upper() in ("DIS6", "DISV6", "DISU6"):
                    grid = words[0].upper()[:-1].lower()
    return nam, have, grid


def _add_package(ws, model: str, nam, ftype: str, package: str, cells: list[tuple[int, ...]], what: str):
    """A list package of Basin Ledger's own, one entry per cell with a value of 0 (set every day), and its line in the
    model's name file."""
    f = f"{model.lower()}.{package}.{ftype.lower()}"
    rows = "\n".join("  " + " ".join(str(i + 1) for i in c) + "  0.0" for c in cells)
    (ws / f).write_text(f"# Basin Ledger: {what} (coupling/existing.py); each value set every day\n"
                        "BEGIN options\nEND options\n\nBEGIN dimensions\n"
                        f"  MAXBOUND  {len(cells)}\nEND dimensions\n\nBEGIN period  1\n{rows}\nEND period  1\n")
    text = nam.read_text().splitlines()
    end = next(i for i, x in enumerate(text) if [w.upper() for w in x.split()[:2]] == ["END", "PACKAGES"])
    nam.write_text("\n".join(text[:end] + [f"  {ftype.upper()}6  {f}  {package}"] + text[end:]) + "\n")


def add_recharge(ws, model: str, package: str, cells: list[tuple[int, ...]]) -> list[str]:
    """A recharge package of Basin Ledger's own for a model without one (``[recharge] add = true``), added to a run's
    copy: an RCH package named ``package``, one entry per cell of the recharge map (0-based, in the grid's form), its
    rate set every day from SWAT+. A package name the model has and cells in the wrong form are refused; a model
    that has a recharge package of its own is refused at the takeover (it would count recharge twice). Returns the
    start-up's notes."""
    from pathlib import Path
    ws = Path(ws)
    nam, have, grid = _model_nam(ws, model)
    size = {"dis": 3, "disv": 2, "disu": 1}[grid]
    if package.lower() in have:
        raise ValueError(f"[recharge] add = true: the model has a package {package!r}; name the one Basin Ledger adds "
                         "otherwise, or leave add out to set the model's own")
    cells = list(dict.fromkeys(tuple(c) for c in cells))
    if not cells or any(len(c) != size for c in cells):
        raise ValueError(f"[recharge] add = true needs the map's cells in {grid.upper()} form")
    _add_package(ws, model, nam, "rch", package.lower(), cells, "recharge the model doesn't have")
    return [(f"{package.lower()}: a recharge package Basin Ledger adds to the run's copy, the model having none: "
             f"{len(cells)} cells of the recharge map")]


def add_wells(ws, model: str, wells: list[tuple[str, str, tuple[int, ...]]]) -> list[str]:
    """Wells a model doesn't have, added to a run's copy (``[[wells]] add = true``): for each package named, a WEL
    package of Basin Ledger's own, one entry per well at its cell (0-based, in the grid's form), its rate 0 in the
    file and set every day by the users' routes, and its line in the model's name file. ``wells``: (name, package,
    cell). A package name the model already has, a cell in the wrong form for the grid and two wells at one cell in
    a package are refused. Returns the start-up's notes."""
    from pathlib import Path
    ws = Path(ws)
    nam, have, grid = _model_nam(ws, model)
    size = {"dis": 3, "disv": 2, "disu": 1}[grid]
    by_package: dict[str, list[tuple[str, tuple[int, ...]]]] = {}
    for well, package, cell in wells:
        if package.lower() in have:
            raise ValueError(f"well {well!r}: the model has a package {package!r}; a well Basin Ledger adds goes in a "
                             "package of its own (another name)")
        if len(cell) != size:
            raise ValueError(f"well {well!r}: cell {list(cell)} isn't a {grid.upper()} cell")
        by_package.setdefault(package.lower(), []).append((well, tuple(cell)))
    notes = []
    for package, added in by_package.items():
        cells = [c for _, c in added]
        if len(set(cells)) != len(cells):
            raise ValueError(f"{package}: two wells Basin Ledger adds share a cell; give them one well")
        _add_package(ws, model, nam, "wel", package, cells, "wells the model doesn't have")
        notes.append(f"{package}: a well package Basin Ledger adds to the run's copy for wells the model doesn't have: "
                     + ", ".join(f"{w} at {list(c)}" for w, c in added))
    return notes


def _package_file(ws, nam, ftypes: tuple[str, ...]):
    """The first package of these types in the model's name file: (its file, its name), or (None, None). A package
    the name file doesn't name is called by its file's name."""
    from .mover import _blocks
    for name, _, lines in _blocks(nam, ws):
        if name == "PACKAGES":
            for words in lines:
                if words[0].upper() in ftypes:
                    f = ws / words[1].strip("'\"")
                    return f, (words[2].lower() if len(words) > 2 else f.name)
    return None, None


def _length_units(ws, nam, stated: str | None) -> str:
    """The model's length units from its DIS file's options (read only as far as they go: a DIS file can be
    hundreds of megabytes), else the config's (``[modflow6] length_units``)."""
    from .modflow6 import METERS
    dis, _ = _package_file(ws, nam, ("DIS6", "DISV6", "DISU6"))
    own = None
    if dis is not None:
        with open(dis, errors="replace") as f:
            inside = False
            for line in f:
                w = line.split()
                if not w or w[0].startswith(("#", "!", "//")):
                    continue
                if w[0].upper() == "BEGIN":
                    inside = len(w) > 1 and w[1].upper() == "OPTIONS"
                elif w[0].upper() == "END":
                    if inside:
                        break
                elif inside and w[0].upper() == "LENGTH_UNITS" and len(w) > 1:
                    own = w[1].lower()
    if own in METERS:
        return own
    if stated in METERS:
        return stated
    raise ValueError("[modflow6] storage: the model's length units are undefined (DIS LENGTH_UNITS); state them, "
                     "[modflow6] length_units")


def steady_periods(ws, model: str) -> str | None:
    """Whether a model solves any day without storage: "no STO" (no storage package: every period steady),
    "steady" (a STO package with a STEADY-STATE period, or no period block, which MODFLOW 6 takes as steady), or
    None (transient throughout)."""
    from pathlib import Path

    from .mover import _blocks
    ws = Path(ws)
    nam, _, _ = _model_nam(ws, model)
    sto, _ = _package_file(ws, nam, ("STO6",))
    if sto is None:
        return "no STO"
    periods = [(words, lines) for name, words, lines in _blocks(sto, ws) if name == "PERIOD"]
    if not periods or any(w and w[0].upper() == "STEADY-STATE" for _, lines in periods for w in lines):
        return "steady"
    return None


def add_storage(ws, model: str, sy: float | None, ss: float | None, iconvert: int = 1,
                length_units: str | None = None) -> list[str]:
    """Storage for a model handed over steady, in a run's copy (``[modflow6] storage``): daily steps need specific
    yield and specific storage, which a steady model never needed (the study's [A]). A model without a STO package
    gets one of Basin Ledger's own, every cell's ``iconvert``, ``sy`` and ``ss`` (given per meter, written in the
    model's length units), transient from the first period; a model with one has its steady periods made transient
    and keeps its own values (``sy`` and ``ss`` given as well are refused: set them in its file). Returns the
    start-up's notes."""
    import re
    from pathlib import Path

    from .modflow6 import METERS
    from .mover import _blocks
    ws = Path(ws)
    nam, have, _ = _model_nam(ws, model)
    sto, pname = _package_file(ws, nam, ("STO6",))
    if sto is None:
        if sy is None or ss is None:
            raise ValueError(f"[modflow6] storage: {model} has no storage package; give sy and ss, the study's values")
        if not (0 <= sy <= 1 and ss >= 0 and iconvert in (0, 1)):
            raise ValueError("[modflow6] storage: sy is a fraction (0 to 1), ss per meter (0 or more), iconvert 0 or 1")
        meters = METERS[_length_units(ws, nam, length_units)]
        name = "sto_basinledger" if "sto_basinledger" not in have else "sto_basinledger_1"
        f = f"{model.lower()}.{name}.sto"
        (ws / f).write_text(f"# Basin Ledger: storage for a model handed over steady (coupling/existing.py)\n"
                            "BEGIN options\nEND options\n\nBEGIN griddata\n"
                            f"  ICONVERT\n    CONSTANT {iconvert}\n  SS\n    CONSTANT {ss * meters:.10g}\n"
                            f"  SY\n    CONSTANT {sy:.10g}\nEND griddata\n\nBEGIN period  1\n  TRANSIENT\nEND period  1\n")
        text = nam.read_text().splitlines()
        end = next(i for i, x in enumerate(text) if [w.upper() for w in x.split()[:2]] == ["END", "PACKAGES"])
        nam.write_text("\n".join(text[:end] + [f"  STO6  {f}  {name}"] + text[end:]) + "\n")
        return [(f"{name}: a storage package Basin Ledger adds to the run's copy, the model having none: ICONVERT "
                 f"{iconvert}, SY {sy:g}, SS {ss:g} per meter, transient from the first day")]
    if sy is not None or ss is not None:
        raise ValueError(f"[modflow6] storage: {model} has its own storage package ({pname}); its values stay "
                         "the model's (set them in its file), and storage = {} only makes its steady periods transient")
    text = sto.read_text()
    periods = [p for name, p, _ in _blocks(sto, ws) if name == "PERIOD"]
    n = len(re.findall(r"(?im)^\s*STEADY-STATE\b", text))
    text = re.sub(r"(?im)^(\s*)STEADY-STATE\b", r"\1TRANSIENT", text)
    if not periods:  # MODFLOW 6 takes a STO package without period blocks as steady
        text = text.rstrip("\n") + "\n\nBEGIN period  1\n  TRANSIENT\nEND period  1\n"
    if n == 0 and periods:
        return []
    sto.write_text(f"# Basin Ledger: steady periods made transient (coupling/existing.py)\n{text}")
    return [f"{pname}: the model's storage package made transient in the run's copy ("
            + (f"{n} steady period{'s' if n != 1 else ''}" if n else "no period block, which is steady")
            + "), its own SY and SS kept"]
