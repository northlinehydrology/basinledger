"""Maps on a MODFLOW 6 grid Basin Ledger didn't write (docs/design.md, "Extending an existing MODFLOW 6 model").

Someone else's model comes with its own grid: structured (DIS) or vertex (DISV), placed in the world by its
XORIGIN, YORIGIN and ANGROT, often with inactive cells, and its own river, drain or SFR entries. These helpers write
the two maps the coupling reads (mapping.py, coupling/exchange.py) from SWAT+'s geometry on that grid, and report what
doesn't fit rather than fixing it:

- **HRUs (or aquifers) to cells** (``hru_map``): each polygon's ground in each column of the grid, on the column's top
  active layer, by area of overlap (flopy's GridIntersect). Ground outside the grid, or over columns with no active
  layer, is left out and reported by object, so the area check (``HruMap.problems``, and ``[recharge] area_rtol`` at
  start-up) says why an object's mapped area falls short of its area in SWAT+.
- **Channels to stream entries** (``channel_map``): each of the model's entries (``stream_entries``: a RIV or DRN
  package's cells over all its periods, or an SFR package's reaches) goes to the SWAT+ channel whose line runs
  furthest through its cell. An entry whose cell no channel crosses is reported as unreached, not given to the
  nearest channel, unless the study sets ``reach_m`` [A], a distance within which a channel reaches an entry; those
  are reported separately. A cell crossed by several channels goes to the longest and is reported as shared. Each
  channel's length is checked (``StreamMap.problems``): its line inside the model's active area against its length
  through the cells it was given, and, for SFR, against its reaches' own lengths.
- **SFR confluences** (``sfr_confluences``): the reaches more than one reach flows into, with their cells and
  places. Building SWAT+ on a model's SFR network with an outlet at each puts every confluence on a channel boundary,
  as the exchange's line-up rule asks.

Channel lines are drawn downstream and polygons are given in the grid's world coordinates, which are in the model's
length units (a grid in feet takes polygons and lines projected to a CRS in feet). Areas and lengths come out in
square meters and meters, converted by the grid's LENGTH_UNITS, or by ``length_units`` where the model leaves them
undefined, as ``[modflow6] length_units`` does for the run. Needs shapely (``pip install basinledger[gis]``). DISU
grids aren't supported.

    sim, gwf = load(workspace)
    hm = hru_map(gwf.modelgrid, hrus)                     # hrus: name -> shapely polygon
    print(hm.problems(swat_areas, rtol=1e-6)); hm.write("maps/hru.csv", gwf.modelgrid)
    sm = channel_map(gwf.modelgrid, stream_entries(gwf, "riv_creek"), lines)
    print(sm.report()); sm.write("maps/river.csv")
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


def load(workspace: str | Path, model: str | None = None):
    """A MODFLOW 6 simulation and its groundwater flow model, read with flopy (the only GWF model, or ``model``)."""
    import flopy
    sim = flopy.mf6.MFSimulation.load(sim_ws=str(workspace), verbosity_level=0)
    names = [m for m in sim.model_names if sim.get_model(m).model_type.startswith("gwf")]
    if model is None and len(names) != 1:
        raise ValueError(f"{workspace}: {len(names)} groundwater flow models {names}; name one")
    gwf = sim.get_model(model or names[0])
    if gwf.modelgrid.grid_type not in ("structured", "vertex"):
        raise ValueError(f"{workspace}: a {gwf.modelgrid.grid_type} grid; maps are for DIS and DISV grids")
    return sim, gwf


def meters_per_unit(grid, length_units: str | None = None) -> float:
    """Meters per unit of the grid's coordinates: its LENGTH_UNITS, or ``length_units`` where it has none. A unit
    stated that differs from the grid's own is refused."""
    from .coupling.modflow6 import METERS
    own = getattr(grid, "units", None)
    own = own if own in METERS else None
    if length_units is not None and length_units not in METERS:
        raise ValueError(f"length_units {length_units!r}: one of {', '.join(METERS)}")
    if own is None and length_units is None:
        raise ValueError("the grid's length units are undefined: give length_units (the model's, as [modflow6] "
                         "length_units states it)")
    if own is not None and length_units is not None and own != length_units:
        raise ValueError(f"length_units {length_units!r}, but the grid is in {own}")
    return METERS[own or length_units]


def _ncpl(grid) -> int:
    return grid.nrow * grid.ncol if grid.grid_type == "structured" else grid.ncpl


def top_active(grid) -> np.ndarray:
    """Per column (row-major on a DIS grid, icell2d on a DISV grid), the first layer with an active cell, -1 if none."""
    idomain = grid.idomain
    if idomain is None:
        return np.zeros(_ncpl(grid), int)
    act = np.asarray(idomain).reshape(grid.nlay, -1) > 0
    return np.where(act.any(axis=0), act.argmax(axis=0), -1)


def _column(grid, cell) -> int:
    """A cell id from GridIntersect or a package, without its layer: (row, col) or icell2d -> column index."""
    cell = tuple(int(c) for c in np.atleast_1d(cell))
    return cell[0] * grid.ncol + cell[1] if grid.grid_type == "structured" else cell[0]


def _cell_columns(grid, column: int) -> dict:
    return ({"row": column // grid.ncol, "col": column % grid.ncol} if grid.grid_type == "structured"
            else {"icell2d": column})


def _polygon(grid, column: int):
    import shapely
    if grid.grid_type == "structured":
        v = grid.get_cell_vertices(column // grid.ncol, column % grid.ncol)
    else:
        v = grid.get_cell_vertices(column)
    return shapely.Polygon(v)


def active_area(grid):
    """The model's active area in world coordinates: the union of the columns with an active cell."""
    import shapely
    return shapely.union_all([_polygon(grid, c) for c in np.flatnonzero(top_active(grid) >= 0)])


# ---------- HRUs to cells ----------

@dataclass
class HruMap:
    rows: list[tuple[str, int, int, float]]  # object, layer, column, area m2
    polygon_area: dict[str, float]           # each object's polygon, m2
    outside: dict[str, float]                # ground outside the grid, m2
    inactive: dict[str, float]               # ground over columns with no active layer, m2
    slivers: dict[str, float] = field(default_factory=dict)  # pieces under ``sliver_m2`` left out, m2
    grid_type: str = "structured"
    ncol: int = 0

    def mapped(self) -> dict[str, float]:
        out = dict.fromkeys(self.polygon_area, 0.0)
        for o, _, _, a in self.rows:
            out[o] += a
        return out

    def problems(self, source_area: dict[str, float] | None = None, rtol: float = 1e-6) -> list[str]:
        """Each object whose mapped area differs from its area in SWAT+ (``source_area``, m2; its polygon's area if
        not given) by more than ``rtol``, with where the difference went."""
        out = []
        for o, m in self.mapped().items():
            want = (source_area or self.polygon_area).get(o)
            if want is None:
                out.append(f"{o}: not an object of the source model")
                continue
            if abs(m - want) > rtol * want:
                why = [f"{v:,.0f} m2 {k}" for k, v in (("outside the grid", self.outside[o]),
                                                       ("over inactive columns", self.inactive[o]),
                                                       ("in slivers", self.slivers.get(o, 0.0))) if v > 0]
                if source_area is not None and abs(self.polygon_area[o] - want) > rtol * want:
                    why.append(f"its polygon is {self.polygon_area[o]:,.0f} m2")
                out.append(f"{o}: mapped {m:,.0f} m2 of {want:,.0f} m2" + (f" ({'; '.join(why)})" if why else ""))
        for o in sorted(set(source_area or {}) - set(self.polygon_area)):
            out.append(f"{o}: in the source model, but no polygon")
        return out

    def write(self, path: str | Path):
        """The map as mapping.py reads it: object, layer, row, col, area_m2 (DIS) or object, layer, icell2d, area_m2
        (DISV), 0-based."""
        cols = ["row", "col"] if self.grid_type == "structured" else ["icell2d"]
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["object", "layer", *cols, "area_m2"])
            for o, lay, c, a in self.rows:
                cell = [c // self.ncol, c % self.ncol] if self.grid_type == "structured" else [c]
                w.writerow([o, lay, *cell, f"{a:.6f}"])


def hru_map(grid, polygons: dict, sliver_m2: float = 1e-3, length_units: str | None = None) -> HruMap:
    """Each polygon's ground in each column of the grid, on the column's top active layer (see the module), in m2.
    Pieces smaller than ``sliver_m2`` are left out and reported (``HruMap.slivers``): where a polygon's edge runs
    along cell edges, as on a rotated grid, rounding leaves pieces of a millionth of a square meter in the next cells.
    ``length_units``: the grid's, where its model leaves them undefined."""
    import shapely
    from flopy.utils.gridintersect import GridIntersect

    m2 = meters_per_unit(grid, length_units) ** 2
    gi = GridIntersect(grid)
    top = top_active(grid)
    footprint = shapely.union_all([_polygon(grid, c) for c in range(_ncpl(grid))])
    rows, area, outside, inactive, slivers = [], {}, {}, {}, {}
    for o, poly in polygons.items():
        area[o] = float(poly.area) * m2
        outside[o] = float(poly.difference(footprint).area) * m2
        inactive[o] = slivers[o] = 0.0
        r = gi.intersect(poly, geo_dataframe=False)
        per: dict[int, float] = {}
        for cellid, a in zip(r["cellids"], r["areas"]):
            c = _column(grid, cellid)
            per[c] = per.get(c, 0.0) + float(a) * m2
        for c, a in sorted(per.items()):
            if a < sliver_m2:
                slivers[o] += a
            elif top[c] < 0:
                inactive[o] += a
            else:
                rows.append((o, int(top[c]), c, a))
    return HruMap(rows, area, outside, inactive, slivers, grid.grid_type,
                  getattr(grid, "ncol", 0) if grid.grid_type == "structured" else 0)


# ---------- channels to stream entries ----------

@dataclass
class Entry:
    cell: tuple[int, ...]       # (layer, row, col) or (layer, icell2d)
    reach: int | None = None    # SFR's 0-based reach number
    length: float | None = None  # SFR's reach length, m


def stream_entries(gwf, package: str, length_units: str | None = None) -> list[Entry]:
    """A package's stream entries: an SFR package's reaches in order (cell, reach number, length in m), or the
    distinct cells of a RIV, DRN or GHB package over all its stress periods. ``length_units``: the model's, where it
    leaves them undefined (needed for SFR's lengths only)."""
    pkg = gwf.get_package(package)
    if pkg is None:
        raise KeyError(f"no package {package!r} in {gwf.name}")
    ptype = pkg.package_type.lower()
    if ptype == "sfr":
        m = meters_per_unit(gwf.modelgrid, length_units)
        pd = pkg.packagedata.get_data()
        number = "ifno" if "ifno" in pd.dtype.names else "rno"
        out = []
        for r in sorted(pd, key=lambda x: int(x[number])):
            if not isinstance(r["cellid"], tuple) or any(c is None for c in r["cellid"]):
                continue  # a reach not connected to the aquifer (NONE) exchanges nothing
            out.append(Entry(tuple(int(c) for c in r["cellid"]), int(r[number]), float(r["rlen"]) * m))
        return out
    spd = pkg.stress_period_data.get_data()
    cells = []
    for per in sorted(spd or {}):
        if spd[per] is None:
            continue
        for c in spd[per]["cellid"]:
            cells.append(tuple(int(x) for x in c))
    return [Entry(c) for c in dict.fromkeys(cells)]


@dataclass
class StreamMap:
    rows: list[tuple[str, Entry, float]]       # channel, entry, the channel's length in the entry's cell, m
    unreached: list[Entry]                     # entries no channel crosses (or comes within reach_m of)
    by_distance: list[tuple[str, Entry, float]]  # channel, entry, distance m: reached only within reach_m
    shared: list[tuple[Entry, dict[str, float]]]  # entries whose cell several channels cross: their lengths, m
    in_grid: dict[str, float]                  # each channel's line inside the model's active area, m
    on_entries: dict[str, float] = field(default_factory=dict)   # its length through the cells it was given, m
    reach_length: dict[str, float] = field(default_factory=dict)  # SFR: its reaches' own lengths, m
    grid_type: str = "structured"

    @property
    def sfr(self) -> bool:
        return any(e.reach is not None for _, e, _ in self.rows) or any(e.reach is not None for e in self.unreached)

    def problems(self, rtol: float = 0.05) -> list[str]:
        """Channels whose length through their entries' cells differs from their line inside the model's active area
        by more than ``rtol`` (a channel the model has no stream entries along, for part of its length), and, for
        SFR, whose reaches' lengths differ from the line through their cells by more than ``rtol``."""
        out = []
        for ch, L in self.in_grid.items():
            got = self.on_entries.get(ch, 0.0)
            if L > 0 and abs(got - L) > rtol * L:
                out.append(f"{ch}: {L:,.0f} m of its line is inside the model, {got:,.0f} m of it through cells with "
                           "stream entries")
            if self.sfr and got > 0 and abs(self.reach_length.get(ch, 0.0) - got) > rtol * got:
                out.append(f"{ch}: its reaches are {self.reach_length.get(ch, 0.0):,.0f} m long, its line through "
                           f"their cells {got:,.0f} m")
        return out

    def report(self) -> str:
        """What the map did, for a reviewer: entries per channel, and every entry it didn't place on a line."""
        name = lambda e: f"reach {e.reach}" if e.reach is not None else f"cell {e.cell}"  # noqa: E731
        lines = []
        for ch in self.in_grid:
            n = sum(1 for c, _, _ in self.rows if c == ch)
            lines.append(f"{ch}: {n} entries, {self.on_entries.get(ch, 0.0):,.0f} m of its line through their cells "
                         f"({self.in_grid[ch]:,.0f} m inside the model)")
        if self.unreached:
            lines.append(f"{len(self.unreached)} entries no channel reaches (left as the model's): "
                         + ", ".join(name(e) for e in self.unreached))
        for ch, e, d in self.by_distance:
            lines.append(f"{name(e)}: to {ch}, {d:.1f} m from its line (reach_m)")
        for e, lengths in self.shared:
            lines.append(f"{name(e)}: crossed by " + ", ".join(f"{c} ({v:.1f} m)" for c, v in lengths.items())
                         + "; given to the longest")
        return "\n".join(lines)

    def write(self, path: str | Path):
        """The channel map as coupling/exchange.py reads it: object, reach, length_m for SFR; object, layer, row,
        col, length_m (DIS) or object, layer, icell2d, length_m (DISV) otherwise. Unreached entries aren't in it."""
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        cols = ["reach"] if self.sfr else ["layer", "row", "col"] if self.grid_type == "structured" else \
            ["layer", "icell2d"]
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["object", *cols, "length_m"])
            for ch, e, L in self.rows:
                w.writerow([ch, *([e.reach] if self.sfr else list(e.cell)), f"{L:.6f}"])


def channel_map(grid, entries: list[Entry], lines: dict, reach_m: float = 0.0,
                length_units: str | None = None) -> StreamMap:
    """Each entry to the channel whose line runs furthest through its cell (see the module). ``lines``: SWAT+ channel
    -> shapely LineString in the grid's world coordinates. ``reach_m`` [A]: a channel within this distance (m) of an
    entry's cell reaches it, when none crosses it; 0 (the default) guesses nothing. Lengths in m. ``length_units``:
    the grid's, where its model leaves them undefined."""
    m = meters_per_unit(grid, length_units)
    domain = active_area(grid)
    in_grid = {ch: float(line.intersection(domain).length) * m for ch, line in lines.items()}
    on = dict.fromkeys(lines, 0.0)
    rl = dict.fromkeys(lines, 0.0)
    rows, unreached, by_distance, shared = [], [], [], []
    seen_cell: set[tuple[str, int]] = set()
    for e in entries:
        col = _column(grid, e.cell[1:])
        poly = _polygon(grid, col)
        lengths = {ch: float(line.intersection(poly).length) * m for ch, line in lines.items()}
        lengths = {ch: v for ch, v in lengths.items() if v > 0}
        if lengths:
            ch = max(lengths, key=lengths.get)
            if len(lengths) > 1:
                shared.append((e, dict(sorted(lengths.items(), key=lambda kv: -kv[1]))))
            L = lengths[ch]
        else:
            near = {c: float(line.distance(poly)) * m for c, line in lines.items()}
            ch = min(near, key=near.get) if near else None
            if ch is None or reach_m <= 0 or near[ch] > reach_m:
                unreached.append(e)
                continue
            by_distance.append((ch, e, near[ch]))
            L = 0.0
        rows.append((ch, e, L))
        if (ch, col) not in seen_cell:  # several reaches or layers in one cell: its line counted once
            on[ch] += L
            seen_cell.add((ch, col))
        if e.length is not None:
            rl[ch] += e.length
    return StreamMap(rows, unreached, by_distance, shared, in_grid, on, rl, grid.grid_type)


# ---------- SFR confluences, where SWAT+ needs an outlet ----------

@dataclass
class Confluence:
    reach: int                    # the reach below the confluence, 0-based
    upstream: list[int]           # the reaches flowing into it, 0-based
    cell: tuple[int, ...] | None  # its cell, (layer, row, col) or (layer, icell2d); None for a reach off the grid
    xy: tuple[float, float] | None  # its cell's center, in the grid's world coordinates (its length units)


def sfr_links(gwf, package: str) -> list[tuple[int, int]]:
    """Each (upstream, downstream) pair of an SFR package's reaches, 0-based, read from its CONNECTIONDATA block as
    written (1-based, a negative number for a reach downstream), so the sign of reach 1 isn't lost."""
    pkg = gwf.get_package(package)
    if pkg is None or pkg.package_type.lower() != "sfr":
        raise KeyError(f"no SFR package {package!r} in {gwf.name}")
    path = Path(gwf.simulation_data.mfpath.get_sim_path()) / pkg.filename
    links, inside = set(), False
    for line in path.read_text().splitlines():
        words = line.split("#")[0].split()
        if not words:
            continue
        if words[0].upper() == "BEGIN" and words[1].upper() == "CONNECTIONDATA":
            inside = True
        elif words[0].upper() == "END" and inside:
            break
        elif inside:
            r = int(words[0]) - 1
            for w in words[1:]:
                c = int(float(w))
                links.add((r, -c - 1) if c < 0 else (c - 1, r))
    return sorted(links)


def sfr_aux(gwf, package: str, name: str) -> list:
    """An auxiliary value of each SFR reach, in reach order: SFRmaker writes ``line_id``, the hydrography line each
    reach was cut from. Where SWAT+'s channels are built on the same lines, it is the channel map (``object,
    reach``) once each line is given its channel's name."""
    pd = gwf.get_package(package).packagedata.get_data()
    if name not in pd.dtype.names:
        raise KeyError(f"{package}: no auxiliary {name!r} (it has {list(pd.dtype.names)})")
    number = "ifno" if "ifno" in pd.dtype.names else "rno"
    return [r[name] for r in sorted(pd, key=lambda x: int(x[number]))]


def sfr_confluences(gwf, package: str) -> list[Confluence]:
    """The reaches more than one reach flows into, with where they are: for a study building SWAT+ on a model's SFR
    network, a SWAT+ outlet at each one's top puts every confluence on a channel boundary, as the exchange's line-up
    rule asks (coupling/exchange.py, ``sfr_problems``)."""
    grid = gwf.modelgrid
    ups: dict[int, list[int]] = {}
    for u, d in sfr_links(gwf, package):
        ups.setdefault(d, []).append(u)
    cells = {e.reach: e.cell for e in stream_entries(gwf, package)}
    out = []
    for r in sorted(r for r, u in ups.items() if len(u) > 1):
        cell = cells.get(r)
        xy = None
        if cell is not None:
            col = _column(grid, cell[1:])
            xs, ys = np.asarray(grid.xcellcenters).ravel(), np.asarray(grid.ycellcenters).ravel()
            xy = (float(xs[col]), float(ys[col]))
        out.append(Confluence(r, sorted(ups[r]), cell, xy))
    return out
