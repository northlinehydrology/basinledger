"""SFR reaches for the two-way exchange, built from SWAT+'s channels (exchange.py, "Through SFR").

Basin Ledger builds and checks the links between models, not the models (docs/design.md, Scope). This module builds
what the SFR exchange depends on and nothing else: each SWAT+ channel's line cut into reaches on the MODFLOW 6 grid,
the reaches chained from upstream to downstream within their channel and never from one channel to the next (SWAT+
routes between channels), Manning's constant for meters and days, and the map of reaches to channels. The
streambed's numbers (top, width, thickness, conductivity, roughness) are the study's: they come in, and the helper
refuses a channel whose streambed rises downstream rather than smoothing it.

    from basinledger.coupling.sfr import channel_reaches, add_sfr, write_reach_map

    reaches = channel_reaches(gwf.modelgrid, lines)        # lines: SWAT+ channel -> shapely line, drawn downstream
    tops = [...]                                           # m, one per reach, from the study's DEM work [A]
    add_sfr(gwf, reaches, top=tops, width=4.0, thickness=1.0, k=0.5, manning=0.035)
    write_reach_map("maps/sfr.csv", reaches)               # [exchange] package = "sfr", map = "maps/sfr.csv"

``width``, ``thickness``, ``k``, ``manning`` and ``layer`` take one value, one per reach (in ``reaches``' order) or
one per channel (a dict). Works on structured (DIS) and vertex (DISV) grids, through flopy's GridIntersect; needs
shapely (``pip install basinledger[gis]``).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

DAY_S = 86400.0  # Manning's constant for meters and days (MODFLOW 6's UNIT_CONVERSION)


@dataclass(frozen=True)
class Reach:
    channel: str
    cell: tuple[int, ...]  # (row, col) on a structured grid, (icell2d,) on a vertex grid; the layer comes later
    length: float          # m
    distance: float        # m along the channel's line from its upstream end to the reach's middle


def channel_reaches(modelgrid, lines: dict, min_length: float = 0.0) -> list[Reach]:
    """Each channel's line cut by the grid's cells, one reach per piece of line in a cell (a line that leaves a cell
    and comes back makes two), in order downstream. ``lines``: SWAT+ channel -> shapely LineString drawn from
    upstream to downstream, in the grid's coordinates. Pieces no longer than ``min_length`` (m) are left out. A
    channel with no reach in the grid is refused."""
    import shapely
    from flopy.utils.gridintersect import GridIntersect

    gi = GridIntersect(modelgrid)
    out = []
    for ch, line in lines.items():
        if line.geom_type != "LineString":
            raise ValueError(f"channel {ch}: one line drawn downstream, not a {line.geom_type}")
        r = gi.intersect(line, geo_dataframe=False)
        mine = []
        for cellid, geom in zip(r["cellids"], r["ixshapes"]):
            cell = tuple(int(c) for c in np.atleast_1d(cellid))
            for part in shapely.get_parts(shapely.line_merge(geom) if geom.geom_type == "MultiLineString" else geom):
                if part.geom_type == "LineString" and part.length > min_length:
                    mid = line.project(part.interpolate(0.5, normalized=True))
                    mine.append(Reach(ch, cell, float(part.length), float(mid)))
        if not mine:
            raise ValueError(f"channel {ch}: its line crosses no cell of the grid")
        out += sorted(mine, key=lambda x: x.distance)
    return out


def _per_reach(name: str, value, reaches: list[Reach]) -> np.ndarray:
    if isinstance(value, dict):
        missing = sorted({r.channel for r in reaches} - set(value))
        if missing:
            raise ValueError(f"{name}: no value for channels {missing}")
        return np.array([value[r.channel] for r in reaches], dtype=float)
    a = np.broadcast_to(np.asarray(value, dtype=float), (len(reaches),)) if np.ndim(value) == 0 else \
        np.asarray(value, dtype=float)
    if a.shape != (len(reaches),):
        raise ValueError(f"{name}: {a.size} values for {len(reaches)} reaches")
    return a


def chains(reaches: list[Reach]) -> dict[str, list[int]]:
    """Each channel's reaches, by index, from upstream to downstream."""
    out: dict[str, list[int]] = {}
    for i, r in enumerate(reaches):
        out.setdefault(r.channel, []).append(i)
    for ch, idx in out.items():
        d = [reaches[i].distance for i in idx]
        if d != sorted(d):
            raise ValueError(f"channel {ch}: its reaches aren't in order downstream (channel_reaches gives them so)")
    return out


def slopes(reaches: list[Reach], top: np.ndarray, min_slope: float | None = None) -> np.ndarray:
    """Each reach's streambed slope from the tops: the fall to the next reach over the distance between their
    middles; a channel's last reach takes its upstream neighbor's. A streambed rising downstream is refused; a flat
    one too, unless ``min_slope`` [A] stands in for it. A channel of one reach needs its slope given."""
    out = np.empty(len(reaches))
    for ch, idx in chains(reaches).items():
        if len(idx) == 1:
            raise ValueError(f"channel {ch}: one reach, so its slope can't come from the tops; give slope")
        for a, b in zip(idx[:-1], idx[1:]):
            fall = top[a] - top[b]
            if fall < 0:
                raise ValueError(f"channel {ch}: the streambed rises {-fall:.3f} m downstream from reach {a} "
                                 f"({top[a]:.3f} m) to reach {b} ({top[b]:.3f} m); the tops are the study's to "
                                 "correct (culverts, ponds and bridges show in a DEM)")
            out[a] = fall / (reaches[b].distance - reaches[a].distance)
        out[idx[-1]] = out[idx[-2]]
        flat = [i for i in idx if out[i] <= 0]
        if flat:
            if min_slope is None:
                raise ValueError(f"channel {ch}: reaches {flat} are flat, and SFR needs a slope; give min_slope [A] "
                                 "or slope")
            out[flat] = min_slope
    return out


def add_sfr(gwf, reaches: list[Reach], *, top, width, thickness, k, manning, layer=0, slope=None,
            min_slope: float | None = None, pname: str = "sfr", **options):
    """An SFR package on ``gwf`` (a flopy GWF model) with a reach for each of ``reaches``, ready for the exchange:
    reaches chained within each channel only, Manning's constant for meters and days, and no inflow of its own (the
    exchange sets each channel's top reach every day). ``top`` (m, streambed top), ``width`` (m), ``thickness`` (m,
    the bed's), ``k`` (m d-1, the bed's vertical hydraulic conductivity), ``manning`` and ``layer`` (0-based) are
    the study's [A]. ``slope``: per reach, or from the tops (``slopes``). ``options`` go to flopy's ModflowGwfsfr
    (e.g. ``save_flows``); setting another Manning's constant is refused. The model must use meters and days."""
    import flopy

    if gwf.modelgrid.lenuni != 2 or str(gwf.simulation.tdis.time_units.get_data()).lower() != "days":
        raise ValueError("SFR for the exchange: the model must use meters and days, as Basin Ledger needs")
    for key in ("unit_conversion", "length_conversion", "time_conversion"):
        if key in options:
            raise ValueError(f"SFR for the exchange: {key} is set here, for meters and days")
    top = _per_reach("top", top, reaches)
    vals = {n: _per_reach(n, v, reaches) for n, v in
            (("width", width), ("thickness", thickness), ("k", k), ("manning", manning))}
    for n, v in vals.items():
        if (v <= 0).any():
            raise ValueError(f"{n}: every reach's must be above 0")
    lay = _per_reach("layer", layer, reaches).astype(int)
    rgrd = _per_reach("slope", slope, reaches) if slope is not None else slopes(reaches, top, min_slope)
    conn = [[i] for i in range(len(reaches))]
    for idx in chains(reaches).values():
        for a, b in zip(idx[:-1], idx[1:]):
            conn[a].append(-b)   # b is downstream of a
            conn[b].append(a)
    pd = [(i, (int(lay[i]), *r.cell), r.length, vals["width"][i], float(rgrd[i]), top[i], vals["thickness"][i],
           vals["k"][i], vals["manning"][i], len(conn[i]) - 1, 1.0, 0) for i, r in enumerate(reaches)]
    return flopy.mf6.ModflowGwfsfr(gwf, pname=pname, nreaches=len(reaches), packagedata=pd, connectiondata=conn,
                                   unit_conversion=DAY_S, **options)


def write_reach_map(path: str | Path, reaches: list[Reach]):
    """The exchange's map of reaches to channels: object, reach (0-based, in ``reaches``' order)."""
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["object", "reach"])
        for i, r in enumerate(reaches):
            w.writerow([r.channel, i])
