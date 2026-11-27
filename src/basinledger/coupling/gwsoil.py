"""Groundwater up into the soil: MODFLOW 6's water to SWAT+'s HRU soil layers (docs/design.md, Groundwater up into the
soil). Built Nov 8-9, 2026 on swatplus-bmi 0.7.0's ``soil_layer_water__groundwater_inflow_rate``. Three rules for how
much water leaves each cell (``rule = "drain"``, chosen Nov 9, ``"gwflow"``, or ``"water_table"``, Nov 27: a drain at
MODFLOW 6's top, which is the profile's bottom, and SWAT+'s layers below the water table held, docs/design.md, A
water-table link) and an optional cap (``cap = "none"``, chosen, or ``"room"``), both required in a config. The
pieces:

- ``cell_shares``: each HRU's share of a cell's mapped ground, which splits a cell's water among the HRUs over it;
- ``to_layers``: water from each (cell, HRU) row into that HRU's layers below the cell's water table, each in
  proportion to its thickness below it (gwflow's split);
- ``drain_elevation``: the bottom of the soil profile under each cell, MODFLOW 6's top less the profile depths of the
  HRUs over it, by their shares, for a drain at the profile's bottom;
- ``hru_table_height``, ``held_layers``: the water table's height above each HRU's profile bottom (MODFLOW 6's heads
  above its top, by the HRU's ground in each cell), and which of its layers are held below it (``water_table``);
- ``top_problems``: MODFLOW 6's top checked against each HRU's land surface less its profile depth (``water_table``).

Soil layers are given flat, as the hly grid has them: ``layer_hru`` (each layer's HRU index) and ``layer_bottom``
(each layer's bottom depth below the surface, m).
"""

from __future__ import annotations

import numpy as np

from ..mapping import AreaMap


def cell_shares(amap: AreaMap, n_cells: int) -> np.ndarray:
    """Per map row, its HRU's share of the cell's mapped ground; a cell's shares add to 1."""
    covered = np.bincount(amap.node, weights=amap.area, minlength=n_cells)
    return amap.area / covered[amap.node]


def profile_depth(layer_hru: np.ndarray, layer_bottom: np.ndarray, n_hru: int) -> np.ndarray:
    """Each HRU's profile depth (its deepest layer's bottom), m."""
    out = np.zeros(n_hru)
    np.maximum.at(out, layer_hru, layer_bottom)
    return out


def drain_elevation(amap: AreaMap, top: np.ndarray, depth: np.ndarray, n_cells: int) -> np.ndarray:
    """Per cell, MODFLOW 6's top less the profile depth of the HRUs over it, weighted by their shares of the cell's
    ground; NaN where no HRU is mapped."""
    shares = cell_shares(amap, n_cells)
    d = np.bincount(amap.node, weights=shares * depth[amap.obj], minlength=n_cells)
    out = top - d
    out[np.bincount(amap.node, minlength=n_cells) == 0] = np.nan
    return out


def layer_split(top_depth: np.ndarray, bottom_depth: np.ndarray, water_depth: float) -> np.ndarray:
    """One HRU's layers (tops and bottoms, m below the surface): the share of water each takes when the water table
    is ``water_depth`` below the surface, by its thickness below the water table. All zero with the water table
    below the profile."""
    below = np.clip(bottom_depth - np.maximum(top_depth, water_depth), 0.0, None)
    total = below.sum()
    return below / total if total > 0 else np.zeros_like(below)


def to_layers(row_volume: np.ndarray, row_water_depth: np.ndarray, amap: AreaMap, layer_hru: np.ndarray,
              layer_bottom: np.ndarray) -> tuple[np.ndarray, float]:
    """Volumes (m3) per map row into the HRU's soil layers, each row split at its own cell's depth to water
    (``layer_split``'s rule, for every row at once). Returns m3 per layer, and the volume of rows whose water table is
    below their HRU's profile: that water still left the cell (a drain under a cell is at the deepest profile over
    it, by share), so it goes to the HRU's bottom layer, where it reaches the soil. Vectorized Nov 13, 2026: the loop
    over rows took half of a coupled day on a study basin (2,371 HRUs, 11,400 rows); ``to_layers_loop`` is the
    reference it is tested against."""
    layer_hru = np.asarray(layer_hru)
    layer_bottom = np.asarray(layer_bottom, float)
    row_volume = np.asarray(row_volume, float)
    out = np.zeros(len(layer_bottom))
    nz = np.flatnonzero(row_volume != 0.0)
    if not len(nz):
        return out, 0.0
    # each HRU's layers in their order, as a padded table (HRU x position) of layer indices, -1 past its last
    srt = np.argsort(layer_hru, kind="stable")
    hs = layer_hru[srt]
    _, first, counts = np.unique(hs, return_index=True, return_counts=True)
    rank = np.arange(len(srt)) - np.repeat(first, counts)
    table = np.full((int(layer_hru.max()) + 1, int(counts.max())), -1)
    table[hs, rank] = srt
    tops = np.empty_like(layer_bottom)
    prev = np.concatenate([[0.0], layer_bottom[srt][:-1]])
    tops[srt] = np.where(rank == 0, 0.0, prev)
    lay = table[np.asarray(amap.obj)[nz]]                 # rows x positions
    ok = lay >= 0
    if not ok.any(axis=1).all():
        raise ValueError("a mapped HRU has no soil layers")
    li = np.where(ok, lay, 0)
    wd = np.asarray(row_water_depth, float)[nz][:, None]
    below = np.where(ok, np.clip(layer_bottom[li] - np.maximum(tops[li], wd), 0.0, None), 0.0)
    total = below.sum(axis=1)
    v = row_volume[nz]
    deep = total <= 0
    share = np.divide(below, total[:, None], out=np.zeros_like(below), where=~deep[:, None])
    np.add.at(out, li[ok & ~deep[:, None]], (v[:, None] * share)[ok & ~deep[:, None]])
    if deep.any():  # the water table below the profile: the HRU's bottom layer
        bottom = np.where(ok, layer_bottom[li], -np.inf)
        np.add.at(out, li[deep, np.argmax(bottom[deep], axis=1)], v[deep])
    return out, float(v[deep].sum())


def to_layers_loop(row_volume: np.ndarray, row_water_depth: np.ndarray, amap: AreaMap, layer_hru: np.ndarray,
                   layer_bottom: np.ndarray) -> tuple[np.ndarray, float]:
    """``to_layers`` one row at a time, as built Nov 8-9, 2026: the reference the vectorized one is tested against."""
    layer_hru = np.asarray(layer_hru)
    layer_bottom = np.asarray(layer_bottom, float)
    order = {h: np.flatnonzero(layer_hru == h) for h in np.unique(layer_hru)}
    tops = np.empty_like(layer_bottom)
    for idx in order.values():
        b = layer_bottom[idx]
        tops[idx] = np.concatenate([[0.0], b[:-1]])
    out = np.zeros(len(layer_bottom))
    at_bottom = 0.0
    for v, wd, h in zip(row_volume, row_water_depth, amap.obj):
        if v == 0.0:
            continue
        idx = order[h]
        s = layer_split(tops[idx], layer_bottom[idx], wd)
        if s.sum() == 0:
            out[idx[np.argmax(layer_bottom[idx])]] += v
            at_bottom += float(v)
            continue
        out[idx] += v * s
    return out, at_bottom


def layer_rates(layer_volume: np.ndarray, layer_hru: np.ndarray, hru_area: np.ndarray, day_s: float = 86400.0):
    """m3 per layer over a day -> m s-1 over each layer's HRU area, as swatplus-bmi takes it."""
    return layer_volume / hru_area[np.asarray(layer_hru)] / day_s


def layer_tops(layer_hru: np.ndarray, layer_bottom: np.ndarray) -> np.ndarray:
    """Each layer's top depth (m below the surface): the bottom of the layer above it in the same HRU, 0 for the
    first. Layers are listed HRU by HRU, top down, as the hly grid has them."""
    layer_hru = np.asarray(layer_hru)
    layer_bottom = np.asarray(layer_bottom, float)
    tops = np.concatenate([[0.0], layer_bottom[:-1]])
    first = np.concatenate([[True], layer_hru[1:] != layer_hru[:-1]])
    tops[first] = 0.0
    return tops


def hru_table_height(amap: AreaMap, head: np.ndarray, top: np.ndarray, n_hru: int) -> np.ndarray:
    """Per HRU, the water table's height (m) above its profile's bottom, MODFLOW 6's top under ``water_table``: each
    cell's head above its top, weighted by the HRU's ground in it (the map before ``conserving``). Negative where the
    table is below the profile. A dry or inactive cell (no finite head below 1e20) counts as far below. NaN for an HRU
    the map doesn't place."""
    h = np.asarray(head, float)
    rel = np.where(np.isfinite(h) & (np.abs(h) < 1e20), h - top, -1e6)[amap.node]
    w = np.bincount(amap.obj, weights=amap.area, minlength=n_hru)
    s = np.bincount(amap.obj, weights=rel * amap.area, minlength=n_hru)
    return np.divide(s, w, out=np.full(n_hru, np.nan), where=w > 0)


def held_layers(height: np.ndarray, depth: np.ndarray, layer_hru: np.ndarray, layer_top: np.ndarray,
                layer_bottom: np.ndarray, before: np.ndarray | None = None) -> np.ndarray:
    """Which layers are held (1) under each HRU's water table (Andrew's decision 4, Nov 27, 2026): held once the table
    rises above the layer's middle, released once it falls below the layer's bottom, and as before in between (that
    half-layer gap keeps a table sitting mid-layer from flipping it day to day). With no ``before``, held above the
    middle. A layer is held only when every layer below it in the HRU is, so the held zone always runs up from the
    profile's bottom, as swatplus-bmi requires. ``height``: the table above each HRU's profile bottom (m);
    ``depth``: each HRU's profile depth (m)."""
    layer_hru = np.asarray(layer_hru)
    height = np.nan_to_num(np.asarray(height, float), nan=-1e6)
    table = (np.asarray(depth, float) - height)[layer_hru]  # m below the surface
    mid = 0.5 * (np.asarray(layer_top, float) + np.asarray(layer_bottom, float))
    held = table < mid
    if before is not None:
        keep = (table >= mid) & (table <= layer_bottom)
        held = np.where(keep, np.asarray(before) > 0.5, held)
    # from the bottom up: a layer held only if every layer below it in its HRU is
    out = held.copy()
    for i in range(len(out) - 2, -1, -1):
        if layer_hru[i + 1] == layer_hru[i] and not out[i + 1]:
            out[i] = False
    return out.astype(float)


def top_problems(amap: AreaMap, top: np.ndarray, land: np.ndarray, depth: np.ndarray, hrus: list[str],
                 tolerance: float) -> list[str]:
    """``water_table`` needs MODFLOW 6's top at the profile's bottom (Andrew's decision 3, Nov 27, 2026): per HRU, the
    mean of its cells' tops, weighted by its ground in each, against its land surface (hru.con's ``elev``) less its
    profile depth. The HRUs more than ``tolerance`` (m) off, worst first."""
    n = len(hrus)
    w = np.bincount(amap.obj, weights=amap.area, minlength=n)
    mean_top = np.divide(np.bincount(amap.obj, weights=top[amap.node] * amap.area, minlength=n), w,
                         out=np.full(n, np.nan), where=w > 0)
    off = mean_top - (np.asarray(land, float) - np.asarray(depth, float))
    bad = np.flatnonzero(np.abs(np.nan_to_num(off)) > tolerance)
    bad = bad[np.argsort(-np.abs(off[bad]))]
    return [f"{hrus[i]}: MODFLOW 6's top {mean_top[i]:.2f} m, its land surface less its profile "
            f"{land[i] - depth[i]:.2f} m ({off[i]:+.2f} m)" for i in bad]


def read_land_surface(path, node_of, n_cells: int) -> np.ndarray:
    """Each cell's land surface (m) from a CSV: a ``land_m`` column and the cell as a map table gives it (``node``, or
    ``layer, row, col`` on DIS, or ``layer, icell2d`` on DISV; ``mapping.cell_node``). NaN where not given."""
    import csv

    from ..mapping import cell_node
    out = np.full(n_cells, np.nan)
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            out[cell_node(r, node_of, str(path))] = float(r["land_m"])
    return out


def profile_bottom_cells(amap: AreaMap, top: np.ndarray, land: np.ndarray, depth: np.ndarray, n_cells: int,
                         tolerance: float) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Under ``water_table`` with each cell's land surface given: the profile's bottom per cell (its land surface less
    the profile depths of the HRUs over it, by ground share). A cell whose top is that bottom (within ``tolerance``,
    m) is a soil-zone cell; one whose top is its land surface (a stream's or lake's cell, which keeps the DEM top) has
    its drain and table a profile below its top. Returns the bottom per cell (NaN off the map), which cells keep the
    land surface as their top, and the cells that are neither, worst first."""
    shares = cell_shares(amap, n_cells)
    d = np.bincount(amap.node, weights=shares * depth[amap.obj], minlength=n_cells)
    bottom = land - d
    mapped = np.bincount(amap.node, minlength=n_cells) > 0
    bottom[~mapped] = np.nan
    to_bottom, to_land = np.abs(top - bottom), np.abs(top - land)
    at_bottom = (to_bottom <= tolerance) & ~(to_land < to_bottom)  # the nearer of the two, within the tolerance
    at_land = (to_land <= tolerance) & (to_land < to_bottom)
    missing = mapped & ~np.isfinite(land)
    bad = np.flatnonzero(mapped & ~missing & ~at_bottom & ~at_land)
    bad = bad[np.argsort(-np.abs(top[bad] - bottom[bad]))]
    msgs = [f"node {i}: no land surface given" for i in np.flatnonzero(missing)]
    msgs += [f"node {i}: top {top[i]:.2f} m, land surface {land[i]:.2f} m, profile's bottom {bottom[i]:.2f} m"
             for i in bad]
    return bottom, mapped & at_land, msgs


def looks_like_profile_bottom(amap: AreaMap, top: np.ndarray, land: np.ndarray, depth: np.ndarray) -> float:
    """The share of the mapped ground whose HRU's mean top is within half its profile depth of its land surface less
    its profile (hru.con's ``elev``), and nearer that than its land surface: a model built for ``water_table``, which
    the ``drain`` and ``gwflow`` rules would drain a second profile down."""
    n = len(land)
    w = np.bincount(amap.obj, weights=amap.area, minlength=n)
    mean_top = np.divide(np.bincount(amap.obj, weights=top[amap.node] * amap.area, minlength=n), w,
                         out=np.full(n, np.nan), where=w > 0)
    off = np.abs(mean_top - (land - depth))
    near = (off < 0.5 * np.asarray(depth)) & (off < np.abs(mean_top - land))
    return float(w[near & (w > 0)].sum() / w.sum()) if w.sum() > 0 else 0.0


def split_back(rows: np.ndarray, amap: AreaMap, hru_volume: np.ndarray, n_cells: int) -> np.ndarray:
    """Per cell (m3): each HRU's volume (``hru_volume``, m3) spread over the map rows its water came from, by each row's
    part of it (``rows``, m3 per map row). An HRU with nothing from its rows spreads it by its ground."""
    obj = np.asarray(amap.obj)
    n_hru = len(hru_volume)
    tot = np.bincount(obj, weights=rows, minlength=n_hru)
    ground = np.bincount(obj, weights=amap.area, minlength=n_hru)
    share = np.where(tot[obj] > 0, np.divide(rows, tot[obj], out=np.zeros(len(rows)), where=tot[obj] > 0),
                     np.divide(amap.area, ground[obj], out=np.zeros(len(rows)), where=ground[obj] > 0))
    return np.bincount(amap.node, weights=np.asarray(hru_volume, float)[obj] * share, minlength=n_cells)


# ---------- the link ----------

INFLOW = "soil_layer_water__groundwater_inflow_rate"  # per HRU soil layer, m s-1, set for the next day (0.7.0)
HRU_INFLOW = "soil_water__groundwater_inflow_rate"    # per HRU, m s-1, gwsoilq
BOTTOM = "soil_layer_bottom__depth"                   # m below the surface
WATER = "soil_layer_water__depth"                     # m above the wilting point (st)
SATURATION = "soil_layer_water__saturation_depth"     # m above the wilting point (ul)
# water_table (docs/design.md, A water-table link; swatplus-bmi's hooks)
FLAG = "soil_layer__saturated_zone_flag"              # per HRU soil layer, 1 held below the water table, a state
NOT_TAKEN = "soil_water__groundwater_inflow_not_taken_rate"  # per HRU, m s-1: what no layer could take at Hook 1
SATURATION_EXCESS = "land_surface_water__saturation_excess_runoff_rate"  # per HRU, m s-1, part of its runoff
RULES = ("drain", "gwflow", "water_table")
CAPS = ("none", "room")
HELD_LATERAL = ("none", "swatplus")
# held_lateral is swatplus-bmi's own option, given to SWAT+ at its start ("none" or "swat"), not a variable
SWATPLUS_BMI_HELD_LATERAL = {"none": "none", "swatplus": "swat"}


def takes_start_file() -> bool:
    """Whether the installed swatplus-bmi takes a TOML start file in place of the project folder (0.17.0, with the
    held soil layers)."""
    import re

    import swatplus_bmi
    m = re.match(r"(\d+)\.(\d+)", getattr(swatplus_bmi, "__version__", "0.0"))
    return m is not None and (int(m.group(1)), int(m.group(2))) >= (0, 17)


def write_start_file(path, project, held_lateral: str):
    """swatplus-bmi's start file for ``rule = "water_table"``: the project and ``held_lateral`` (Basin Ledger's
    "none" or "swatplus", swatplus-bmi's "none" or "swat"). Paths written POSIX-style, so Windows' backslashes aren't
    read as TOML escapes."""
    from pathlib import Path
    path = Path(path)
    path.write_text(f'project = "{Path(project).resolve().as_posix()}"\n'
                    f'held_lateral = "{SWATPLUS_BMI_HELD_LATERAL[held_lateral]}"\n')
    return path
DAY = 86400.0


class GroundwaterToSoil:
    """MODFLOW 6's water up into SWAT+'s HRU soil layers, through one MODFLOW 6 package with an entry on every mapped
    cell (``[groundwater_to_soil]``; docs/design.md, Groundwater up into the soil).

    Each day, once MODFLOW 6 has prepared it (``before_solve``), the package is set for the day: ``rule = "drain"``,
    a drain (DRN) whose elevation is the bottom of the soil profile under the cell, MODFLOW 6 computing the outflow
    inside its solve at the conductance in its input [A]; ``rule = "gwflow"``, SWAT+ gwflow's rule, a well (WEL)
    taking specific yield x the water table's height above the profile's bottom x each HRU's ground in the cell, at
    the day's starting heads. With ``cap = "room"`` no cell gives more than the room left (saturation less water, as
    SWAT+ had it at the end of its day) in the layers of its HRUs below its water table, at the day's starting heads:
    the drain's elevation is raised so its outflow at those heads fits (a head that rises during the day takes a
    little more, which SWAT+ sends to the surface), or the well's rate is cut. The drain isn't moved within the solve:
    moving its elevations between outer iterations stopped MODFLOW 6 converging, and aborted it on some runs. After the solve
    (``after_solve``), what left each cell is split among the HRUs over it by their ground and among each HRU's
    layers below the water table, and set in SWAT+ the next day (``before_swat``). What SWAT+ reports went in
    (``received``) is what was due."""

    def __init__(self, sw, amap: AreaMap, hru_area: np.ndarray, rule: str, cap: str, package: str,
                 due: np.ndarray | None = None, held_lateral: str = "none", land: np.ndarray | None = None,
                 top_tolerance: float = 1.0, saved: dict | None = None, n_cells: int | None = None,
                 land_cells: np.ndarray | None = None, ponded: np.ndarray | None = None):
        """SWAT+'s side, before either model steps; ``attach`` adds MODFLOW 6's once its package has its entries
        (after the first ``prepare_day``). ``amap``: each HRU's ground in each cell, before ``conserving``.
        ``due``: what a saved run had due for SWAT+'s next day. ``water_table`` only: ``held_lateral``, ``land`` (each
        HRU's land surface, m, hru.con's ``elev``) and ``top_tolerance`` (m) for the start-up check of MODFLOW 6's top,
        ``saved`` (a saved run's flags and the rows its due came from) and ``n_cells`` (MODFLOW 6's nodes).
        ``land_cells``: each cell's land surface (m), ``water_table``'s check cell by cell in place of hru.con's.
        ``ponded``: per HRU, whether it is a paddy or wetland (``swatplus.ponded_hrus``), whose saturation excess
        SWAT+ puts into its ponded water, not its runoff."""
        if rule not in RULES:
            raise ValueError(f"[groundwater_to_soil] rule must be one of {RULES}")
        if cap not in CAPS:
            raise ValueError(f"[groundwater_to_soil] cap must be one of {CAPS}")
        if not sw.has(INFLOW):
            raise ValueError(f"[groundwater_to_soil] needs swatplus-bmi's {INFLOW} (0.7.0)")
        if rule == "water_table":
            if cap != "none":
                raise ValueError('[groundwater_to_soil] rule = "water_table" takes cap = "none": the held layers '
                                 "take what they have room for, and what no layer takes goes back as recharge")
            if held_lateral not in HELD_LATERAL:
                raise ValueError(f"[groundwater_to_soil] held_lateral must be one of {HELD_LATERAL}")
            missing = [v for v in (FLAG, NOT_TAKEN, SATURATION_EXCESS) if not sw.has(v)]
            if missing:
                raise ValueError(f'[groundwater_to_soil] rule = "water_table" needs swatplus-bmi\'s held soil layers '
                                 f"({', '.join(missing)}): a newer swatplus-bmi")
        self.sw, self.map, self.rule, self.cap, self.package = sw, amap, rule, cap, package
        self.mf6 = None
        self.hru_area = np.asarray(hru_area, float)
        hrus = list(amap.objects)
        index = {h: i for i, h in enumerate(hrus)}
        names = sw.names(INFLOW)
        self.layer_hru = np.array([index[n.split(":")[0]] for n in names])
        self.bottom = sw.get(BOTTOM).astype(float)
        self.depth = profile_depth(self.layer_hru, self.bottom, len(hrus))
        self.due = np.zeros(len(names)) if due is None else np.asarray(due, float)  # m s-1 per layer, SWAT+'s next day
        self.raised = 0
        self.held_lateral, self.land, self.top_tolerance = held_lateral, land, top_tolerance
        self.land_cells = None if land_cells is None else np.asarray(land_cells, float)
        self.ponded = np.zeros(len(hrus), bool) if ponded is None else np.asarray(ponded, bool)
        self.notes = []  # what attach found worth saying (a warning, or the cells keeping the land surface)
        self.n_cells = int(np.max(amap.node)) + 1 if n_cells is None else int(n_cells)
        self.not_taken = 0.0
        if rule == "water_table":
            self.hrus = hrus
            self.top_depth = layer_tops(self.layer_hru, self.bottom)
            saved = saved or {}
            self.flags = None if saved.get("flags") is None else np.asarray(saved["flags"], float)
            self.rows = None if saved.get("rows") is None else np.asarray(saved["rows"], float)
            self.returned = None  # per cell, m3: what SWAT+ didn't take today, recharge in today's solve

    def attach(self, mf6):
        n = mf6.nodes
        self.mf6 = mf6
        self.shares = cell_shares(self.map, n)
        self.top = mf6.node_values("TOP", "DIS")
        nodes = mf6.package_nodes(self.package)
        mapped = np.unique(self.map.node)
        if sorted(nodes.tolist()) != mapped.tolist():
            raise ValueError(f"MODFLOW 6's {self.package} package must have one entry on each mapped cell and no "
                             f"other ({len(nodes)} entries, {len(mapped)} mapped cells)")
        self.nodes = nodes
        if self.rule == "water_table":  # MODFLOW 6's top is the profile's bottom: the drain sits at the top
            need = (f'[groundwater_to_soil] rule = "water_table" needs MODFLOW 6\'s top at the soil profile\'s '
                    f"bottom (its land surface less its profile depth, within top_tolerance = {self.top_tolerance} m)")
            self.bottom_elev = self.top.copy()  # the profile's bottom per cell: the drain's elevation, the table's zero
            if self.land_cells is not None:  # cell by cell; a stream's or lake's cell may keep its land surface
                self.bottom_elev, kept, bad = profile_bottom_cells(self.map, self.top, self.land_cells, self.depth,
                                                                   n, self.top_tolerance)
                if bad:
                    raise ValueError(f"{need}, or at its land surface (a stream's or lake's cell); {len(bad)} cells "
                                     "are neither: " + "; ".join(bad[:5]))
                if kept.any():
                    self.notes.append(f"{int(kept.sum())} cells keep the land surface as MODFLOW 6's top (streams' "
                                      "and lakes' cells): their drain and water table are a profile below the top, "
                                      "and their soil zone is in both models")
            elif self.land is not None:  # HRU by HRU, from hru.con's elev
                bad = top_problems(self.map, self.top, self.land, self.depth, self.hrus, self.top_tolerance)
                if bad:
                    raise ValueError(f"{need}; {len(bad)} HRUs are not: " + "; ".join(bad[:5]))
            self.elev0 = self.bottom_elev[nodes].copy()
            if self.rows is None:
                self.rows = np.zeros(len(self.map.obj))
        else:
            self.elev0 = drain_elevation(self.map, self.top, self.depth, n)[nodes]
            if self.land is not None:
                share = looks_like_profile_bottom(self.map, self.top, self.land, self.depth)
                if share > 0.5:
                    self.notes.append(
                        f"MODFLOW 6's top looks like the soil profile's bottom under {share:.0%} of the mapped "
                        f"ground (hru.con's elev less the profile), a model built for rule = \"water_table\": "
                        f'rule = "{self.rule}" puts its soil link a profile below the top, a second profile down')
        if self.rule == "gwflow":
            self.sy = mf6.node_values("SY", "STO")
        self.planned = np.zeros(len(nodes))

    # -- SWAT+'s side --
    def due_volume(self) -> float:
        return float(np.sum(self.due * DAY * self.hru_area[self.layer_hru]))

    def before_swat(self) -> float:
        """Set what's due in SWAT+'s layers for its coming day; returns its volume, m3. With ``water_table``, the
        held layers too (from yesterday's heads)."""
        if self.rule == "water_table":
            if self.flags is not None:
                self.sw.set(FLAG, self.flags)
        self.sw.set(INFLOW, self.due)
        return self.due_volume()

    def held_count(self) -> int:
        """The layers held for SWAT+'s coming day."""
        return 0 if getattr(self, "flags", None) is None else int(np.sum(self.flags > 0.5))

    def received(self) -> float:
        """What SWAT+ reports went into its layers on the day just done, m3. With ``water_table``, what no layer could
        take (``not_taken``) is spread back over the cells it came from, for today's recharge (``returned``)."""
        got = float(np.sum(self.sw.get(INFLOW) * DAY * self.hru_area[self.layer_hru]))
        if self.rule == "water_table":
            nt = np.clip(self.sw.get(NOT_TAKEN).astype(float), 0.0, None) * DAY * self.hru_area
            self.not_taken = float(nt.sum())
            rows = self.rows if self.rows is not None else np.zeros(len(self.map.obj))
            self.returned = split_back(rows, self.map, nt, self.n_cells) if self.not_taken > 0 else None
        return got

    def saturation_excess(self) -> tuple[float, float]:
        """What SWAT+'s full profiles spilled at the surface on the day just done (``swr_satexcess``), m3: on most
        HRUs a part of its surface runoff, not on top of it; on a paddy or wetland, into its ponded water. Returns
        the two."""
        v = self.sw.get(SATURATION_EXCESS).astype(float) * DAY * self.hru_area
        return float(v[~self.ponded].sum()), float(v[self.ponded].sum())

    def state(self) -> dict:
        """What a restart needs beyond ``due``: the held layers, and the rows the due came from."""
        if self.rule != "water_table":
            return {}
        return {"flags": None if self.flags is None else self.flags.tolist(),
                "rows": None if self.rows is None else np.asarray(self.rows).tolist()}

    # -- MODFLOW 6's side --
    def _room_rows(self, head: np.ndarray) -> np.ndarray:
        """Per map row, the room (m3) in its HRU's layers below the cell's water table, at SWAT+'s end of day."""
        room = np.clip(self.sw.get(SATURATION) - self.sw.get(WATER), 0.0, None)
        wd = self.top - head
        out = np.zeros(len(self.map.obj))
        for r, (h, node) in enumerate(zip(self.map.obj, self.map.node)):
            idx = np.flatnonzero(self.layer_hru == h)
            tops = np.concatenate([[0.0], self.bottom[idx][:-1]])
            below = np.clip(self.bottom[idx] - np.maximum(tops, wd[node]), 0.0, None)
            thick = self.bottom[idx] - tops
            out[r] = self.map.area[r] * float(np.sum(room[idx] * np.divide(below, thick, out=np.zeros_like(below),
                                                                              where=thick > 0)))
        return out

    def _cell_cap(self, head: np.ndarray) -> np.ndarray:
        return np.bincount(self.map.node, weights=self._room_rows(head), minlength=self.mf6.nodes)[self.nodes]

    def before_solve(self):
        """Set the package for the day, after ``prepare_day`` and before the solve."""
        mf = self.mf6
        h0 = mf.head()
        self.h0 = h0
        if self.rule in ("drain", "water_table"):  # water_table: a drain at the top, never raised
            elev = self.elev0.copy()
            self.raised = 0
            if self.cap == "room":
                # the drain raised so it gives no more than the room at the day's starting heads: q = C (h - elev)
                cond = mf.entry_values("COND", self.package)
                cap = self._cell_cap(h0)
                want = h0[self.nodes] - cap / np.where(cond > 0, cond, np.inf)
                elev = np.maximum(self.elev0, want)
                self.raised = int(np.sum(elev > self.elev0))
            mf.set_entries("ELEV", self.package, elev)
            return
        # gwflow: Sy x the water table's height above each HRU's profile bottom x its ground in the cell
        bottom_elev = self.top[self.map.node] - self.depth[self.map.obj]
        row = self.sy[self.map.node] * np.clip(h0[self.map.node] - bottom_elev, 0.0, None) * self.map.area
        if self.cap == "room":
            row = np.minimum(row, self._room_rows(h0))
        self.row_planned = row
        self.planned = np.bincount(self.map.node, weights=row, minlength=mf.nodes)[self.nodes]
        mf.set_entries("Q", self.package, -self.planned)
        return None

    def after_solve(self) -> tuple[float, float]:
        """What left the cells on the day just solved (m3), split for SWAT+'s next day; returns it and the part that
        went to bottom layers because the water table was below that HRU's profile."""
        out = np.clip(-self.mf6.simvals(self.package), 0.0, None)  # m3 per entry, out of the aquifer
        cell = np.zeros(self.mf6.nodes)
        cell[self.nodes] = out
        if self.rule in ("drain", "water_table"):
            row = cell[self.map.node] * self.shares
            head = self.mf6.head()
        else:  # the well's rate as taken (MODFLOW 6 may reduce it in a drying cell), by each row's plan
            planned = np.zeros(self.mf6.nodes)
            planned[self.nodes] = self.planned
            frac = np.divide(cell, planned, out=np.zeros_like(cell), where=planned > 0)
            row = self.row_planned * frac[self.map.node]
            head = self.h0
        self.rows, self.row_depth = row, (self.top - head)[self.map.node]  # for the mass with this water
        if self.rule == "water_table":
            # the top is the profile's bottom: the depth to water below the land surface is the HRU's profile less
            # the head's height above the top; the held layers for tomorrow, from today's end-of-day heads
            self.row_depth = self.depth[self.map.obj] - (head - self.bottom_elev)[self.map.node]
            height = hru_table_height(self.map, head, self.bottom_elev, len(self.depth))
            self.flags = held_layers(height, self.depth, self.layer_hru, self.top_depth, self.bottom, self.flags)
        vol, at_bottom = to_layers(row, self.row_depth, self.map, self.layer_hru, self.bottom)
        self.due = layer_rates(vol, self.layer_hru, self.hru_area, DAY)
        return float(out.sum()), at_bottom

    def mass_split(self, x: np.ndarray) -> np.ndarray:
        """The mass (kg) the day's water just solved carried out of each cell at its concentration ``x`` (kg m-3 per
        node, the transport model's at the end of the day, as its outflow carries it), split among the HRUs' layers as
        the water was (``to_layers``, at the same depths): kg per layer, for SWAT+'s next day."""
        mass, _ = to_layers(self.rows * x[self.map.node], self.row_depth, self.map, self.layer_hru, self.bottom)
        return mass
