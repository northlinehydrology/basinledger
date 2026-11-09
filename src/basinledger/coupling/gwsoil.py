"""Groundwater up into the soil: MODFLOW 6's water to SWAT+'s HRU soil layers (docs/design.md, Groundwater up into the
soil). Built Nov 8-9, 2026 on swatplus-bmi 0.7.0's ``soil_layer_water__groundwater_inflow_rate``. Two rules for how
much water leaves each cell (``rule = "drain"``, chosen Nov 9, or ``"gwflow"``) and an optional cap (``cap = "none"``,
chosen, or ``"room"``), both required in a config. The pieces:

- ``cell_shares``: each HRU's share of a cell's mapped ground, which splits a cell's water among the HRUs over it;
- ``to_layers``: water from each (cell, HRU) row into that HRU's layers below the cell's water table, each in
  proportion to its thickness below it (gwflow's split);
- ``drain_elevation``: the bottom of the soil profile under each cell, MODFLOW 6's top less the profile depths of the
  HRUs over it, by their shares, for a drain at the profile's bottom.

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
    """Volumes (m3) per map row into the HRU's soil layers, each row split at its own cell's depth to water.
    Returns m3 per layer, and the volume of rows whose water table is below their HRU's profile: that water still
    left the cell (a drain under a cell is at the deepest profile over it, by share), so it goes to the HRU's bottom
    layer, where it reaches the soil."""
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


# ---------- the link ----------

INFLOW = "soil_layer_water__groundwater_inflow_rate"  # per HRU soil layer, m s-1, set for the next day (0.7.0)
HRU_INFLOW = "soil_water__groundwater_inflow_rate"    # per HRU, m s-1, gwsoilq
BOTTOM = "soil_layer_bottom__depth"                   # m below the surface
WATER = "soil_layer_water__depth"                     # m above the wilting point (st)
SATURATION = "soil_layer_water__saturation_depth"     # m above the wilting point (ul)
RULES = ("drain", "gwflow")
CAPS = ("none", "room")
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
                 due: np.ndarray | None = None):
        """SWAT+'s side, before either model steps; ``attach`` adds MODFLOW 6's once its package has its entries
        (after the first ``prepare_day``). ``amap``: each HRU's ground in each cell, before ``conserving``.
        ``due``: what a saved run had due for SWAT+'s next day."""
        if rule not in RULES:
            raise ValueError(f"[groundwater_to_soil] rule must be one of {RULES}")
        if cap not in CAPS:
            raise ValueError(f"[groundwater_to_soil] cap must be one of {CAPS}")
        if not sw.has(INFLOW):
            raise ValueError(f"[groundwater_to_soil] needs swatplus-bmi's {INFLOW} (0.7.0)")
        self.sw, self.map, self.rule, self.cap, self.package = sw, amap, rule, cap, package
        self.mf6 = None
        self.hru_area = np.asarray(hru_area, float)
        hrus = list(amap.objects)
        names = sw.names(INFLOW)
        self.layer_hru = np.array([hrus.index(n.split(":")[0]) for n in names])
        self.bottom = sw.get(BOTTOM).astype(float)
        self.depth = profile_depth(self.layer_hru, self.bottom, len(hrus))
        self.due = np.zeros(len(names)) if due is None else np.asarray(due, float)  # m s-1 per layer, SWAT+'s next day
        self.raised = 0

    def attach(self, mf6):
        n = mf6.nodes
        self.mf6 = mf6
        self.shares = cell_shares(self.map, n)
        self.top = mf6.get("TOP", mf6.model, "DIS").astype(float)[:n]
        nodes = mf6.package_nodes(self.package)
        mapped = np.unique(self.map.node)
        if sorted(nodes.tolist()) != mapped.tolist():
            raise ValueError(f"MODFLOW 6's {self.package} package must have one entry on each mapped cell and no "
                             f"other ({len(nodes)} entries, {len(mapped)} mapped cells)")
        self.nodes = nodes
        self.elev0 = drain_elevation(self.map, self.top, self.depth, n)[nodes]
        if self.rule == "gwflow":
            self.sy = mf6.get("SY", mf6.model, "STO").astype(float)[:n]
        self.planned = np.zeros(len(nodes))

    # -- SWAT+'s side --
    def due_volume(self) -> float:
        return float(np.sum(self.due * DAY * self.hru_area[self.layer_hru]))

    def before_swat(self) -> float:
        """Set what's due in SWAT+'s layers for its coming day; returns its volume, m3."""
        self.sw.set(INFLOW, self.due)
        return self.due_volume()

    def received(self) -> float:
        """What SWAT+ reports went into its layers on the day just done, m3."""
        return float(np.sum(self.sw.get(INFLOW) * DAY * self.hru_area[self.layer_hru]))

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
        mf, P = self.mf6, self.package.upper()
        h0 = mf.head()
        self.h0 = h0
        if self.rule == "drain":
            elev = self.elev0.copy()
            self.raised = 0
            if self.cap == "room":
                # the drain raised so it gives no more than the room at the day's starting heads: q = C (h - elev)
                cond = mf.get("COND", mf.model, P).astype(float)[: len(self.nodes)]
                cap = self._cell_cap(h0)
                want = h0[self.nodes] - cap / np.where(cond > 0, cond, np.inf)
                elev = np.maximum(self.elev0, want)
                self.raised = int(np.sum(elev > self.elev0))
            mf.set("ELEV", self.package, elev)
            return
        # gwflow: Sy x the water table's height above each HRU's profile bottom x its ground in the cell
        bottom_elev = self.top[self.map.node] - self.depth[self.map.obj]
        row = self.sy[self.map.node] * np.clip(h0[self.map.node] - bottom_elev, 0.0, None) * self.map.area
        if self.cap == "room":
            row = np.minimum(row, self._room_rows(h0))
        self.row_planned = row
        self.planned = np.bincount(self.map.node, weights=row, minlength=mf.nodes)[self.nodes]
        mf.set("Q", self.package, -self.planned)
        return None

    def after_solve(self) -> tuple[float, float]:
        """What left the cells on the day just solved (m3), split for SWAT+'s next day; returns it and the part that
        went to bottom layers because the water table was below that HRU's profile."""
        out = np.clip(-self.mf6.simvals(self.package), 0.0, None)  # m3 per entry, out of the aquifer
        cell = np.zeros(self.mf6.nodes)
        cell[self.nodes] = out
        if self.rule == "drain":
            row = cell[self.map.node] * self.shares
            head = self.mf6.head()
        else:  # the well's rate as taken (MODFLOW 6 may reduce it in a drying cell), by each row's plan
            planned = np.zeros(self.mf6.nodes)
            planned[self.nodes] = self.planned
            frac = np.divide(cell, planned, out=np.zeros_like(cell), where=planned > 0)
            row = self.row_planned * frac[self.map.node]
            head = self.h0
        vol, at_bottom = to_layers(row, (self.top - head)[self.map.node], self.map, self.layer_hru, self.bottom)
        self.due = layer_rates(vol, self.layer_hru, self.hru_area, DAY)
        return float(out.sum()), at_bottom
