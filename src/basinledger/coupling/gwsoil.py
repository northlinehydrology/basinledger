"""Groundwater up into the soil: MODFLOW 6's water to SWAT+'s HRU soil layers (docs/design.md, Groundwater up into the
soil). Designed Nov 8, 2026; the parts here don't depend on the rule still to be chosen (how much water leaves each
cell) or on swatplus-bmi 0.7.0 (``soil_layer_water__groundwater_inflow_rate``):

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
    Returns m3 per layer. A row whose water table is below its HRU's profile puts nothing in; its volume is
    returned as the second value, so a caller can refuse it or leave it in MODFLOW 6."""
    layer_hru = np.asarray(layer_hru)
    layer_bottom = np.asarray(layer_bottom, float)
    order = {h: np.flatnonzero(layer_hru == h) for h in np.unique(layer_hru)}
    tops = np.empty_like(layer_bottom)
    for idx in order.values():
        b = layer_bottom[idx]
        tops[idx] = np.concatenate([[0.0], b[:-1]])
    out = np.zeros(len(layer_bottom))
    unplaced = 0.0
    for v, wd, h in zip(row_volume, row_water_depth, amap.obj):
        if v == 0.0:
            continue
        idx = order.get(h)
        s = layer_split(tops[idx], layer_bottom[idx], wd) if idx is not None else np.zeros(0)
        if s.sum() == 0:
            unplaced += float(v)
            continue
        out[idx] += v * s
    return out, unplaced


def layer_rates(layer_volume: np.ndarray, layer_hru: np.ndarray, hru_area: np.ndarray, day_s: float = 86400.0):
    """m3 per layer over a day -> m s-1 over each layer's HRU area, as swatplus-bmi takes it."""
    return layer_volume / hru_area[np.asarray(layer_hru)] / day_s
