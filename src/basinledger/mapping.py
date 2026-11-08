"""Maps between SWAT+ objects and MODFLOW 6 cells, as areas of overlap.

A map is a list of rows (object, cell, area in m2): the part of a cell's area that a SWAT+ object (an aquifer, an HRU,
a catchment) covers. Depth rates move through it by area, so volume is kept exactly: an object's rate times its
mapped area is what the cells receive. Maps come from an example's config, never from code:

    cells = "all"             every cell of the top layer, its whole area, to one object
    cells = "map.csv"         a table: object, layer, row, col, area_m2 (or object, node, area_m2)

A map from GIS (NextGen catchments or SWAT+ HRU polygons intersected with the grid, e.g. with flopy's GridIntersect)
is written to the same table, so the coupling code never sees geometry. ``pixel_rows`` writes it from a raster (HRUs
built from land cover, soil and slope grids, as QSWAT+ and scripted builds do), splitting each pixel exactly among
the cells it overlaps.

HRU by HRU (``[recharge] source = "hru"``), the table gives each HRU's ground in each cell, and ``conserving`` scales
each HRU's rows so they add to its area in SWAT+: the cells then take each HRU's share of its percolation by the share
of its ground in them, and the volume handed to MODFLOW 6 is SWAT+'s percolation times its HRU areas, exactly.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class AreaMap:
    objects: list[str]   # object names, in the order of the source model's grid
    obj: np.ndarray      # row -> index into objects
    node: np.ndarray     # row -> MODFLOW 6 node, 0-based, reduced numbering
    area: np.ndarray     # row -> m2 of overlap

    def mapped_area(self) -> np.ndarray:
        """m2 per object."""
        return np.bincount(self.obj, weights=self.area, minlength=len(self.objects))

    def to_cells(self, rate: np.ndarray, cell_area: np.ndarray) -> np.ndarray:
        """Depth rates per object -> depth rates per cell (same units), conserving volume."""
        vol = np.bincount(self.node, weights=rate[self.obj] * self.area, minlength=len(cell_area))
        return vol / cell_area

    def to_objects(self, values: np.ndarray) -> np.ndarray:
        """Values per cell (heads, drawdown) -> the overlap-area-weighted mean per object: the map run in reverse."""
        w = np.bincount(self.obj, weights=self.area, minlength=len(self.objects))
        return np.bincount(self.obj, weights=values[self.node] * self.area, minlength=len(self.objects)) / w

    def outside(self, n_cells: int) -> np.ndarray:
        """Cells no object covers (reported together as outside the basin)."""
        return np.setdiff1d(np.arange(n_cells), self.node)

    def volume(self, rate: np.ndarray) -> float:
        """Depth rate per object times its mapped area, summed: what leaves the source model."""
        return float(np.sum(rate * self.mapped_area()))

    def conserving(self, source_area: np.ndarray) -> "AreaMap":
        """The map with each object's rows scaled to add to its area in the source model (``source_area``, m2), so
        a depth rate times the object's own area is what its cells receive: placement by the map's shares, volume
        from the source model. An object the source model has but the map doesn't place is refused."""
        source_area = np.asarray(source_area, float)
        m = self.mapped_area()
        missing = (m <= 0) & (source_area > 0)
        if missing.any():
            names = [self.objects[i] for i in np.flatnonzero(missing)]
            raise ValueError(f"{len(names)} objects have area in the source model but no cells in the map: "
                             f"{names[:5]}")
        scale = np.divide(source_area, m, out=np.zeros_like(m), where=m > 0)
        return AreaMap(self.objects, self.obj, self.node, self.area * scale[self.obj])

    def check(self, cell_area: np.ndarray, source_area: np.ndarray | None = None, rtol: float | None = None):
        """Cells may not be over-covered; with ``rtol``, each object's mapped area must equal its area in the
        source model (``source_area``, m2) within that tolerance."""
        covered = np.bincount(self.node, weights=self.area, minlength=len(cell_area))
        if np.any(covered > cell_area * (1 + 1e-9)):
            bad = int(np.argmax(covered / cell_area))
            raise ValueError(f"cell {bad} is mapped {covered[bad]:.1f} m2 of {cell_area[bad]:.1f} m2")
        if rtol is not None and source_area is not None:
            m = self.mapped_area()
            off = np.abs(m - source_area) > rtol * source_area
            if off.any():
                i = int(np.argmax(off))
                raise ValueError(f"{self.objects[i]}: mapped area {m[i]:.0f} m2, source model {source_area[i]:.0f} m2")

    @classmethod
    def all_cells(cls, objects: list[str], name: str, nodes: np.ndarray, cell_area: np.ndarray) -> "AreaMap":
        nodes = np.asarray(nodes)
        return cls(list(objects), np.full(len(nodes), objects.index(name)), nodes, cell_area[nodes].astype(float))

    @classmethod
    def from_csv(cls, path: str | Path, objects: list[str], node_of) -> "AreaMap":
        """``node_of(layer, row, col)`` -> 0-based reduced node; rows may give ``node`` instead."""
        obj, node, area = [], [], []
        with open(path, newline="") as f:
            for r in csv.DictReader(f):
                if r["object"] not in objects:
                    raise KeyError(f"{path}: {r['object']!r} is not an object of the source model")
                obj.append(objects.index(r["object"]))
                if r.get("node"):
                    node.append(int(r["node"]))
                else:
                    node.append(node_of(int(r["layer"]), int(r["row"]), int(r["col"])))
                area.append(float(r["area_m2"]))
        return cls(list(objects), np.array(obj), np.array(node), np.array(area))


def pixel_rows(obj: np.ndarray, weight: np.ndarray, x0: np.ndarray, y0: np.ndarray, size: float,
               col_edges: np.ndarray, row_edges: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Square pixels split exactly among the cells of a structured grid with no rotation, by area of overlap.

    ``obj`` and ``weight`` give each pixel's object and the share of the pixel it holds (1 for a pixel wholly one
    object; a pixel folded into several objects, as small classes are when HRUs are defined, appears once per object
    with weights adding to 1). ``x0``, ``y0`` are each pixel's lower-left corner and ``size`` its side, in the grid's
    coordinates. ``col_edges`` rise west to east and ``row_edges`` fall north to south (row 0 at the top), as flopy's
    ``xvertices``/``yvertices`` give them. Returns (object, row, col, area m2), summed by object and cell. Pixel parts
    outside the grid are left out (``area`` then adds to less than the pixels' area, which the map's area check sees).
    A pixel may be no larger than a cell, so it touches at most two columns and two rows."""
    obj, weight = np.asarray(obj), np.asarray(weight, float)
    x0, y0 = np.asarray(x0, float), np.asarray(y0, float)
    col_edges, row_edges = np.asarray(col_edges, float), np.asarray(row_edges, float)
    if size > min(np.diff(col_edges).min(), -np.diff(row_edges).max()) + 1e-9:
        raise ValueError("pixels larger than a cell: map them as polygons instead")
    out_o, out_r, out_c, out_a = [], [], [], []
    ys = -row_edges  # rising, so searchsorted works north to south as rows do
    for dx in (0, 1):
        c = np.searchsorted(col_edges, x0, side="right") - 1 + dx
        okc = (c >= 0) & (c < len(col_edges) - 1)
        cc = np.clip(c, 0, len(col_edges) - 2)
        wx = np.clip(np.minimum(x0 + size, col_edges[cc + 1]) - np.maximum(x0, col_edges[cc]), 0, None) * okc
        for dy in (0, 1):
            r = np.searchsorted(ys, -(y0 + size), side="right") - 1 + dy
            okr = (r >= 0) & (r < len(row_edges) - 1)
            rr = np.clip(r, 0, len(row_edges) - 2)
            wy = np.clip(np.minimum(y0 + size, row_edges[rr]) - np.maximum(y0, row_edges[rr + 1]), 0, None) * okr
            a = wx * wy * weight
            k = a > 0
            out_o.append(obj[k]); out_r.append(rr[k]); out_c.append(cc[k]); out_a.append(a[k])
    o, r, c, a = (np.concatenate(v) for v in (out_o, out_r, out_c, out_a))
    key = np.stack([o, r, c], axis=1)
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    area = np.bincount(inv.ravel(), weights=a, minlength=len(uniq))
    return uniq[:, 0], uniq[:, 1], uniq[:, 2], area
