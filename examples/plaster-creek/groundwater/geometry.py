"""Hydrofabric geometry for the groundwater model: catchments and flowpaths from the NGIAB run's GeoPackage, and the
maps between them and a MODFLOW 6 grid.

The hydrofabric v2.2 subset is in EPSG:5070 (CONUS Albers, meters), and the model grid is built in the same
projection, so nothing is reprojected here. Maps are written as the CSVs Basin Ledger reads:

- catchment_cells.csv: object (cat- id), layer, row, col, area_m2: the overlap of each catchment with each cell;
- river_cells.csv: object (wb- id), layer, row, col, length_m: the length of each flowpath in each cell.

On a vertex grid (DISV) the helpers take GridIntersect's icell2d instead, and the maps' cell columns are layer,
icell2d.
"""

from __future__ import annotations

import csv
import sqlite3
import struct
from pathlib import Path

import numpy as np
import shapely
from shapely import wkb


def _gpkg_geometry(blob: bytes):
    """A GeoPackage geometry blob (header, then WKB) as a shapely geometry."""
    if blob[:2] != b"GP":
        raise ValueError("not a GeoPackage geometry")
    flags = blob[3]
    envelope = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}[(flags >> 1) & 0b111]
    return wkb.loads(blob[8 + envelope:])


def read_layer(gpkg: Path, table: str, id_col: str) -> dict[str, object]:
    with sqlite3.connect(gpkg) as db:
        srs = db.execute("select srs_id from gpkg_contents where table_name = ?", (table,)).fetchone()[0]
        if srs != 5070:
            raise ValueError(f"{table} is in EPSG:{srs}; expected 5070")
        return {i: _gpkg_geometry(g) for i, g in db.execute(f'select "{id_col}", geom from "{table}"')}


def catchments(gpkg: Path) -> dict[str, object]:
    return read_layer(gpkg, "divides", "divide_id")


def flowpaths(gpkg: Path) -> dict[str, object]:
    return read_layer(gpkg, "flowpaths", "id")


def grid_intersect(modelgrid):
    from flopy.utils.gridintersect import GridIntersect
    return GridIntersect(modelgrid)


def _cell(cellid) -> tuple[int, ...]:
    """GridIntersect's cell id: (row, col) on a structured grid, icell2d on a vertex grid."""
    return tuple(int(c) for c in np.atleast_1d(cellid))


def catchment_cells(gi, polys: dict[str, object]) -> list[tuple]:
    """(catchment, row, col, area m2) for every overlap of a catchment with a top-layer cell; (catchment, icell2d,
    area m2) on a vertex grid."""
    rows = []
    for cid, poly in polys.items():
        r = gi.intersect(poly, geo_dataframe=False)
        for cell, area in zip(r["cellids"], r["areas"]):
            if area > 0:
                rows.append((cid, *_cell(cell), float(area)))
    return rows


def flowpath_cells(gi, lines: dict[str, object]) -> list[tuple]:
    """(flowpath, row, col, length m) for every cell a flowpath crosses; (flowpath, icell2d, length m) on a vertex
    grid."""
    rows = []
    for fid, line in lines.items():
        r = gi.intersect(line, geo_dataframe=False)
        for cell, length in zip(r["cellids"], r["lengths"]):
            if length > 0:
                rows.append((fid, *_cell(cell), float(length)))
    return rows


def write_map(path: Path, rows, value: str, layer: int = 0):
    """A Basin Ledger map (mapping.py): object, layer, row, col, value, or object, layer, icell2d, value for rows
    from a vertex grid."""
    rows = list(rows)
    cell = ["row", "col"] if not rows or len(rows[0]) == 4 else ["icell2d"]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["object", "layer", *cell, value])
        for obj, *c, v in rows:
            w.writerow([obj, layer, *c, repr(float(v))])


def basin_outline(polys: dict[str, object]):
    return shapely.union_all(list(polys.values()))


def extent(outline, buffer_m: float, cell_m: float) -> tuple[float, float, int, int]:
    """Grid origin (lower left x, y) and size (rows, cols) covering the outline plus a buffer, snapped to cell_m."""
    x0, y0, x1, y1 = outline.buffer(buffer_m).bounds
    x0, y0 = np.floor(x0 / cell_m) * cell_m, np.floor(y0 / cell_m) * cell_m
    x1, y1 = np.ceil(x1 / cell_m) * cell_m, np.ceil(y1 / cell_m) * cell_m
    return float(x0), float(y0), int(round((y1 - y0) / cell_m)), int(round((x1 - x0) / cell_m))
