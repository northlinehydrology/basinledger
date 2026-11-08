"""SWAT+ recharge into MODFLOW 6, one day at a time (one-way).

Two sources, chosen by ``[recharge] source``:

- **aquifer** (the default): each SWAT+ aquifer's recharge as a depth rate (``soil_water_table__recharge_rate``, m
  s-1), spread over the cells its map gives it. Every cell under an aquifer takes the same rate.
- **hru**: each HRU's percolation out of the bottom of its soil profile (``soil_water__percolation_rate``, m s-1,
  swatplus-bmi 0.6.0), spread over the cells under that HRU's ground, so recharge follows land cover, soil and slope
  within a catchment. The map is scaled to each HRU's area in SWAT+ (``AreaMap.conserving``), so the volume sent is
  SWAT+'s percolation times its HRU areas, exactly, and nothing is booked ``by_the_map``. SWAT+ 62.0.1 puts the same
  day's percolation into the HRU's aquifer (``hru_hyds.f90``'s recharge hydrograph, summed by ``aqu_1d_control.f90``
  with its delay commented out), so where every HRU recharges an aquifer the total equals the aquifer path's, day by
  day: the driver books SWAT+'s own aquifer recharge beside it (``recharge:swat_aquifers``) and holds them equal.

Both are read through swatplus-bmi after the day, and the cells' rates go into the recharge package for the same day. ``sent`` is the volume leaving SWAT+ (rate x mapped area); after the solve,
``received`` is what MODFLOW 6's recharge package took in. The ledger holds them equal.

``made`` is what the surface model made: rate x each object's own area in it. Where a map's area differs from the
source's, sent differs from made, and the difference is booked as ``recharge:by_the_map``: water the map added
(positive) or couldn't place (negative). The area check (``[recharge] area_rtol``, 1e-6 unless the example turns it
off with ``area_check = false`` and says why) holds it to a tolerance.
"""

from __future__ import annotations

import numpy as np

from ..mapping import AreaMap
from ..units import DAY
from .modflow6 import Modflow6

RECHARGE = "soil_water_table__recharge_rate"
PERCOLATION = "soil_water__percolation_rate"  # per HRU, swatplus-bmi 0.6.0
SOURCES = {"aquifer": RECHARGE, "hru": PERCOLATION}
# HRU path: the percolation sent against SWAT+'s own aquifer recharge, each read as float32 and each area rounded to
# 1e-5 ha in the project's files [A]
AQUIFER_RTOL = 1e-5


class RechargeLink:
    def __init__(self, swat, mf6: Modflow6, amap: AreaMap, package: str = "rch", source_area: np.ndarray | None = None,
                 var: str = RECHARGE):
        self.swat, self.mf6, self.map, self.package, self.var = swat, mf6, amap, package, var
        self.source_area = amap.mapped_area() if source_area is None else np.asarray(source_area, float)
        self.made = 0.0
        grid = swat.get_var_grid(var)
        self._buf = np.empty(swat.get_grid_size(grid), dtype=np.float32)
        names = swat.get_grid_object_names(grid)
        if list(names) != amap.objects:
            raise ValueError(f"area map objects {amap.objects} differ from SWAT+'s {list(names)}")
        if swat.get_var_units(var) != "m s-1":
            raise ValueError(f"unexpected recharge units {swat.get_var_units(var)}")

    def transfer(self) -> float:
        """Read the day's recharge from SWAT+ (after its update) and set it in MODFLOW 6 (after prepare_day).
        Returns the volume sent, m3."""
        self.swat.get_value(self.var, self._buf)
        rate = self._buf.astype(np.float64) * DAY  # m d-1 per aquifer
        if not getattr(self, "_checked", False):
            self.check_package()
            self._checked = True
        self.mf6.set_recharge(self.package, self.map.to_cells(rate, self.mf6.cell_area))
        self.made = float(np.sum(rate * self.source_area))  # what the surface model made, one day
        return self.map.volume(rate)  # one day

    def check_package(self):
        """The recharge package Basin Ledger sets must give each cell one entry with no auxiliary multiplier: an
        entry's rate is its cell's rate, so a second entry or a multiplier would change the volume sent."""
        nodes = self.mf6.package_nodes(self.package)
        unmapped = np.setdiff1d(np.unique(self.map.node), nodes)
        nodes = nodes[nodes >= 0]  # array recharge (RCHA) lists inactive columns as node 0
        if len(unmapped):
            raise ValueError(f"{len(unmapped)} mapped cells have no entry in MODFLOW 6's {self.package} package")
        if len(np.unique(nodes)) != len(nodes):
            raise ValueError(f"MODFLOW 6's {self.package} package has {len(nodes) - len(np.unique(nodes))} cells with "
                             "more than one entry; Basin Ledger sets one rate per cell")
        if int(self.mf6.get("IAUXMULTCOL", self.mf6.model, self.package.upper())[0]) > 0:
            raise ValueError(f"MODFLOW 6's {self.package} package has an auxiliary multiplier; Basin Ledger sets "
                             "the rate itself")

    def received(self) -> float:
        """Volume the recharge package took in on the day just solved, m3."""
        return float(self.mf6.simvals(self.package).sum())
