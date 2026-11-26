"""Shared set-up for the land route's tests: Osu_1hru made upland (its HRU is a rice paddy as shipped, and SWAT+ puts
a paddy's irrigation into its ponded water), and irrigators reported day by day."""

from __future__ import annotations

import datetime as dt
import shutil

from conftest import ROOT, RUNS

from basinledger import config as C
from basinledger.config import LandConfig, UserConfig
from basinledger.coupling import swatplus as S

DAY = 86400.0
WELL = (((0, 5, 6), 0.75), ((0, 7, 6), 0.25))


def upland(src, dest):
    """A copy of a SWAT+ project with Osu_1hru's paddy taken off its HRU, on the agrl land use, as swatplus-bmi's
    own tests make it."""
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(src, dest)
    p = dest / "hru-data.hru"
    p.write_text(p.read_text().replace("rice140_lum", "agrl_lum").replace("paddy0001", "null"))
    return dest


def irrigator(tmp, name, start, days, m3_s, consumed, withdrawal, land):
    w = m3_s * DAY
    rows = ["date,withdrawal_m3,consumption_m3,return_m3"]
    rows += [f"{start + dt.timedelta(days=d)},{w},{consumed * w},{(1 - consumed) * w}" for d in range(days)]
    (tmp / f"{name}.csv").write_text("\n".join(rows) + "\n")
    (tmp / f"{name}.toml").write_text(f'series = "{name}.csv"\nstart = "{start}"\n')
    return UserConfig(name=name, component="basinledger.generic:ReportedUser", config=tmp / f"{name}.toml",
                      withdrawal=withdrawal, returns=f"land:{land}")


def land_config(example, tmp, days, golf_runoff=0.2, stated=None):
    """The example with two irrigators on Osu_1hru's one HRU: a golf course on an intake on the channel, 30% of its
    water lost in the application and a share running off, and a farm on a well over two cells."""
    cfg = C.load(ROOT / "examples" / "synthetic" / example)
    cfg.run_folder = RUNS
    cfg.swat_project = upland(cfg.swat_project, tmp / "Osu_1hru_upland")
    cfg.places_of_use = list(cfg.places_of_use) + [("golf_intake", "cha01")]
    cfg.wells = list(cfg.wells) + [C.WellConfig("farm", "wel_farm", None, True, cells=WELL)]
    cfg.land = [LandConfig("golf", (("hru0001", 1.0),), golf_runoff, dict(stated or {})),
                LandConfig("farm", (("hru0001", 1.0),), 0.0, {})]
    start = S._period(cfg.swat_project)[0]
    cfg.users = list(cfg.users) + [
        irrigator(tmp, "golf", start, days, 0.0004, 0.3, "swat:place_of_use:golf_intake", "golf"),
        irrigator(tmp, "farm", start, days, 0.0003, 0.1, "well:farm", "farm")]
    for settings in cfg.salt.values():  # sources by land use follow the HRU onto agrl
        for src in settings.sources:
            if src.by == "landuse" and isinstance(src.rate, dict):
                src.rate = {("agrl_lum" if k == "rice140_lum" else k): v for k, v in src.rate.items()}
            if src.group == "septic":
                src.outside_routes = True
    cfg.outputs = {"channels": ["cha01"]}
    return cfg
