"""An example's configuration: one TOML file that wires real models for one basin and its water users.

Paths are relative to the config file. ``{swatplus_refdata}`` stands for the example projects that ship with SWAT+
(``refdata/`` in the source that ``swatplus-bmi build`` fetched), for the synthetic example.

    name = "synthetic"
    description = "..."

    [run]
    folder = "../../runs/synthetic"     # each run copies the models here
    days = 0                            # 0: the SWAT+ project's whole period

    [swatplus]                          # the surface model: SWAT+, stepped live through swatplus-bmi
    project = "{swatplus_refdata}/Osu_1hru"
    places_of_use = [{name = "intake", channel = "cha01"}]    # added to the run's copy, set daily
    point_sources = [{name = "wwtp", channel = "cha01"}]

    # or [nextgen]                      # the surface model: a finished NGIAB run, replayed (coupling/nextgen.py)
    # run = "{nextgen_michigan}/runs/scenarios/baseline"

    [modflow6]                          # optional: the groundwater model
    builder = "model:build"             # model.py beside this file: build(workspace, days, start_date)
    builder_args = {}                   # optional keyword arguments for the builder (a scenario's variant)
    # workspace = "mf6"                 # or an existing simulation folder, copied into the run
    # length_units = "feet"             # the model's, only where its DIS leaves LENGTH_UNITS undefined
    # time_units = "days"               # the model's, only where its TDIS leaves TIME_UNITS undefined
    model = "gwf"
    # storage = {sy = 0.15, ss = 1e-5}  # for a model handed over steady: Basin Ledger adds a storage package to the
                                        #   run's copy (or makes the model's own transient, with storage = {}); sy a
                                        #   fraction, ss per meter, the study's values [A] (existing.add_storage)
    # models = ["parent", "inset"]      # or several groundwater models of the simulation coupled together (an LGR
                                        #   inset and its parent): maps then give each row's model in a model
                                        #   column, a well its model, and the exchange, a lake and head outputs
                                        #   theirs (each the first if not given)

    [recharge]                          # SWAT+ -> MODFLOW 6 recharge (one-way)
    package = "rch"
    source = "aquifer"                  # each aquifer's recharge (default), or "hru": each HRU's percolation, placed
                                        #   on the cells under its ground (swatplus-bmi 0.6.0; coupling/recharge.py)
    map = "all:aqu001"                  # every top-layer cell to one aquifer (or HRU), or a CSV (see mapping.py)
    # add = true                        # a model without a recharge package: Basin Ledger adds one, named by package,
                                        #   on the map's cells (a CSV of cells), to the run's copy (existing.py)
    area_rtol = 1e-6                    # mapped area must match each object's area in the surface model (default 1e-6);
    area_check = false                  #   or no check, where an example maps more or less area on purpose (say why)

    [weather]                           # optional, SWAT+ only: series written into its gage records (weather.py)
    precipitation = "weather/pcp.csv"   # a date column, then one column per gage file, mm/day (temperature_max,
    change = {temperature = 2.0}        #   temperature_min, solar, humidity, wind); a scenario's change, on the
                                        #   series or the project's own records
    [groundwater_to_soil]               # optional: MODFLOW 6's water up into SWAT+'s HRU soil layers (coupling/gwsoil.py;
    rule = "drain"                      #   swatplus-bmi 0.7.0). Required, no default: "drain" (a drain at the profile's
    cap = "none"                        #   bottom, solved by MODFLOW 6) or "gwflow" (SWAT+ gwflow's rule, a well);
    package = "soil"                    #   cap "none" or "room" (no more than the layers' room). The package (DRN for
                                        #   drain, WEL for gwflow) has one entry on each mapped cell. Needs recharge by
                                        #   HRU and the two-way [exchange]

    [[wells]]
    name = "pw1"
    package = "wel"
    cell = [0, 20, 30]                  # layer, row, col (0-based) on a structured grid (DIS); layer, icell2d on a
                                        #   vertex grid (DISV), e.g. [0, 412]
    # add = true                        # a well the model doesn't have (a new user's): Basin Ledger adds a WEL
                                        #   package of its own, named by package, to the run's copy (existing.py)
    # well = 3                          # or, for a multi-aquifer well (package a MAW package), its number there,
                                        #   from 1, in place of cell: its rate is set and what it pumped read

    [irrigation]                        # optional: SWAT+'s irrigation pumped from MODFLOW 6's wells (coupling/
    wells = [{hru = "hru0001", well = "irr1", share = 1.0}]  # irrigation.py; swatplus-bmi 0.11.0): each day, the
                                        #   water irrigation put on each HRU, from the [[wells]] named, in their
                                        #   shares; or a CSV with columns hru, well, share. The HRU must irrigate
                                        #   from an outside source in SWAT+

    [[lakes]]                           # optional: a SWAT+ reservoir holding a MODFLOW 6 lake's water (coupling/
    reservoir = "res01"                 #   lakes.py; swatplus-bmi 0.11.0): the lake set from the reservoir each day,
    package = "lak"                     #   its leakage the reservoir's groundwater exchange the next day. lake: its
    lake = 1                            #   1-based number in the package; form "active" (A, within the day from the
    form = "active"                     #   reservoir's start) or "fixed_stage" (B, at the stage of its volume)

    [exchange]                          # optional: two-way stream-aquifer exchange (MODFLOW 6 RIV or SFR <-> SWAT+)
    package = "riv"                     #   a river (RIV) or streamflow routing (SFR) package
    map = "all:cha01"                   # every RIV entry to one channel, or a CSV: object, layer, row, col (or
                                        #   object, layer, icell2d on a DISV grid, or object, node; mapping.py; or
                                        #   object, reach for SFR's 0-based reach numbers)
    limit = true                        # RIV: a losing reach takes no more than its channel carries (default; SFR
                                        #   does this itself)
    limit_tolerance = 0.001             # m: the limit's stages agree with the heads once none moves more than this
                                        #   in an outer iteration (default 1 mm; a large model's stages can chase
                                        #   its heads for dozens of iterations at 1 mm: 1 cm took a 1,157-channel
                                        #   basin from 43 outer iterations a day to 11)
    limit_warm_start = true             # each day's limit starts from yesterday's lowered stages (default); the
                                        #   answer is the same, found in fewer outer iterations
    sfr = "route"                       # SFR: "route" (each channel's reaches route its water, A; a channel that
                                        #   doesn't line up is refused), "route_or_stage" (A where a channel lines up,
                                        #   a fixed stage from SWAT+'s depth where it doesn't, B) or "stage" (B for all)
    drains = {package = "drn", map = "maps/drains.csv"}   # optional: seepage cells, to their channels the next day
    lakes = {package = "lakes", map = "maps/lakes.csv"}   # optional: lake cells at a fixed stage (a RIV package of
                                        #   their own), each to its outlet channel, one of the rivers' or drains'

    [transport.chloride]                # optional: a MODFLOW 6 transport model per species (coupling/transport.py),
    recharge = 0.02                     #   built by the example's builder (add_transport); kg m-3 of recharge [A]
    stream = 0.03                       #   kg m-3 of water a coupled river loses to the aquifer [A] (with [exchange])
    # model = "gwt_chloride"            #   the GWT model's name, if not gwt_<species>

    [transport.nitrate]                 # a species SWAT+ carries (coupling/surface_mass.py), stage 2:
    surface = "no3"                     #   SWAT+'s species through swatplus-bmi: "no3" or "solp" (default the name)
    recharge = "swatplus"               #   SWAT+'s leached mass, HRU by HRU with [recharge] source = "hru"
    stream = "swatplus"                 #   gains, drains and the soil link carry the cells' concentrations into SWAT+,
                                        #   losses SWAT+'s own (needs add_transport's stream_mass on rivers, drains)

    [salt.chloride]                     # optional: a conservative species through the land, the aquifer's
    precipitation = 0.0003              #   [transport.<species>] (its recharge and stream then come from here) and the
                                        #   streams, each day (salt.py; needs SWAT+, MODFLOW 6 in meters, [recharge]
                                        #   source = "hru" and the two-way [exchange]); kg m-3 in precipitation [A]
    road_salt = {annual = "road_salt_hru.csv", repeat_last = 5, repeat_first = 1}  # kg m-2 of the HRU by calendar
                                        #   year (hru, year, kg_m2) over the year's snow and ice days [A]; or rate and
                                        #   share, or series and share
    sources = [{name = "septic", rate = "septic.csv", timing = "uniform"},          # land sources, kg m-2 a year:
               {name = "potash", rate = "potash.csv", timing = "fertilizer", group = "agricultural"}]  # a number,
                                        #   hru, year, kg_m2, hru, kg_m2_yr or (by = "landuse") landuse, kg_m2_yr;
                                        #   timing uniform, fertilizer or 12 monthly weights [A]
    point_sources = {wwtp = 0.08}       # a point source's concentration, for water no user reports [A]
    multipliers = {road_salt = 1.0, agricultural = 1.0}  # on each input group's mass (calibration: model salt,
                                        #   target <species>:<group>)
    surface_store = true                # a snowpack and road store in front of the soil's [A]
    soil_start = 0.002                  # starting concentrations, kg m-3 [A]
    stream_start = 0.002
    [salt.weathering]                   # a second species: weathering's ions, from GWT's weathering source by zone
    precipitation = 0.0
    [salt.conductance]                  # specific conductance (uS/cm) from the species' mg/L, per study [A]
    base = 20.0
    slopes = {chloride = 3.0, weathering = 1.5}

    [rivers]                            # optional: MODFLOW 6 river cells by flowpath, for stage and base flow out
    package = "riv"
    map = "maps/river_cells.csv"        # object (flowpath or channel), layer, row, col (or layer, icell2d)
    stage = "nextgen_depth"             # or "fixed": the stage in MODFLOW 6's input
    bed_thickness = 1.0                 # m [A]: stage = RBOT + bed thickness + the day's depth

    [[utilities]]                       # optional: a water system (see chain.py)
    name = "city"
    sources = [{route = "well:pw1", share = 1.0}]
    outfall = "swat:point_source:wwtp"

    [[rules]]                           # optional: feedback rules (see rules.py)
    name = "low-flow cut"
    kind = "curtailment"
    user = "plant"
    channel = "cha01"
    threshold = 0.001

    [outputs]                           # optional: daily series written with the ledger
    channels = ["cha01"]                # outflow, m3 s-1 (default: every channel)
    heads = [{name = "pw1", cell = [0, 20, 28]}]
    temperature = true                  # SWAT+ prints what the reach heat balance needs (temperature.py)
    balance = "cha01"                   # the outlet: SWAT+ prints what the whole-basin balance needs (balance.py)

    [calibration]                       # optional: the example's calibration record
    record = "calibration/record.toml"  # without one, the example is a screening run
    parameters = "calibration/parameters.csv"  # optional: applied on every run (calibration.py's table)
    zones = "calibration/zones.csv"     # optional: MODFLOW 6 cells' zones for parameters by zone
    pilot_points = {file = "calibration/pilot_points.csv", range = 2000.0, nugget = 0.0}  # optional: pp=<set>

    [[users]]
    name = "plant"
    component = "basinledger.testing:SyntheticUser"
    config = "user.toml"
    withdrawal = "well:pw1"             # source: a MODFLOW 6 well (later: a SWAT+ place of use, a utility)
    return = "outfall:creek"            # sink: an outfall account (later: a SWAT+ point source, a sewer)
    # variables = {withdrawal = "...", gain = "..."}   for components without the common names
    # forcing = {source = "nextgen", catchment = "cat-652059"}   hourly weather inputs (forcing.py)
    # forcing = {source = "swatplus", station = "sta035"}         or from a SWAT+ station's daily files
    # inputs = {water_user_intake_water_chloride__mass_concentration = 0.012}   constant inputs [A]
    # report = {unremoved_heat_MWh = "data_center__unremoved_heat_power"}      booked as user:<name>:<key> daily

Routes also include ``offsite:<name>``: water from or to somewhere outside the basin (a Great Lakes intake, a plant on
another river), entered in the ledger only.

``{nextgen_michigan}`` stands for a checkout of github.com/northlinehydrology/nextgen-michigan (``$NEXTGEN_MICHIGAN``,
or ``~/repos/nextgen-michigan``).
"""

from __future__ import annotations

import os
import tomllib
from dataclasses import dataclass, field
from pathlib import Path


def cache_dir() -> Path:
    """Where examples keep downloaded and generated data, outside git: $BASINLEDGER_CACHE or ~/.cache/basinledger."""
    return Path(os.environ.get("BASINLEDGER_CACHE", Path.home() / ".cache" / "basinledger"))


def nextgen_michigan() -> Path:
    return Path(os.environ.get("NEXTGEN_MICHIGAN", Path.home() / "repos" / "nextgen-michigan"))


def swatplus_refdata() -> Path:
    if os.environ.get("SWATPLUS_REFDATA"):
        return Path(os.environ["SWATPLUS_REFDATA"])
    from swatplus_bmi.build import DEFAULT_TAG, home
    return home() / "build" / DEFAULT_TAG / "swatplus" / "refdata"


@dataclass
class UserConfig:
    name: str
    component: str
    config: Path | None
    withdrawal: str
    returns: str
    variables: dict = field(default_factory=dict)
    constituents: dict | None = None
    forcing: dict | None = None
    inputs: dict = field(default_factory=dict)   # constant inputs, e.g. intake concentrations [A]
    report: dict = field(default_factory=dict)   # other outputs booked each day (users.WaterUser)


@dataclass
class WellConfig:
    name: str
    package: str
    cell: tuple[int, ...] | None  # (layer, row, col) on a DIS grid, (layer, icell2d) on a DISV grid; None for MAW
    add: bool = False       # a well the model doesn't have, in a WEL package Basin Ledger adds (existing.add_wells)
    model: str | None = None  # with several models coupled, the one it's in (the first if not given)
    well: int | None = None   # a multi-aquifer well (MAW): its number in the package, 1-based, in place of a cell


def _flag(path: Path, where: str, v) -> bool:
    if not isinstance(v, bool):
        raise ValueError(f"{path}: {where} must be true or false; got {v!r}")
    return v


def _pilot_points(path, base, pp) -> dict | None:
    """``[calibration] pilot_points``: {file, range, nugget}, the file resolved."""
    if pp is None:
        return None
    if not isinstance(pp, dict) or "file" not in pp or "range" not in pp or set(pp) - {"file", "range", "nugget"}:
        raise ValueError(f"{path}: [calibration] pilot_points is {{file, range, nugget}} (an exponential variogram's "
                         "range in the grid's length units; nugget optional, 0 to 1)")
    return {"file": _path(base, pp["file"]), "range": float(pp["range"]), "nugget": float(pp.get("nugget", 0.0))}


def _well(path: Path, w: dict) -> WellConfig:
    """A [[wells]] entry: a cell for a WEL package, or ``well``, a multi-aquifer well's number in a MAW package."""
    where = f"well {w.get('name')!r}"
    if "well" in w:
        if "cell" in w or w.get("add"):
            raise ValueError(f"{path}: {where}: a MAW well is named by its number (well), not a cell, and isn't added")
        if not isinstance(w["well"], int) or isinstance(w["well"], bool) or w["well"] < 1:
            raise ValueError(f"{path}: {where}: well is the MAW well's number in its package, from 1; got {w['well']!r}")
        return WellConfig(w["name"], w.get("package", "maw"), None, False, w.get("model"), w["well"])
    return WellConfig(w["name"], w.get("package", "wel"), _cell(path, where, w.get("cell")),
                      _flag(path, f"{where} add", w.get("add", False)), w.get("model"))


def _cell(path: Path, where: str, cell) -> tuple[int, ...]:
    """A cell as configured: 0-based layer, row, col (DIS) or layer, icell2d (DISV). Which one the model's grid takes
    is checked when MODFLOW 6 starts (``Modflow6.node``)."""
    if not isinstance(cell, (list, tuple)) or len(cell) not in (2, 3) or \
            not all(isinstance(c, int) and not isinstance(c, bool) for c in cell):
        raise ValueError(f"{path}: {where} cell must be [layer, row, col] (DIS) or [layer, icell2d] (DISV), "
                         f"0-based integers; got {cell!r}")
    return tuple(cell)


@dataclass
class ExampleConfig:
    path: Path
    name: str
    description: str
    run_folder: Path
    days: int
    swat_project: Path | None
    mf6_builder: str | None
    mf6_workspace: Path | None
    mf6_model: str
    recharge_package: str
    recharge_map: str
    area_rtol: float | None
    wells: list[WellConfig]
    users: list[UserConfig]
    calibration_record: Path | None = None
    calibration_parameters: Path | None = None
    calibration_zones: Path | None = None
    calibration_pilot_points: dict | None = None  # {file, range, nugget}: MODFLOW 6 parameters by pilot point
    places_of_use: list[tuple[str, str]] = field(default_factory=list)
    point_sources: list[tuple[str, str]] = field(default_factory=list)
    exchange: dict | None = None
    utilities: list[dict] = field(default_factory=list)
    rules: list[dict] = field(default_factory=list)
    outputs: dict = field(default_factory=dict)
    nextgen_run: Path | None = None
    rivers: dict | None = None
    mf6_builder_args: dict = field(default_factory=dict)
    transport: dict[str, dict] = field(default_factory=dict)  # species -> recharge, stream (kg m-3) [A], model
    salt: dict = field(default_factory=dict)  # species -> salt.Settings: chloride through land, aquifer and streams
    salt_conductance: dict | None = None  # {base, slopes: {species: uS/cm per mg/L}} (salt.conductance)
    recharge_source: str = "aquifer"   # or "hru"
    recharge_add: bool = False          # a recharge package Basin Ledger adds, for a model without one
    weather: dict | None = None        # kind -> CSV path, and change (weather.py)
    groundwater_to_soil: dict | None = None  # rule, cap, package (coupling/gwsoil.py)
    mf6_length_units: str | None = None  # the model's, where its DIS leaves them undefined (coupling/modflow6.py)
    mf6_time_units: str | None = None
    mf6_models: list[str] = field(default_factory=list)  # several groundwater models coupled together, the first first
    mf6_storage: dict | None = None  # sy, ss (per meter), iconvert: storage for a model handed over steady
    irrigation: dict[str, list[tuple[str, float]]] | None = None  # HRU -> its wells and shares (coupling/irrigation.py)
    lakes: list = field(default_factory=list)  # LakeConfig: a SWAT+ reservoir coupled to a LAK lake (coupling/lakes.py)

    @property
    def surface(self) -> str:
        return "nextgen" if self.nextgen_run else "swatplus"

    @property
    def has_modflow6(self) -> bool:
        return bool(self.mf6_builder or self.mf6_workspace)

    @property
    def models(self) -> list[str]:
        """The hydrologic models a run of this example uses, for its calibration status."""
        return [self.surface] + (["modflow6"] if self.has_modflow6 else [])

    @property
    def folder(self) -> Path:
        return self.path.parent

    @property
    def calibration_status(self) -> str:
        """'calibrated' or 'screening' (no record, or a record that doesn't say calibrated)."""
        from .calibration import record_status
        return record_status(self.calibration_record, self.models)


def _path(base: Path, p: str) -> Path:
    if "{swatplus_refdata}" in p:
        p = p.replace("{swatplus_refdata}", str(swatplus_refdata()))
    p = p.replace("{nextgen_michigan}", str(nextgen_michigan()))
    p = p.replace("{cache}", str(cache_dir()))
    q = Path(p).expanduser()
    return q if q.is_absolute() else (base / q).resolve()


def load(path: str | Path) -> ExampleConfig:
    path = Path(path).resolve()
    with open(path, "rb") as f:
        c = tomllib.load(f)
    base = path.parent
    known = {"name", "description", "run", "swatplus", "modflow6", "recharge", "wells", "users", "calibration",
             "exchange", "utilities", "rules", "outputs", "nextgen", "rivers", "transport", "weather",
             "groundwater_to_soil", "irrigation", "lakes", "salt"}
    if set(c) - known:
        raise ValueError(f"{path}: unknown sections {sorted(set(c) - known)}")
    mf = c.get("modflow6", {})
    if "model" in mf and "models" in mf:
        raise ValueError(f"{path}: [modflow6] takes model (one) or models (several coupled together), not both")
    models = mf.get("models", [])
    if "models" in mf and (not isinstance(models, list) or not models or not all(isinstance(m, str) for m in models)):
        raise ValueError(f"{path}: [modflow6] models is a list of the groundwater models' names; got {models!r}")
    st = mf.get("storage")
    if st is not None and (not isinstance(st, dict) or set(st) - {"sy", "ss", "iconvert"}):
        raise ValueError(f"{path}: [modflow6] storage takes sy, ss (per meter) and iconvert; got {st!r}")
    if "modflow6" in c and bool(mf.get("builder")) == bool(mf.get("workspace")):
        raise ValueError(f"{path}: [modflow6] needs one of builder or workspace")
    from .coupling.modflow6 import DAYS, METERS
    for key, known in (("length_units", METERS), ("time_units", DAYS)):
        if key in mf and mf[key] not in known:
            raise ValueError(f"{path}: [modflow6] {key} = {mf[key]!r}; one of {', '.join(known)}")
    if ("swatplus" in c) == ("nextgen" in c):
        raise ValueError(f"{path}: needs one surface model, [swatplus] or [nextgen]")
    sw = c.get("swatplus", {})
    rch = c.get("recharge", {})
    if rch.get("source", "aquifer") not in ("aquifer", "hru"):
        raise ValueError(f"{path}: [recharge] source must be aquifer or hru")
    if rch.get("source") == "hru" and "nextgen" in c:
        raise ValueError(f"{path}: [recharge] source = hru needs SWAT+'s HRUs; a NextGen run gives catchments")
    weather = None
    if "weather" in c:
        from .weather import CHANGES, KINDS
        wx = c["weather"]
        if "nextgen" in c:
            raise ValueError(f"{path}: [weather] writes SWAT+'s gage records; a replayed NextGen run has its own forcing")
        if set(wx) - set(KINDS) - {"change"}:
            raise ValueError(f"{path}: [weather] unknown settings {sorted(set(wx) - set(KINDS) - {'change'})}")
        if set(wx.get("change", {})) - CHANGES:
            raise ValueError(f"{path}: [weather.change] unknown settings {sorted(set(wx['change']) - CHANGES)} "
                             f"(one of {sorted(CHANGES)})")
        weather = {k: _path(base, v) for k, v in wx.items() if k != "change"}
        weather["change"] = dict(wx.get("change", {}))
    gs = c.get("groundwater_to_soil")
    if gs is not None:
        from .coupling.gwsoil import CAPS, RULES
        missing = [k for k in ("rule", "cap", "package") if k not in gs]
        if missing:
            raise ValueError(f"{path}: [groundwater_to_soil] needs {missing}: the rule and the cap have no default")
        if set(gs) - {"rule", "cap", "package"}:
            raise ValueError(f"{path}: [groundwater_to_soil] unknown settings {sorted(set(gs) - {'rule', 'cap', 'package'})}")
        if gs["rule"] not in RULES or gs["cap"] not in CAPS:
            raise ValueError(f"{path}: [groundwater_to_soil] rule is one of {RULES}, cap one of {CAPS}")
        if rch.get("source") != "hru" or "exchange" not in c or "modflow6" not in c or "swatplus" not in c:
            raise ValueError(f"{path}: [groundwater_to_soil] needs SWAT+, MODFLOW 6, [recharge] source = \"hru\" (the "
                             "HRU map) and the two-way [exchange] (SWAT+'s aquifers passing recharge through)")
    transport = c.get("transport", {})
    salt = {}
    salt_conductance = None
    for sp, st in c.get("salt", {}).items():
        from .salt import Settings
        if sp == "conductance":
            salt_conductance = st
            continue
        if sp not in transport:
            raise ValueError(f"{path}: [salt.{sp}] carries {sp} through the aquifer's transport model: it needs "
                             f"[transport.{sp}]")
        if set(transport[sp]) & {"recharge", "stream"}:
            raise ValueError(f"{path}: [transport.{sp}]'s recharge and stream come from [salt.{sp}]; remove them")
        if "swatplus" not in c or "modflow6" not in c or rch.get("source") != "hru" or "exchange" not in c:
            raise ValueError(f"{path}: [salt.{sp}] needs SWAT+, MODFLOW 6, [recharge] source = \"hru\" (the "
                             "recharge's mass is placed HRU by HRU) and the two-way [exchange]")
        salt[sp] = Settings.from_config(sp, st, base)
    if salt_conductance is not None:
        if not isinstance(salt_conductance, dict) or set(salt_conductance) != {"base", "slopes"} \
                or not isinstance(salt_conductance["slopes"], dict):
            raise ValueError(f"{path}: [salt] conductance is {{base = uS/cm, slopes = {{<species> = uS/cm per mg/L}}}}")
        unknown = sorted(set(salt_conductance["slopes"]) - set(salt))
        if unknown or not salt_conductance["slopes"]:
            raise ValueError(f"{path}: [salt] conductance slopes name species with no [salt.<species>]: {unknown}")
        salt_conductance = {"base": float(salt_conductance["base"]),
                            "slopes": {k: float(v) for k, v in salt_conductance["slopes"].items()}}
    for sp, t in transport.items():
        if sp in salt:
            continue
        if not isinstance(t, dict) or "recharge" not in t:
            raise ValueError(f"{path}: [transport.{sp}] needs recharge (kg m-3, or \"swatplus\")")
        if set(t) - {"recharge", "stream", "model", "surface"}:
            raise ValueError(f"{path}: [transport.{sp}] unknown settings "
                             f"{sorted(set(t) - {'recharge', 'stream', 'model', 'surface'})}")
        if "stream" in t and "exchange" not in c:
            raise ValueError(f"{path}: [transport.{sp}] stream needs the two-way [exchange]")
        if "exchange" in c and "stream" not in t:
            raise ValueError(f"{path}: [transport.{sp}] needs stream (kg m-3, or \"swatplus\"): the coupled rivers' "
                             "losses carry it")
        for k in ("recharge", "stream"):
            v = t.get(k)
            if isinstance(v, str) and v != "swatplus" or isinstance(v, bool):
                raise ValueError(f"{path}: [transport.{sp}] {k} is a concentration (kg m-3) or \"swatplus\"")
        from_swat = "swatplus" in (t.get("recharge"), t.get("stream"))
        if from_swat or "surface" in t:
            t.setdefault("surface", sp)
            from .coupling.surface_mass import SPECIES
            if t["surface"] not in SPECIES:
                raise ValueError(f"{path}: [transport.{sp}] surface is SWAT+'s species, one of {SPECIES} "
                                 "(swatplus-bmi's names)")
            if "swatplus" not in c:
                raise ValueError(f"{path}: [transport.{sp}] carries SWAT+'s {t['surface']}: it needs [swatplus]")
        if any(isinstance(x, dict) and x.get("surface") == t.get("surface") and o != sp and t.get("surface")
               for o, x in transport.items()):
            raise ValueError(f"{path}: [transport.{sp}]: two species carry SWAT+'s {t['surface']}")
    irrigation = None
    if "irrigation" in c:
        from .coupling.irrigation import read_map
        ir = c["irrigation"]
        if set(ir) != {"wells"}:
            raise ValueError(f"{path}: [irrigation] takes wells (a list of {{hru, well, share}}, or a CSV file "
                             "with those columns)")
        if "swatplus" not in c or "modflow6" not in c:
            raise ValueError(f"{path}: [irrigation] pumps SWAT+'s irrigation from MODFLOW 6's wells; it needs both")
        irrigation = read_map(base, ir["wells"])
        named = {w["name"] for w in c.get("wells", [])}
        missing = sorted({w for ws in irrigation.values() for w, _ in ws} - named)
        if missing:
            raise ValueError(f"{path}: [irrigation] wells {missing} aren't in [[wells]]")
    lakes = []
    if c.get("lakes"):
        from .coupling.lakes import FORMS, LakeConfig
        if "swatplus" not in c or "modflow6" not in c:
            raise ValueError(f"{path}: [[lakes]] couple SWAT+'s reservoirs to MODFLOW 6's lakes; they need both")
        for lk in c["lakes"]:
            if set(lk) - {"reservoir", "package", "lake", "form", "model"} or "reservoir" not in lk or \
                    "package" not in lk:
                raise ValueError(f"{path}: each [[lakes]] has reservoir and package, and may have lake (1-based) and "
                                 f"form; got {sorted(lk)}")
            if lk.get("form", "active") not in FORMS:
                raise ValueError(f"{path}: [[lakes]] {lk['reservoir']}: form is one of {', '.join(FORMS)}")
            lakes.append(LakeConfig(lk["reservoir"], lk["package"], int(lk.get("lake", 1)), lk.get("form", "active"),
                                    lk.get("model")))
    ex_lakes = c.get("exchange", {}).get("lakes")
    if ex_lakes is not None and (not isinstance(ex_lakes, dict) or set(ex_lakes) != {"package", "map"}):
        raise ValueError(f"{path}: [exchange] lakes is {{package = \"<a RIV package>\", map = \"<a map>\"}}")
    for h in c.get("outputs", {}).get("heads", []):
        _cell(path, f"[outputs] heads {h.get('name')!r}", h.get("cell"))
    for r in c.get("rules", []):
        if r.get("kind") == "drawdown_limit":
            _cell(path, f"rule {r.get('name', r['kind'])!r}", r.get("cell"))
    users = []
    for u in c.get("users", []):
        users.append(UserConfig(
            name=u["name"], component=u["component"],
            config=_path(base, u["config"]) if u.get("config") else None,
            withdrawal=u["withdrawal"], returns=u["return"],
            variables=u.get("variables", {}), constituents=u.get("constituents"), forcing=u.get("forcing"),
            inputs=u.get("inputs", {}), report=u.get("report", {})))
    return ExampleConfig(
        path=path, name=c["name"], description=c.get("description", ""),
        run_folder=_path(base, c.get("run", {}).get("folder", f"runs/{c['name']}")),
        days=int(c.get("run", {}).get("days", 0)),
        swat_project=_path(base, sw["project"]) if sw else None,
        mf6_builder=mf.get("builder"),
        mf6_workspace=_path(base, mf["workspace"]) if mf.get("workspace") else None,
        mf6_model=models[0] if models else mf.get("model", "gwf"),
        mf6_models=list(models),
        mf6_storage=dict(st) if st is not None else None,
        recharge_package=rch.get("package", "rch"),
        recharge_map=rch.get("map", ""),
        recharge_add=_flag(path, "[recharge] add", rch.get("add", False)),
        area_rtol=None if rch.get("area_check", True) is False else float(rch.get("area_rtol", 1e-6)),
        wells=[_well(path, w) for w in c.get("wells", [])],
        users=users,
        calibration_record=_path(base, c["calibration"]["record"]) if c.get("calibration", {}).get("record") else None,
        calibration_parameters=(_path(base, c["calibration"]["parameters"])
                                if c.get("calibration", {}).get("parameters") else None),
        calibration_zones=_path(base, c["calibration"]["zones"]) if c.get("calibration", {}).get("zones") else None,
        calibration_pilot_points=_pilot_points(path, base, c.get("calibration", {}).get("pilot_points")),
        places_of_use=[(x["name"], x["channel"]) for x in sw.get("places_of_use", [])],
        point_sources=[(x["name"], x["channel"]) for x in sw.get("point_sources", [])],
        exchange=c.get("exchange"),
        utilities=c.get("utilities", []),
        rules=c.get("rules", []),
        outputs=c.get("outputs", {}),
        nextgen_run=_path(base, c["nextgen"]["run"]) if "nextgen" in c else None,
        rivers=c.get("rivers"),
        mf6_builder_args=mf.get("builder_args", {}),
        mf6_length_units=mf.get("length_units"),
        mf6_time_units=mf.get("time_units"),
        transport={sp: dict(t) for sp, t in transport.items()},
        salt=salt,
        salt_conductance=salt_conductance,
        recharge_source=rch.get("source", "aquifer"),
        weather=weather,
        groundwater_to_soil=dict(gs) if gs is not None else None,
        irrigation=irrigation,
        lakes=lakes,
    )
