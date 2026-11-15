"""The driver: steps every model of an example one day at a time, moves water between them, and keeps the ledger.

Day d, in order:

1. **Rules** read yesterday's state (SWAT+ channel flow, MODFLOW 6 heads) and set controls on users for today.
2. **Users** are stepped through the day (hourly) and their volumes integrated.
3. **Routing.** Each user's withdrawal and return go to their routes, directly or through a utility chain. Yesterday's
   stream-aquifer exchange is routed to SWAT+ too.
4. **SWAT+** takes today's demands and point sources and runs the day. Then what it actually withdrew and received is
   read back.
5. **MODFLOW 6** prepares the day, takes SWAT+'s recharge and today's well pumping, and solves. Then the wells and the
   river exchange are read back.
6. **The ledger** enters both sides of every transfer and checks them. The day's flows and heads go to the outputs.

**What Basin Ledger takes over in MODFLOW 6** (coupling/existing.py): on the first day, once the packages have their
entries, each package with flows is listed with the entries a link sets (recharge, the exchange, the soil link,
users' wells) and the rest, the model's, booked as water from or to outside the coupled models
(``Basin.takeover``). A model's own package that would count water twice (a second recharge package, UZF, EVT with
groundwater up into the soil) stops the start-up.

**Transport** (``[transport.<species>]``, coupling/transport.py): after MODFLOW 6 prepares the day, the recharge's and
the coupled rivers' losses' concentrations are set [A]; after the solve, each species' mass is booked (kg) and the
concentration at every configured well becomes the intake concentration
(``water_user_intake_water_<species>__mass_concentration``) of each user withdrawing from it, for the next day. The
users run before MODFLOW 6 each day, so their intake is the aquifer's at the end of the day before: a one-day lag, the
same as the exchange's. On the first day it is the starting concentration.

``save(folder)`` writes the state at the pause after a day: SWAT+'s saved state (swatplus-bmi), MODFLOW 6's heads and
concentrations, the ledger, and the exchange due. ``Basin(..., restart=folder)`` carries on from there. SWAT+ and
MODFLOW 6 continue exactly, and the users are stepped forward to the restart day.
"""

from __future__ import annotations

import datetime as dt
import importlib.util
import json
import math
import shutil
import warnings
from pathlib import Path

import numpy as np
from swatplus_bmi import SwatPlusBmi
from swatplus_bmi.exchange import ExchangeError

from . import config as C
from . import rules as R
from .chain import UtilityChain
from .coupling import existing as EX
from .coupling import irrigation as IR
from .coupling import lakes as LK
from .coupling import recharge as RC
from .coupling import swatplus as S
from .coupling.exchange import ChannelMap, StreamExchange, StreamExchanges, package_type
from .coupling.exchange import requests as exchange_requests
from .coupling.gwsoil import GroundwaterToSoil
from .coupling.modflow6 import Models, Modflow6, remaining_simulation
from .coupling.recharge import RechargeLink
from .coupling.timesteps import daily_steps, start_date
from .ledger import FLOAT32, Ledger
from .mapping import AreaMap
from .routes import Router, kind
from .units import DAY
from .users import WaterUser, load_component


def _example_function(cfg: C.ExampleConfig, spec: str, default: str):
    mod, _, fn = spec.partition(":")
    name = f"basinledger_example_{cfg.name}_{mod}".replace("/", "_").replace("-", "_")
    spec_ = importlib.util.spec_from_file_location(name, cfg.folder / f"{mod}.py")
    m = importlib.util.module_from_spec(spec_)
    spec_.loader.exec_module(m)
    return getattr(m, fn or default)


class Basin:
    def __init__(self, cfg: C.ExampleConfig, run: str = "base", users: bool = True, days: int | None = None,
                 restart: str | Path | None = None, prepare=None, heads: np.ndarray | None = None,
                 wrap_surface=None, parameters: list | None = None, _problems: list | None = None):
        """``prepare(model, folder)``, if given, is called on the run's copy of each model's files ("swatplus",
        "modflow6") before that model starts, after the parameter set is written into them: the example's own
        (``[calibration] parameters``) with ``parameters``' values in place of any it names (calibration.py's
        ``apply_parameters``). ``heads``: MODFLOW 6's
        starting heads (per reduced node), in place of its input's, as from ``spin_up``. ``wrap_surface(model)``,
        if given, returns what the driver steps in place of the surface model once it has started: a study's way to
        hold or replace one of its values through BMI each day (a process substitution), with everything else
        passed through."""
        self.cfg = cfg
        self._problems = _problems  # a check (Basin.check): problems listed here, each part carrying on
        self.status = cfg.calibration_status  # 'screening' unless the example's record says 'calibrated'
        from .calibration import apply_parameters, merged, read_parameters
        own = read_parameters(cfg.calibration_parameters) if cfg.calibration_parameters else []
        self.parameters = merged(own, parameters or [])
        unused = sorted({q.model for q in self.parameters} - set(cfg.models) - {"user"})
        if unused:
            self._problem(ValueError(f"parameters for {', '.join(unused)}, which this example doesn't run "
                                     f"({cfg.models})"))
        user_prepare = prepare

        def prepare(model, folder):
            apply_parameters(model, folder, self.parameters, cfg.calibration_zones, cfg.mf6_model)
            if user_prepare:
                user_prepare(model, folder)

        self.folder = cfg.run_folder / run
        if self.folder.exists():
            shutil.rmtree(self.folder)
        self.folder.mkdir(parents=True)
        saved = json.loads((Path(restart) / "basin.json").read_text()) if restart else None

        # the surface model: SWAT+ from a copy of the project with the objects the run needs, or a NextGen run
        ex = cfg.exchange
        swat_dir = None
        self.direct_exchange = False
        if cfg.surface == "swatplus":
            swat_dir = self.folder / "swatplus"
            shutil.copytree(cfg.swat_project, swat_dir)
            channels = S._channels(swat_dir)
            ex_channels = self._exchange_channels(ex) if ex else []
            self.direct_exchange = bool(ex) and S.has_groundwater_exchange()
            via_routes = [] if self.direct_exchange else ex_channels  # before 0.2.0: a point source and a seep per channel
            places = list(cfg.places_of_use) + [(f"seep_{c}", c) for c in via_routes]
            points = list(cfg.point_sources) + [(f"gw_{c}", c) for c in via_routes]
            if places:
                S.add_places_of_use(swat_dir, places)
            if points:
                S.add_point_sources(swat_dir, points)
            if ex:
                S.aquifer_pass_through(swat_dir)
                S.no_transmission_loss(swat_dir)
            from .temperature import print_daily
            if cfg.outputs.get("temperature"):  # the reach heat balance reads the routing units' daily water
                print_daily(swat_dir, ["lsunit_wb", "channel_sd"])
            if cfg.outputs.get("balance"):      # the whole-basin balance reads SWAT+'s daily land budget
                print_daily(swat_dir, ["basin_wb", "channel_sd", "reservoir"])
            prepare("swatplus", swat_dir)
            problems = S.coupling_problems(swat_dir, two_way=bool(ex))
            if problems:
                self._problem(ValueError("SWAT+ isn't set up to be coupled:\n- " + "\n- ".join(problems)))
            for w in ([] if saved else S.start_warnings(swat_dir)):  # a restart's soils are the saved state's
                warnings.warn(w, S.StartWarning, stacklevel=2)
            if cfg.weather and not saved:  # the days SWAT+ reads before it first pauses, into its files (weather.py)
                from . import weather as W
                W.write_first_days(swat_dir, cfg.weather, S._period(swat_dir)[0])
            self.swat = SwatPlusBmi()
            try:
                self.swat.initialize(str(swat_dir))
            except ExchangeError as e:
                if "protocol" not in str(e):
                    raise
                raise ExchangeError(f"{e}: SWAT+ was built by another swatplus-bmi; build it again with this one "
                                    f"(uv run swatplus-bmi build --jobs 2)") from e
            if wrap_surface is not None:
                self.swat = wrap_surface(self.swat)
        else:
            from .coupling.nextgen import NgenReplay
            if ex or cfg.places_of_use or cfg.point_sources:
                raise ValueError("a replayed NextGen run takes no surface routes or two-way exchange")
            self.swat = NgenReplay()
            self.swat.initialize(str(cfg.nextgen_run))
            channels = list(self.swat.ids)
            ex_channels = []
        self.sw = S.Swat(self.swat)
        units = self.swat.get_time_units()  # days since YYYY-MM-DD hh:mm:ss
        self.start = dt.date.fromisoformat(units.split(" since ")[1].split()[0])
        self.days = days or cfg.days or int(self.swat.get_end_time() - self.swat.get_start_time())
        self.day = 0
        if saved:
            self.swat.load_state(str(Path(restart) / "surface"))
            self.day = saved["day"]
        # weather into SWAT+'s gage records, from the first day it hasn't read (weather.py): a fresh run's first days
        # are already in its files; a restart keeps the warm-up's records up to then. A change is on the project's
        # own records, not the run's copy, whose first days hold it already
        self.weather = {}
        if cfg.weather:
            from . import weather as W
            last = self.start + dt.timedelta(days=min(self.days, int(self.swat.get_end_time())) - 1)
            self.weather = W.apply(self.swat, cfg.swat_project, cfg.weather, last,
                                   written=None if saved else self.start)

        # MODFLOW 6 (optional): built (or copied) for the run, or for the days left after a restart
        self.mf6 = None
        if cfg.has_modflow6:
            if cfg.surface == "nextgen" and "soil_water_table__recharge_rate" not in self.swat.get_output_var_names():
                raise NotImplementedError("MODFLOW 6 with NextGen needs CFE's percolation (SOIL_TO_GW_FLUX) in the "
                                          "run's catchment outputs; see coupling/nextgen.py")
            ws = self.folder / "mf6"
            if cfg.mf6_builder:
                _example_function(cfg, cfg.mf6_builder, "build")(ws, self.days - self.day,
                                                                 self.start + dt.timedelta(days=self.day),
                                                                 **cfg.mf6_builder_args)
            else:
                shutil.copytree(cfg.mf6_workspace, ws)
            # every step one day, each stress period kept (coupling/timesteps.py), then a restart's days cut off
            for note in daily_steps(ws, cfg.mf6_time_units):
                warnings.warn(note, EX.TakeoverNote, stacklevel=2)
            several = len(cfg.mf6_models) > 1  # several groundwater models coupled together (modflow6.Models)
            if several:
                for what, on in (("transport", cfg.transport), ("[rivers]", cfg.rivers),
                                 ("[groundwater_to_soil]", cfg.groundwater_to_soil), ("[recharge] add", cfg.recharge_add)):
                    if on:
                        raise ValueError(f"{what} with several groundwater models coupled isn't built")
            by_model: dict[str, list] = {}  # new users' wells, a package of its own in each model that has them
            for w in cfg.wells:
                if w.add:
                    by_model.setdefault((w.model or cfg.mf6_model).lower(), []).append((w.name, w.package, w.cell))
            for model, added in by_model.items():
                for note in EX.add_wells(ws, model, added):
                    warnings.warn(note, EX.TakeoverNote, stacklevel=2)
            if cfg.recharge_add:  # a model without recharge: a package of Basin Ledger's own on the map's cells
                for note in EX.add_recharge(ws, cfg.mf6_model, cfg.recharge_package, self._map_cells(cfg)):
                    warnings.warn(note, EX.TakeoverNote, stacklevel=2)
            if saved and not cfg.mf6_builder:
                remaining_simulation(ws, ws, self.day, cfg.mf6_time_units)
            elif not saved:
                first = start_date(ws)
                if first is not None and first != self.start:
                    warnings.warn(f"MODFLOW 6's time starts on {first} (TDIS START_DATE_TIME) and the run on "
                                  f"{self.start}: the model's first day is taken as the run's, so its stress "
                                  f"periods fall {abs((first - self.start).days):,} days "
                                  f"{'earlier' if first > self.start else 'later'} than its dates",
                                  EX.TakeoverNote, stacklevel=2)
            prepare("modflow6", ws)
            self.mf6 = Modflow6(ws, cfg.mf6_model, length_units=cfg.mf6_length_units,
                                time_units=cfg.mf6_time_units)
            if not self.mf6.si:
                warnings.warn(f"MODFLOW 6 in {self.mf6.units_summary()}", EX.TakeoverNote, stacklevel=2)
            if several:
                try:
                    self.mf6 = Models(self.mf6, cfg.mf6_models)
                except ValueError:
                    self.mf6.finalize()
                    raise
            if saved:
                self.mf6.set_heads(np.load(Path(restart) / "heads.npy"))
                # what a restart doesn't carry: packages holding water from day to day start from their files
                held = {"maw": "its wells' water", "uzf": "the unsaturated zone's water", "lak": "the lakes' stages",
                        "csub": "compaction"}
                coupled_lakes = {lk.package.upper() for lk in cfg.lakes}  # their stages come from SWAT+ each day
                for v in self._views():
                    for p in EX.flow_packages(v, v.model):
                        t = str(v.get("PACKAGE_TYPE", v.model, p)[0]).strip().lower()
                        if t in held and p.upper() not in coupled_lakes:
                            warnings.warn(f"{p.lower()}: a restart carries MODFLOW 6's heads but not {held[t]} "
                                          f"({t.upper()}), which start from the model's files", EX.TakeoverNote,
                                          stacklevel=2)
            elif heads is not None:
                self.mf6.set_heads(np.asarray(heads, float))
            if self.mf6.end < self.days - self.day - 1e-9:
                self._problem(ValueError(f"MODFLOW 6 simulation is {self.mf6.end} days; the run needs "
                                         f"{self.days - self.day}"))
        elif cfg.wells or ex or cfg.outputs.get("heads") or cfg.rivers or cfg.transport:
            raise ValueError("wells, rivers, the two-way exchange, transport and head outputs need a [modflow6] model")
        self.transport = {}
        if self.mf6 is not None:
            self.transport = self._transport_models(cfg)
            for sp, t in self.transport.items():
                if saved:
                    t.set_concentration(np.load(Path(restart) / f"concentration_{sp}.npy"))
        self.rivers = None  # built on the first day, once MODFLOW 6's river package has its entries
        if cfg.rivers and cfg.rivers.get("stage", "fixed") not in ("fixed", "surface_depth", "nextgen_depth"):
            self._problem(ValueError("[rivers] stage must be fixed or surface_depth"))
        if cfg.rivers and cfg.rivers.get("stage") in ("surface_depth", "nextgen_depth") and not self.sw.has(S.DEPTH):
            self._problem(ValueError("[rivers] stage = surface_depth needs the surface model's channel depth "
                                     "(NextGen's t-route, or swatplus-bmi 0.2.0)"))

        self.ledger = Ledger()
        L = self.ledger
        # recharge: SWAT+ aquifers (or HRUs) -> MODFLOW 6 cells
        self.aquifer_check = None
        self.gw_soil = None
        self.recharge = None
        if self.mf6 is not None:
            self._part("[recharge]", lambda: self._recharge(cfg, swat_dir, saved))
        wells = {}
        self._wells_left_out = set()  # a check's wells with a problem: listed once, not again by the users on them
        for w in cfg.wells:
            if self._part(f"well {w.name}", lambda w=w: wells.__setitem__(w.name, (w.package, self._well_node(w)))
                          or True) is None:
                self._wells_left_out.add(w.name)
        self.router = Router(self.sw, self.mf6, wells, L)

        # irrigation pumped from MODFLOW 6's wells (coupling/irrigation.py): the HRUs mapped and their wells
        self.irrigation = None
        if cfg.irrigation:
            self._part("[irrigation]", lambda: self._irrigation(cfg, swat_dir, wells))
        # lakes: SWAT+'s reservoirs holding MODFLOW 6's lakes' water (coupling/lakes.py)
        self.lakes = None
        if cfg.lakes:
            self._part("[[lakes]]", lambda: self._lakes(cfg, swat_dir, saved))

        # stream-aquifer exchange (built on the first day, once MODFLOW 6's river package has its entries)
        self.exchange = None
        self.due = {c: 0.0 for c in ex_channels}  # yesterday's exchange, m3 per channel, applied in SWAT+ today
        if saved and saved.get("exchange_due"):
            self.due = dict(saved["exchange_due"])
        if ex:
            for c in ex_channels:
                if not self.direct_exchange:
                    self._part("[exchange]", lambda c=c: self.router.check(f"swat:point_source:gw_{c}"))
                    self._part("[exchange]", lambda c=c: self.router.check(f"swat:place_of_use:seep_{c}"))
                k = f"exchange:{c}:"
                L.pair(k + "due", [k + "applied", k + "shortfall"], rtol=FLOAT32)
            L.zero("swatplus:aquifer_outflow", atol=1e-9)

        # utilities, users, rules
        self.utilities = {u["name"]: UtilityChain.from_config(u) for u in cfg.utilities}
        for u in self.utilities.values():
            u.declare(L)
            for route in [s.route for s in u.sources] + [u.outfall]:
                if route.removeprefix("well:") not in self._wells_left_out:
                    self._part(f"utility {u.name}", lambda r=route: self.router.check(r))
        self.users = []
        if users:
            for u in cfg.users:
                self._part(f"user {u.name}", lambda u=u: self._add_user(u, L))
        # intakes from the aquifer: each user on a configured well, for each species its component takes
        self.intakes = []  # (user, WaterUser, input variable, species, node)
        for u, wu in self.users:
            k, name = kind(u.withdrawal)
            if k != "well" or name not in wells:
                continue
            takes = set(wu.bmi.get_input_var_names())
            for sp in self.transport:
                var = f"water_user_intake_water_{sp}__mass_concentration"
                if var not in takes:
                    continue
                if var in (u.inputs or {}):
                    self._problem(ValueError(f"user {u.name}: {var} is set from the aquifer ([transport.{sp}]); "
                                             "remove its constant input"))
                    continue
                self.intakes.append((u.name, wu, var, sp, wells[name][1]))
        self._set_intakes()
        users_by_name = {u.name: wu for u, wu in self.users}
        self.rules = []
        for r in cfg.rules:
            self._part(f"rule {r.get('name', r['kind'])}", lambda r=r: self._add_rule(r, users, users_by_name))
        self.rule_memory = {}

        # outputs
        out = cfg.outputs
        self.out_channels = out.get("channels", channels)
        self._ch_index = {c: channels.index(c) for c in self.out_channels}
        self.out_heads = []
        for h in out.get("heads", []):
            self._part(f"head output {h['name']}",
                       lambda h=h: self.out_heads.append((h["name"], self.mf6.node(*h["cell"], model=h.get("model")))))
        self.series: dict[str, list[float]] = {f"flow:{c}": [] for c in self.out_channels}
        # NextGen: flow at each nexus, from the flowpaths draining to it (the hydrofabric's topology)
        self._nexus: dict[str, list[int]] = {}
        for c, nx in getattr(self.swat, "nexus", {}).items():
            if c in self._ch_index and nx:
                self._nexus.setdefault(nx, []).append(self._ch_index[c])
        self.series.update({f"flow:{nx}": [] for nx in self._nexus})
        if cfg.rivers:
            for c in self._river_objects(cfg.rivers):
                self.series[f"baseflow:{c}"] = []
        self.series.update({f"head:{n}": [] for n, _ in self.out_heads})

        if saved:
            L.load(saved["ledger"])
            self.series = saved["series"]
            self.rule_memory = saved["rule_memory"]
            for _, wu in self.users:
                wu.fast_forward(self.day)

    # ---------- start-up, and checks ----------
    def _problem(self, e: Exception):
        """A problem at start-up: raised in a run, listed in a check, which carries on."""
        if self._problems is None:
            raise e
        self._problems.append(e.args[0] if isinstance(e, KeyError) and e.args else str(e))

    def _part(self, where: str, fn):
        """One part of the start-up: in a run it is just called, and a problem stops the run; in a check (Basin.check)
        a problem is listed under ``where`` and the start-up goes on without that part."""
        if self._problems is None:
            return fn()
        try:
            return fn()
        except Exception as e:  # noqa: BLE001 (a check lists every kind of problem)
            msg = e.args[0] if isinstance(e, KeyError) and e.args else str(e)
            self._problems.append(msg if msg.startswith(where) else f"{where}: {msg}")
            return None

    def _recharge(self, cfg, swat_dir, saved):
        """Recharge from the surface model's aquifers or HRUs into MODFLOW 6's cells, and groundwater back up into the
        soil (its SWAT+ side; MODFLOW 6's on the first day)."""
        L = self.ledger
        var = RC.SOURCES[cfg.recharge_source]
        if not self.sw.has(var):
            raise ValueError(f"[recharge] source = {cfg.recharge_source} needs the surface model's {var} "
                             "(swatplus-bmi 0.6.0 for HRU percolation)")
        objects = self.sw.names(var)
        self.map = self._area_map(cfg.recharge_map, objects)
        if cfg.recharge_source == "hru":
            areas = S.hru_areas(swat_dir)
        else:
            areas = S.aquifer_areas(swat_dir) if swat_dir else self.swat.areas
        source = np.array([areas[a] for a in objects])
        self.map.check(self.mf6.cell_area, source, cfg.area_rtol)
        if cfg.recharge_source == "hru":
            self.ground_map = self.map  # each HRU's ground in each cell, for groundwater up into the soil
            # each HRU's rows scaled to its area in SWAT+: the map places, never adds or loses (mapping.py)
            self.map = self.map.conserving(source)
            aquifers = self.sw.names(RC.RECHARGE) if self.sw.has(RC.RECHARGE) else []
            if aquifers:  # SWAT+'s own aquifer recharge, the same percolation summed by aquifer (recharge.py)
                a = S.aquifer_areas(swat_dir)
                self.aquifer_check = np.array([a[x] for x in aquifers])
                L.pair("recharge:made", "recharge:swat_aquifers", rtol=RC.AQUIFER_RTOL)
        self.recharge = RechargeLink(self.swat, self.mf6, self.map, cfg.recharge_package, source, var)
        # groundwater up into the soil: SWAT+'s side now, MODFLOW 6's on the first day (coupling/gwsoil.py)
        if cfg.groundwater_to_soil:
            gs = cfg.groundwater_to_soil
            self.gw_soil = GroundwaterToSoil(self.sw, self.ground_map, source, gs["rule"], gs["cap"],
                                             gs["package"], saved.get("gw_soil_due") if saved else None)
            L.pair("gw_soil:due", "gw_soil:received", rtol=FLOAT32)
        L.pair("recharge:sent", "recharge:received")
        # what the surface model made, and what the map added to it (+) or couldn't place (-): water from outside
        L.pair("recharge:sent", ["recharge:made", "recharge:by_the_map"])

    def _add_user(self, u, L):
        comp = load_component(u.component)
        wu = WaterUser(u.name, comp, str(u.config) if u.config else None,
                       variables=u.variables or None, constituents=u.constituents,
                       report=u.report or None)
        wu.inputs = self._forcing(u, comp)  # after initialize, when a component's inputs are known
        for route in (u.withdrawal, u.returns):
            k, name = kind(route)
            if k == "utility":
                if name not in self.utilities:
                    raise KeyError(f"user {u.name}: no utility named {name!r}")
            elif not (k == "well" and name in self._wells_left_out):
                self.router.check(route)
        self.users.append((u, wu))
        k = f"user:{u.name}:"
        L.pair(k + "withdrawal", [k + "consumption", k + "return", k + "storage_less_gain"])

    def _add_rule(self, r, users, users_by_name):
        settings = {k: v for k, v in r.items() if k not in ("name", "kind")}
        fn = R.BUILT_IN.get(r["kind"]) or _example_function(self.cfg, r["kind"], "rule")
        if users and settings.get("user") not in users_by_name:
            raise KeyError(f"rule {r.get('name', r['kind'])}: no user {settings.get('user')!r}")
        need = {"curtailment": "limit", "drawdown_limit": "limit", "dry_mode": "dry_mode"}.get(r["kind"])
        if users and need and need not in users_by_name[settings["user"]].control_vars:
            raise KeyError(f"rule {r.get('name')}: user {settings['user']!r} takes no {need!r} input")
        if "cell" in settings and self.mf6 is not None:
            self.mf6.node(*settings["cell"])  # the cell is checked now, not on the first day it's read
        self.rules.append(R.Rule(r.get("name", r["kind"]), fn, settings))

    def _irrigation(self, cfg, swat_dir, wells):
        if not self.sw.has(IR.IRRIGATION):
            raise ValueError("[irrigation] reads each HRU's irrigation (land_irrigation_water__volume_flow_rate, "
                             "swatplus-bmi 0.11.0)")
        bad, notes = IR.problems(swat_dir, cfg.irrigation, set(wells))
        if bad:
            raise ValueError("\n- ".join(["irrigation can't be pumped as mapped:"] + bad))
        for n in notes:
            warnings.warn(n, EX.TakeoverNote, stacklevel=3)
        for ws in cfg.irrigation.values():
            for w, _ in ws:
                self.router.check(f"well:{w}")
        used = sorted({w for ws in cfg.irrigation.values() for w, _ in ws})
        self.ledger.pair("irrigation:pumped", [f"irrigation:well:{w}" for w in used])
        names = list(self.sw.names(IR.IRRIGATION))
        self.irrigation = [(names.index(h), ws) for h, ws in cfg.irrigation.items()]
        self._irrigation_rest = np.ones(len(names), bool)
        for i, _ in self.irrigation:
            self._irrigation_rest[i] = False

    def _lakes(self, cfg, swat_dir, saved):
        if not self.sw.has(LK.VOLUME) or not self.sw.has(LK.GW):
            raise ValueError("[[lakes]] read and set SWAT+'s reservoirs (reservoir_water__volume and its groundwater "
                             "exchange, swatplus-bmi 0.11.0)")
        if self.mf6 is None:
            raise ValueError("[[lakes]] need a [modflow6] model")
        bad = LK.problems(self.mf6, cfg.lakes, list(self.sw.names(LK.VOLUME)), swat_dir)
        if bad:
            raise ValueError("\n- ".join(["the lakes can't be coupled as configured:"] + bad))
        self.lakes = LK.LakeLink(self.mf6, self.sw, cfg.lakes)
        if saved and saved.get("lakes_due"):
            self.lakes.due = dict(saved["lakes_due"])
        for lk in cfg.lakes:
            k = f"lake:{lk.reservoir}:"
            self.ledger.pair(k + "due", [k + "applied", k + "shortfall"], rtol=FLOAT32)

    def _lakes_day(self, L):
        """Each lake's day: yesterday's leakage as SWAT+'s reservoir took it, and today's leakage, due tomorrow."""
        form = {lk.reservoir: lk.form for lk in self.cfg.lakes}
        for res, d in self.lakes.after_solve().items():
            k = f"lake:{res}:"
            short = d["due"] - d["applied"] if d["due"] < 0 and d["applied"] > d["due"] else 0.0  # a loss it lacked
            L.enter(k + "due", d["due"])
            L.enter(k + "applied", d["applied"])
            L.enter(k + "shortfall", short)
            L.enter(k + "to_aquifer", d["to_aquifer"])       # today's, out of the lake into MODFLOW 6
            L.enter(k + "stage_m", d["stage"])
            L.enter(k + "volume", d["volume"])
            if form[res] == "active":  # where LAK's stage ended the day; not kept, tomorrow starts from SWAT+'s
                L.enter(k + "mf6_end_stage_m", d["mf6_end_stage"])

    def _irrigate(self, L, rt):
        """The day's irrigation on the mapped HRUs, requested from their wells; the rest stays water from outside."""
        v = self.sw.get(IR.IRRIGATION) * DAY
        pumped = 0.0
        for i, ws in self.irrigation:
            for w, share in ws:
                rt.request(f"well:{w}", float(v[i]) * share)
                L.enter(f"irrigation:well:{w}", float(v[i]) * share)
            pumped += float(v[i])
        L.enter("irrigation:pumped", pumped)
        L.enter("irrigation:from_outside", float(v[self._irrigation_rest].sum()))

    def _exchange_channels(self, ex) -> list[str]:
        """The channels the exchange's rivers and drains couple, in order; needed before SWAT+ starts."""
        import csv
        out = []
        for spec in [ex.get("map", "")] + ([ex["drains"]["map"]] if ex.get("drains") else []):
            if spec.startswith("all:"):
                out.append(spec[4:])
                continue
            with open(self.cfg.folder / spec, newline="") as f:
                out += [r["object"] for r in csv.DictReader(f)]
        return list(dict.fromkeys(out))

    def _river_objects(self, rv) -> list[str]:
        spec = rv.get("map", "")
        if spec.startswith("all:"):
            return [spec[4:]]
        import csv
        with open(self.cfg.folder / spec, newline="") as f:
            return list(dict.fromkeys(r["object"] for r in csv.DictReader(f)))

    def _rivers_day(self):
        """Set today's river stages (if from NextGen's depth), before the solve."""
        rv = self.cfg.rivers
        if self.rivers is None:
            self._build_rivers()
        pkg, cmap, rbot, stage0 = self.rivers
        if rv.get("stage") in ("surface_depth", "nextgen_depth"):
            depth = dict(zip(self.sw.names(S.DEPTH), self.sw.get(S.DEPTH)))  # yesterday's (SWAT+) or today's (replay)
            e = cmap.entry_channel
            stage = stage0.copy()
            for i, c in enumerate(cmap.channels):
                stage[e == i] = rbot[e == i] + float(rv.get("bed_thickness", 1.0)) + depth[c]
            self.mf6.set_stage(pkg, stage)

    def _build_rivers(self):
        """The river cells by flowpath, once MODFLOW 6 has read its first day."""
        rv = self.cfg.rivers
        pkg = rv.get("package", "riv")
        if package_type(self.mf6, pkg) != "riv":
            raise ValueError(f"[rivers] sets a river (RIV) package's stage; {pkg!r} is "
                             f"{package_type(self.mf6, pkg).upper()}")
        cmap = ChannelMap.from_spec(rv.get("map", ""), self.cfg.folder, self.mf6, pkg, self.sw.names(S.QOUT))
        self.rivers = (pkg, cmap, self.mf6.get("RBOT", self.mf6.model, pkg.upper()).astype(float),
                       self.mf6.get("STAGE", self.mf6.model, pkg.upper()).astype(float))

    def _rivers_out(self):
        pkg, cmap, _, _ = self.rivers
        q = self.mf6.simvals(pkg)  # m3 d-1, positive into the aquifer
        for i, c in enumerate(cmap.channels):
            self.series[f"baseflow:{c}"].append(-float(q[cmap.entry_channel == i].sum()) / DAY)

    def _forcing(self, u: C.UserConfig, comp) -> dict:
        """Hourly weather inputs for a user, from its forcing source, for the inputs its component takes, and its
        constant inputs."""
        takes = set(comp.get_input_var_names())
        out = {}
        if u.forcing:
            start = dt.datetime.combine(self.start, dt.time())
            src = u.forcing.get("source")
            if src == "nextgen" and self.cfg.surface == "nextgen":
                from .forcing import nextgen_inputs
                series = nextgen_inputs(self.cfg.nextgen_run, u.forcing["catchment"], start)
            elif src == "swatplus" and self.cfg.surface == "swatplus":
                from .forcing import swatplus_inputs
                series = swatplus_inputs(self.folder / "swatplus", u.forcing["station"], start, self.weather)
            else:
                raise NotImplementedError(f"user {u.name}: forcing from {src!r} with a {self.cfg.surface} surface "
                                          "model (nextgen with NextGen, swatplus with SWAT+; AORC is planned)")
            out.update({k: v for k, v in series.items() if k in takes})
        for k, v in (u.inputs or {}).items():
            if k not in takes:
                raise KeyError(f"user {u.name}: {k!r} is not an input of {u.component}")
            out[k] = (lambda value: (lambda t: value))(float(v))
        return out

    @staticmethod
    def _map_cells(cfg: C.ExampleConfig) -> list[tuple[int, ...]]:
        """The recharge map's cells, 0-based, in the grid's form, for a recharge package Basin Ledger adds."""
        import csv
        if cfg.recharge_map.startswith("all:") or not cfg.recharge_map:
            raise ValueError("[recharge] add = true needs a map of cells (a CSV), not \"all:\"")
        out = []
        with open(cfg.folder / cfg.recharge_map, newline="") as f:
            for r in csv.DictReader(f):
                if (r.get("row") or "").strip():
                    out.append((int(r["layer"]), int(r["row"]), int(r["col"])))
                elif (r.get("icell2d") or "").strip():
                    out.append((int(r["layer"]), int(r["icell2d"])))
                else:
                    raise ValueError(f"{cfg.recharge_map}: [recharge] add = true needs cells, not nodes")
        return out

    def _area_map(self, spec: str, aquifers: list[str]) -> AreaMap:
        if spec.startswith("all:"):
            nodes = self.mf6.package_nodes(self.cfg.recharge_package)
            return AreaMap.all_cells(aquifers, spec[4:], nodes, self.mf6.cell_area)
        return AreaMap.from_csv(self.cfg.folder / spec, aquifers, self.mf6.node)

    def _transport_models(self, cfg: C.ExampleConfig) -> dict:
        """The configured species' transport models, and a check that the simulation has no others: a transport
        model the driver doesn't know would be solved with its files' concentrations and its mass left out."""
        from .coupling.transport import Transport
        names = set(self.mf6.api.get_input_var_names()) | set(self.mf6.api.get_output_var_names())
        found = {n.split("/")[0] for n in names if n.endswith("/SSM/IAUXPAK")}
        out = {sp: Transport(self.mf6, sp, t.get("model")) for sp, t in cfg.transport.items()}
        extra = found - {t.model for t in out.values()}
        if extra:
            raise ValueError(f"MODFLOW 6 has transport models with no [transport] section: {sorted(extra)}")
        return out

    def _set_intakes(self):
        """Each intake user's concentration for the coming day: the aquifer's at its well now."""
        self.intake = {}
        x = {sp: t.concentration() for sp, t in self.transport.items()}
        for user, wu, var, sp, node in self.intakes:
            c = float(x[sp][node])
            wu.inputs[var] = (lambda value: (lambda t: value))(c)
            self.intake[(user, sp)] = c

    def _build_mass_links(self):
        """The transport's links, once MODFLOW 6 has read its first day: recharge, the wells and the streams."""
        rp = self.cfg.recharge_package
        wells = {w.name: (w.package, self.mf6.entry(w.package, self.mf6.node(*w.cell, model=w.model)))
                 for w in self.cfg.wells}
        streams = []
        if self.exchange is not None:
            streams.append((self.exchange.package, self.exchange.map.entry_channel, self.exchange.map.channels))
            if self.exchange.drains is not None:
                dp, dmap = self.exchange.drains
                streams.append((dp, dmap.entry_channel, dmap.channels))
        if self.exchange is not None and self.exchange.kind == "sfr":
            raise NotImplementedError("transport with an SFR exchange needs MODFLOW 6's SFT package for the "
                                      "streams' concentrations; not built")
        for sp, t in self.transport.items():
            if not t.takes_inflow(rp):
                raise ValueError(f"{t.model.lower()}: the recharge package {rp!r} isn't an AUX source in SSM "
                                 "(add_transport's sources)")
            if self.exchange is not None and not t.takes_inflow(self.exchange.package):
                raise ValueError(f"{t.model.lower()}: the coupled river package {self.exchange.package!r} isn't an "
                                 "AUX source in SSM, so its losses would carry nothing")
        self._mass_links = (rp, wells, streams)
        L = self.ledger
        for sp in self.transport:
            k = f"transport:{sp}:"
            L.pair(k + "recharge:sent", k + "recharge:received")
            for c in (self.exchange.map.channels if self.exchange is not None else []):
                L.pair(k + f"stream:{c}:leakage_sent", k + f"stream:{c}:leakage_received")
            L.zero(k + "closure_share", atol=1e-6)

    def _transport_day(self):
        """The day's inflow concentrations, after MODFLOW 6 prepares the day [A]. The couplings are found on the first
        day, once the packages have their entries."""
        if not hasattr(self, "_mass_links"):
            self._build_mass_links()
        rp, _, streams = self._mass_links
        for sp, t in self.transport.items():
            t.set_inflow(rp, float(self.cfg.transport[sp]["recharge"]))
            if self.exchange is not None:
                t.set_inflow(self.exchange.package, float(self.cfg.transport[sp]["stream"]),
                             entries=self.exchange.map.entry_channel >= 0)

    def _transport_accounts(self, L):
        """Each species' mass for the day, kg: in with recharge (as sent: the stated concentration times the volume
        sent; as received: what the transport model took in), out with each well, to and from each coupled channel,
        across the boundaries not coupled (net, positive into the aquifer), the release from storage, and the
        closure share, checked every day. Then the wells' concentrations, for tomorrow's intakes."""
        rp, wells, streams = self._mass_links
        for sp, t in self.transport.items():
            b = t.budget(rp, wells, streams)
            x = t.concentration()
            k = f"transport:{sp}:"
            L.enter(k + "recharge:sent", float(self.cfg.transport[sp]["recharge"]) * L.today("recharge:sent"))
            L.enter(k + "recharge:received", b.recharge)
            for w, m in b.wells.items():
                L.enter(k + f"well:{w}:withdrawn", m)
                L.enter(k + f"well:{w}:concentration", float(x[self.mf6.node(*self._well_cell(w))]))
            for c in b.to_stream:
                L.enter(k + f"stream:{c}:to_stream", b.to_stream[c])
                L.enter(k + f"stream:{c}:leakage_received", b.from_stream[c])
                L.enter(k + f"stream:{c}:leakage_sent", float(self.cfg.transport[sp]["stream"]) * b.leakage_m3[c])
            for p, m in b.boundary.items():
                L.enter(k + f"boundary:{p}", m)
            L.enter(k + "storage_release", b.storage_release)
            L.enter(k + "closure_share", b.closure_share)
        self._set_intakes()

    def _well_cell(self, name: str):
        return next(w.cell for w in self.cfg.wells if w.name == name)

    def _well_node(self, w: C.WellConfig):
        """A well's node (a WEL well's cell), or ("maw", its 0-based number, its model) for a multi-aquifer well."""
        if w.well is None:
            return self.mf6.node(*w.cell, model=w.model)
        v = self._view(w.model)
        if not v.has_package(w.package):
            raise ValueError(f"well {w.name}: the model has no package {w.package!r}")
        t = package_type(v, w.package)
        if t != "maw":
            raise ValueError(f"well {w.name}: well = {w.well} names a multi-aquifer well, but {w.package} is a "
                             f"{t.upper() or 'missing'} package")
        if not 1 <= w.well <= v.maw_wells(w.package):
            raise ValueError(f"well {w.name}: {w.package} has wells 1 to {v.maw_wells(w.package)}")
        if self.cfg.transport:
            raise ValueError(f"well {w.name}: transport with a multi-aquifer well isn't built")
        return ("maw", w.well - 1, v.model)

    def _views(self) -> list:
        """Each coupled groundwater model's view: several (``[modflow6] models``), or the one."""
        return self.mf6.views if isinstance(self.mf6, Models) else [self.mf6]

    def _view(self, model: str | None):
        """A link's model: the one named among several coupled, or the one."""
        return self.mf6.view_of(model) if isinstance(self.mf6, Models) else self.mf6

    # ---------- one day ----------
    def _state(self) -> R.State:
        q = self.sw.get(S.QOUT)
        heads = self.mf6.head() if self.mf6 is not None else None

        def head_at(c):
            if heads is None:
                raise ValueError("a drawdown rule needs a [modflow6] model")
            return float(heads[self.mf6.node(*c)])

        return R.State(day=self.day, flows=dict(zip(self.sw.names(S.QOUT), q)), head_at=head_at,
                       memory=self.rule_memory)

    def step(self):
        d = self.day
        L = self.ledger
        rt = self.router
        for wx, w in self.weather.items():  # a series that ends before the run would leave SWAT+ on its files
            if self.start + dt.timedelta(days=d) > w["to"]:
                raise ValueError(f"[weather] {wx}: the series ends {w['to']}; the run needs "
                                 f"{self.start + dt.timedelta(days=d)}")

        # 1. rules
        if self.rules:
            state = self._state()
            decisions = []
            for rule in self.rules:
                dec = rule(state)
                L.enter(f"rule:{rule.name}:active", 1.0 if dec else 0.0)
                decisions.append(dec)
            controls = R.combine(decisions)
            for u, wu in self.users:
                wu.set_controls(**controls.get(u.name, {}))

        # 2-3. users and routing
        for u, wu in self.users:
            use = wu.advance_day(d)
            k = f"user:{u.name}:"
            L.enter(k + "withdrawal", use.withdrawal)
            L.enter(k + "gain", use.gain)
            L.enter(k + "storage_less_gain", use.storage_change - use.gain)
            L.enter(k + "consumption", use.consumption)
            L.enter(k + "return", use.returned)
            for c, kg in use.constituents.items():
                L.enter(k + f"{c}:kg", kg)
            for sp in self.transport:  # the intake from the aquifer: its concentration and the mass it brought in
                if (u.name, sp) in self.intake:
                    L.enter(k + f"intake:{sp}", self.intake[(u.name, sp)])
                    L.enter(k + f"intake:{sp}:kg", use.withdrawal * self.intake[(u.name, sp)])
            if use.energy is not None:
                L.enter(k + "MWh", use.energy)
            for name, v in use.reported.items():
                L.enter(k + name, v)
            for route, vol, is_return in ((u.withdrawal, use.withdrawal, False), (u.returns, use.returned, True)):
                kd, name = kind(route)
                if kd == "utility":
                    self.utilities[name].add_user(0.0 if is_return else vol, vol if is_return else 0.0)
                else:
                    rt.request(route, vol)
                    if is_return and use.return_temperature is not None:  # heat, m3 C, for the stream temperature
                        L.enter(f"route:{route}:heat_m3C", vol * use.return_temperature)
        for chain in self.utilities.values():
            for route, vol in chain.close_day(L).items():
                rt.request(route, vol)
        if not self.direct_exchange:
            for route, vol in exchange_requests(self.due).items():
                rt.request(route, vol)

        # 4. the surface model (SWAT+ live, or the next day of a NextGen run)
        rt.apply_swat()
        gw_soil = getattr(self, "gw_soil", None)
        if gw_soil is not None:  # yesterday's groundwater, into the soil layers today
            L.enter("gw_soil:due", gw_soil.before_swat())
        if self.direct_exchange:
            names = self.sw.names(S.GW_EXCHANGE)
            self.sw.set(S.GW_EXCHANGE, np.array([self.due.get(c, 0.0) / DAY for c in names]))
        if self.lakes is not None:  # yesterday's lake leakage, on the reservoirs today
            self.lakes.before_swat()
        self.swat.update()
        rt.read_swat()
        if self.lakes is not None:
            self.lakes.after_swat()
        if self.irrigation is not None:  # the day's irrigation, pumped from MODFLOW 6's wells today
            self._irrigate(L, rt)
        if gw_soil is not None:
            L.enter("gw_soil:received", gw_soil.received())
        if self.direct_exchange:
            got = dict(zip(self.sw.names(S.GW_EXCHANGE), self.sw.get(S.GW_EXCHANGE) * DAY))
            self._exchanged = {c: got[c] for c in self.due}

        # 5. MODFLOW 6
        if self.mf6 is not None:
            self._modflow_day(L, rt)
        rt.close_day()

        # 6. outputs and checks
        q = self.sw.get(S.QOUT) if self._ch_index or self._nexus else None  # a project of HRUs only has no channels
        for c, i in self._ch_index.items():
            self.series[f"flow:{c}"].append(float(q[i]))
        for nx, idx in self._nexus.items():
            self.series[f"flow:{nx}"].append(float(q[idx].sum()))
        if self.out_heads:
            h = self.mf6.head()
            for n, node in self.out_heads:
                self.series[f"head:{n}"].append(float(h[node]))
        L.close_day()
        self.day += 1

    def _modflow_day(self, L, rt):
        self.mf6.prepare_day()
        if self.cfg.exchange and self.exchange is None:
            self._build_exchange()
        if not hasattr(self, "_gw_packages"):
            self._start_takeover()
        if self.transport:
            self._transport_day()
        sent = self.recharge.transfer()
        L.enter("recharge:sent", sent)
        L.enter("recharge:made", self.recharge.made)
        L.enter("recharge:by_the_map", sent - self.recharge.made)
        if self.aquifer_check is not None:
            L.enter("recharge:swat_aquifers", float(np.sum(self.sw.get(RC.RECHARGE) * DAY * self.aquifer_check)))
        rt.apply_wells()
        if self.lakes is not None:
            self.lakes.before_solve()
        if self.cfg.rivers:
            self._rivers_day()
        limit = None
        sfr = self.exchange is not None and self.exchange.kind == "sfr"
        if self.exchange is not None and (sfr or self.cfg.exchange.get("limit", True)):
            # what each channel has without today's exchange, which won't be there tomorrow when this one is taken:
            # SFR's inflow, or RIV's limit (a losing reach takes no more than its channel carries, exchange.limit)
            flow = dict(zip(self.sw.names(S.QOUT), self.sw.get(S.QOUT) * DAY))
            if self.direct_exchange:
                applied = self._exchanged
            else:
                applied = {c: L.today(f"route:swat:point_source:gw_{c}:delivered")
                           - L.today(f"route:swat:place_of_use:seep_{c}:delivered") for c in self.due}
            available = {c: max(flow[c] - applied[c], 0.0) for c in self.due}
            if sfr:
                depth = None
                if self.exchange._stage:  # B: the stage from SWAT+'s depth for the day just run
                    depth = dict(zip(self.sw.names(S.DEPTH), self.sw.get(S.DEPTH)))
                    for c in self.exchange._stage:
                        L.enter(f"exchange:{c}:sfr_depth_m", float(depth[c]))
                self.exchange.set_inflow(available, depth)
                for c, v in available.items():
                    L.enter(f"exchange:{c}:sfr_inflow", v)
            else:
                limit = self.exchange.limit(available)
        if self.gw_soil is not None:
            if self.gw_soil.mf6 is None:
                self.gw_soil.attach(self.mf6)
            self.gw_soil.before_solve()
        self.mf6.solve_day(limit)
        if limit is not None:
            for c, d in self.exchange.lowered.items():
                L.enter(f"exchange:{c}:stage_lowered_m", d)
            self.exchange.end_day()
        L.enter("recharge:received", self.recharge.received())
        if isinstance(self.mf6, Models):  # each coupled model's share
            for v in self.mf6.views:
                if v.has_package(self.cfg.recharge_package):
                    L.enter(f"recharge:received:{v.model.lower()}", float(v.simvals(self.cfg.recharge_package).sum()))
        if self.lakes is not None:
            self._lakes_day(L)
        if self.gw_soil is not None:
            sent, at_bottom = self.gw_soil.after_solve()
            L.enter("gw_soil:sent", sent)                  # out of MODFLOW 6 today, into SWAT+'s soil tomorrow
            L.enter("gw_soil:to_bottom_layer", at_bottom)  # of it, from cells whose water table is below the HRU's profile
            if self.gw_soil.rule == "drain" and self.gw_soil.cap == "room":
                L.enter("gw_soil:drains_raised", float(self.gw_soil.raised))
        rt.read_wells()
        if self.cfg.rivers:
            self._rivers_out()
        if self.cfg.exchange:
            self._enter_exchange(L)
            L.enter("swatplus:aquifer_outflow", float(self.sw.get(S.AQU_OUT).sum()) * DAY)
        self._groundwater_accounts(L)
        if self.transport:
            self._transport_accounts(L)

    def _build_exchange(self):
        """The stream-aquifer exchange, once MODFLOW 6 has read its first day (the packages have their entries)."""
        ex = self.cfg.exchange
        if self._exchange_models():
            return self._build_exchanges()
        mf = self._view(ex.get("model"))  # the model whose package the channels exchange through
        names = self.sw.names(S.QOUT)
        cmap = ChannelMap.from_spec(ex.get("map", ""), self.cfg.folder, mf, ex.get("package", "riv"),
                                    names, channels=list(self.due))
        if cmap.channels != list(self.due):
            raise ValueError(f"exchange map channels {cmap.channels} differ from {list(self.due)}")
        drains = None
        if ex.get("drains"):
            dp = ex["drains"].get("package", "drn")
            drains = (dp, ChannelMap.from_spec(ex["drains"]["map"], self.cfg.folder, mf, dp, names,
                                               channels=list(self.due)))
        if "sfr" in ex and package_type(mf, ex.get("package", "riv")) != "sfr":
            raise ValueError("[exchange] sfr chooses how an SFR package's channels run; this exchange is on a "
                             f"{package_type(mf, ex.get('package', 'riv')).upper()} package")
        self.exchange = StreamExchange(mf, ex.get("package", "riv"), cmap, drains, sfr=ex.get("sfr", "route"))
        if any(m == "stage" for m, _ in self.exchange.modes.values()) and not self.sw.has(S.DEPTH):
            raise ValueError("SFR channels at a fixed stage take SWAT+'s channel depth (swatplus-bmi 0.2.0)")
        if self.exchange.kind == "sfr" and ex.get("limit") is False:
            raise ValueError("[exchange] limit = false: SFR limits a reach's loss to what it carries itself; the "
                             "limit is for a river (RIV) package")

    def _exchange_models(self) -> list[str]:
        """With several groundwater models coupled, the models the exchange map's rows name (a row without one is the
        first model's), if more than one: the streams then run through each model's package (``StreamExchanges``)."""
        import csv
        ex = self.cfg.exchange
        spec = ex.get("map", "")
        if not isinstance(self.mf6, Models) or not spec or spec.startswith("all:"):
            return []
        with open(self.cfg.folder / spec, newline="") as f:
            named = {((r.get("model") or "").strip() or self.mf6.model).upper() for r in csv.DictReader(f)}
        models = [m for m in self.mf6.names if m in named]
        if named - set(models):
            raise ValueError(f"exchange map: models {sorted(m.lower() for m in named - set(models))} aren't coupled")
        return models if len(models) > 1 else []

    def _build_exchanges(self):
        """The stream exchange through each coupled model's package, its channels mapped there; the mover between
        the models lined up and cut as a connection between channels is (exchange.sfr_problems, ``others``)."""
        ex = self.cfg.exchange
        if ex.get("drains"):
            raise ValueError("[exchange] drains with streams in several groundwater models isn't built")
        pkg, names = ex.get("package", "riv"), self.sw.names(S.QOUT)
        views = [self.mf6.view_of(m) for m in self._exchange_models()]
        maps = [ChannelMap.from_spec(ex["map"], self.cfg.folder, v, pkg, names, channels=list(self.due),
                                     first_model=self.mf6.model) for v in views]
        mapped = [c for m in maps for c in m.channels]
        if sorted(mapped) != sorted(self.due) or len(set(mapped)) != len(mapped):
            twice = sorted({c for c in mapped if mapped.count(c) > 1})
            raise ValueError(f"exchange map: channels {twice} have reaches in more than one model; a channel's "
                             "reaches are in one model" if twice else
                             f"exchange map channels {mapped} differ from {list(self.due)}")
        parts = []
        for v, cmap in zip(views, maps):
            others = {(w.model.lower(), pkg.lower()): m for w, m in zip(views, maps) if w is not v}
            parts.append(StreamExchange(v, pkg, cmap, None, sfr=ex.get("sfr", "route"), others=others))
        self.exchange = StreamExchanges(parts)
        if any(m == "stage" for m, _ in self.exchange.modes.values()) and not self.sw.has(S.DEPTH):
            raise ValueError("SFR channels at a fixed stage take SWAT+'s channel depth (swatplus-bmi 0.2.0)")

    def _coupled_entries(self, mf) -> dict[str, tuple[np.ndarray, list[str]]]:
        """Each package with flows (SIMVALS) of one coupled groundwater model (``mf``, its view), by upper-case name,
        with a mask of the entries a link sets and the links that set them. Entries no link sets are the model's
        (coupling/existing.py)."""
        pkgs = EX.flow_packages(mf)
        ex = self.exchange
        several = isinstance(self.mf6, Models)
        out = {}
        for p in pkgs:
            n = mf.nbound(p)
            m = np.zeros(n, bool)
            links: list[str] = []

            def take(sel, what, m=m, links=links):
                sel = np.asarray(sel, bool)
                if sel.any():
                    m[sel] = True
                    links.append(what)

            if p == self.cfg.recharge_package.upper() and self.recharge is not None:
                take(np.ones(n, bool), "recharge from SWAT+")
            for part in (getattr(ex, "parts", [ex]) if ex is not None else []):
                if part.mf6 is mf and p == part.package.upper():
                    modes = [m for m, _ in part.modes.values()]
                    how = (f" (routed: {modes.count('route')} channels, at SWAT+'s stage: {modes.count('stage')})"
                           if "stage" in modes else "")
                    take(part.map.entry_channel[:n] >= 0, "stream exchange with SWAT+'s channels" + how)
                if part.mf6 is mf and part.drains and p == part.drains[0].upper():
                    take(part.drains[1].entry_channel[:n] >= 0, "seepage to SWAT+'s channels")
            if self.gw_soil is not None and p == self.gw_soil.package.upper():
                take(np.ones(n, bool), "groundwater up into SWAT+'s soil")
            for lk, s, v in zip(self.cfg.lakes, self.lakes.conn if self.lakes is not None else [],
                                self.lakes.views if self.lakes is not None else []):
                if v is mf and p == lk.package.upper():
                    sel = np.zeros(n, bool)
                    sel[s] = True
                    how = "active within the day (A)" if lk.form == "active" else "at a fixed stage (B)"
                    take(sel, f"lake {lk.lake}'s exchange with SWAT+'s reservoir {lk.reservoir}, {how}")
            wells = np.zeros(n, bool)
            named = []
            for name, (pkg, node) in self.router.wells.items():
                if pkg.upper() != p:
                    continue
                if isinstance(node, tuple):  # a multi-aquifer well: its connections
                    if self._view(node[2]) is mf:
                        wells[mf.maw_connections(p, node[1])] = True
                        named.append(name)
                    continue
                v, local = self.mf6.model_of(node) if several else (self.mf6, node)
                if v is not mf:
                    continue
                try:
                    wells[mf.entry(p, local)] = True
                except ValueError:  # listed under the well by the start-up; refused by the router in a run
                    continue
                named.append(name)
            take(wells, "users' wells: " + ", ".join(named))
            out[p] = (m, links)
        return out

    def _start_takeover(self):
        """What Basin Ledger takes over in MODFLOW 6 and what it leaves as the model's (coupling/existing.py), once the
        packages have their entries: ``self.takeover``, one row per package (``<model>/<package>`` with several
        models coupled). What would count water twice stops the start-up; what a reviewer should know is a warning
        (a check's notes)."""
        views = self._views()
        several = len(views) > 1
        coupled_names = {v.model for v in views}
        self._gw_packages = []  # (view, {package: mask of the entries a link sets})
        self.takeover = []
        notes, bad = [], []
        for v in views:
            coupled = self._coupled_entries(v)
            self._gw_packages.append((v, {p: m for p, (m, _) in coupled.items()}))
            rows = EX.roles(v, coupled)
            if several:
                for r in rows:
                    r.package = f"{v.model.lower()}/{r.package}"
            self.takeover += rows
        bad, notes = EX.problems(self.takeover, self.gw_soil is not None)
        # other groundwater models of the simulation: what crosses to one not coupled is water from or to outside;
        # between coupled models it stays inside the groundwater
        for v in views:
            for name, other, sign in v.exchanges():
                if other.upper() not in coupled_names:
                    notes.append(f"{name.lower()}: groundwater to and from {other}, another model of the simulation "
                                 "that Basin Ledger doesn't couple, booked as from or to outside")
                elif sign > 0:
                    notes.append(f"{name.lower()}: groundwater between {v.model.lower()} and {other}, both coupled: "
                                 f"inside the groundwater (groundwater:exchange:{name.lower()}, positive into "
                                 f"{v.model.lower()})")
        if self.exchange is not None and any(m == "stage" for m, _ in self.exchange.modes.values()):
            routed = [c for c, (m, _) in self.exchange.modes.items() if m == "route"]
            notes.append(f"{self.exchange.package}: channels routed by SFR (A): {', '.join(routed) or 'none'}")
            for c, (m, why) in self.exchange.modes.items():
                if m == "stage":
                    notes.append(f"{self.exchange.package}: channel {c} at a fixed stage from SWAT+'s depth (B): "
                                 + "; ".join(w.split(": ", 1)[-1] for w in why))
        if self.cfg.rivers:
            notes.append(f"{self.cfg.rivers.get('package', 'riv')}: [rivers] sets its stage from the surface model's "
                         "depth; its flow stays the model's")
        for n in notes:
            warnings.warn(n, EX.TakeoverNote, stacklevel=2)
        if bad:
            raise ValueError("MODFLOW 6's own packages would count water twice:\n- " + "\n- ".join(bad))

    def _groundwater_accounts(self, L):
        """Every flow in or out of the groundwater that no link carries, by package (water from or to outside the
        coupled models: fixed or general heads, drains and rivers not coupled, recharge packages Basin Ledger
        doesn't set, and what crosses to a groundwater model of the simulation that isn't coupled, by its GWF-GWF
        exchange), the release from storage, and the closure, all in m3 for the day. Positive is into the aquifer.
        With several models coupled, each model's accounts are under ``<model>/``, a GWF-GWF exchange between two of
        them is booked once, inside the groundwater (``groundwater:exchange:<name>``), and each model's closure is
        checked as well as the whole's. The closure is checked every day, against the day's gross flow."""
        several = len(self._gw_packages) > 1
        coupled_names = {v.model for v, _ in self._gw_packages}
        if not getattr(self, "_closure_declared", False):
            L.zero("groundwater:closure_share", atol=1e-6)
            if several:
                for v, _ in self._gw_packages:
                    L.zero(f"groundwater:{v.model.lower()}:closure_share", atol=1e-6)
            self._closure_declared = True
        net_all = gross_all = sto_all = 0.0
        for mf, pkgs in self._gw_packages:
            pre = f"{mf.model.lower()}/" if several else ""
            net = gross = 0.0
            for p, m in pkgs.items():
                q = mf.simvals(p)
                L.enter(f"groundwater:boundary:{pre}{p.lower()}", float(q[~m[: len(q)]].sum()))
                net += float(q.sum())
                gross += float(np.abs(q).sum())
            for name, other, sign in mf.exchanges():
                q = float(mf.exchange_flow(name, sign).sum())
                if other.upper() not in coupled_names:
                    L.enter(f"groundwater:boundary:{pre}{name.lower()}", q)
                elif sign > 0:
                    L.enter(f"groundwater:exchange:{name.lower()}", q)
                net += q
                gross += abs(q)
            sto = 0.0
            for v in ("STRGSS", "STRGSY"):
                try:
                    x = mf.get(v, mf.model, "STO")
                except Exception:  # noqa: BLE001, S112  (a steady model has no storage)
                    continue
                sto += float(np.sum(x))
            if several:
                L.enter(f"groundwater:{mf.model.lower()}:storage_release", sto)
                L.enter(f"groundwater:{mf.model.lower()}:closure_share", (net + sto) / max(gross + abs(sto), 1.0))
            net_all, gross_all, sto_all = net_all + net, gross_all + gross, sto_all + sto
        L.enter("groundwater:storage_release", sto_all)  # MODFLOW 6's STRGSS + STRGSY: positive is released into the flow
        # against the day's gross flow, with a floor of 1 m3 so a day with almost no flow isn't noise over noise
        L.enter("groundwater:closure_share", (net_all + sto_all) / max(gross_all + abs(sto_all), 1.0))

    def _enter_exchange(self, L):
        """Today's SWAT+ side of yesterday's exchange, then today's MODFLOW 6 exchange as tomorrow's due."""
        for c in self.due:
            k = f"exchange:{c}:"
            if self.direct_exchange:
                applied = self._exchanged[c]                 # what SWAT+ actually exchanged
                short = max(applied - self.due[c], 0.0) if self.due[c] < 0 else 0.0  # a loss it couldn't supply
            else:
                applied = (L.today(f"route:swat:point_source:gw_{c}:delivered")
                           - L.today(f"route:swat:place_of_use:seep_{c}:delivered"))
                short = L.today(f"route:swat:place_of_use:seep_{c}:shortfall")
            L.enter(k + "due", self.due[c])
            L.enter(k + "applied", applied)
            L.enter(k + "shortfall", -short)
        self.due = self.exchange.after_solve()
        for c, net in self.due.items():
            L.enter(f"exchange:{c}:to_stream", net)
            if self.exchange.drains is not None:
                L.enter(f"exchange:{c}:seepage", self.exchange.seepage[c])  # the drains' part of to_stream

    def run(self, days: int | None = None):
        n = self.days - self.day if days is None else days
        for _ in range(n):
            self.step()
        return self

    # ---------- restarts ----------
    def save(self, folder: str | Path):
        """The state at the pause after the last day run, for ``Basin(..., restart=folder)``."""
        folder = Path(folder)
        if folder.exists():
            shutil.rmtree(folder)
        folder.mkdir(parents=True)
        self.swat.save_state(str(folder / "surface"))
        if self.mf6 is not None:
            np.save(folder / "heads.npy", self.mf6.head())
            for sp, t in self.transport.items():
                np.save(folder / f"concentration_{sp}.npy", t.concentration())
        (folder / "basin.json").write_text(json.dumps({
            "day": self.day, "ledger": self.ledger.state(), "series": self.series,
            "rule_memory": self.rule_memory,
            "exchange_due": self.due,
            "lakes_due": self.lakes.due if self.lakes is not None else None,
            "gw_soil_due": self.gw_soil.due.tolist() if getattr(self, "gw_soil", None) is not None else None,
        }))

    # ---------- outputs ----------
    def dates(self) -> list[str]:
        return [(self.start + dt.timedelta(days=i)).isoformat() for i in range(self.ledger.days)]

    def write_outputs(self, folder: str | Path | None = None) -> Path:
        """ledger.csv (every account, m3 per day unless named otherwise; ``transport:`` accounts in kg, and
        concentrations in kg m-3), series.csv (flows m3 s-1, heads m) and
        summary.json (status, totals, what's in transit)."""
        out = Path(folder) if folder else self.folder / "outputs"
        out.mkdir(parents=True, exist_ok=True)
        dates = self.dates()
        self.ledger.write_csv(out / "ledger.csv", dates)
        with open(out / "series.csv", "w") as f:
            cols = sorted(self.series)
            f.write(",".join(["day", "date", *cols]) + "\n")
            for i in range(len(dates)):
                f.write(",".join([str(i), dates[i], *(repr(self.series[c][i]) for c in cols)]) + "\n")
        from .calibration import model_statuses, record_targets
        targets = record_targets(self.cfg.calibration_record)
        if self.status == "calibrated":
            note = ("calibrated: " + "; ".join(f"{m} to {t}" for m, t in targets.items()) +
                    ". Water users' settings are assumptions or applicant data, not calibrated.") if targets else None
        else:
            note = ("screening run: not every model is calibrated; results show how the engine and scenarios work, "
                    "not what would happen in this basin")
        summary = {
            "example": self.cfg.name, "status": self.status,
            "models": model_statuses(self.cfg.calibration_record, self.cfg.models),
            "calibrated_to": targets,
            "status_note": note,
            "start": self.start.isoformat(), "days": self.ledger.days,
            "totals": {k: v for k, v in self.ledger.totals().items() if not math.isnan(v)},
            "in_transit": self.due,
        }
        if self.weather:
            summary["weather"] = {k: {"from": w["from"].isoformat(), "to": w["to"].isoformat(), "gages": w["series"].gages}
                                  for k, w in self.weather.items()}
        if getattr(self, "gw_soil", None) is not None:
            summary["in_transit_gw_soil"] = self.gw_soil.due_volume()
        (out / "summary.json").write_text(json.dumps(summary, indent=2))
        return out

    def finalize(self):
        for _, wu in self.users:
            wu.finalize()
        if self.mf6 is not None:
            self.mf6.finalize()
        self.swat.finalize()

    @classmethod
    def check(cls, cfg: C.ExampleConfig, users: bool = True, run: str = "check",
              takeover: list | None = None, units: list | None = None) -> tuple[list[str], list[str]]:
        """Every problem the example's start-up can find, without running a day: (problems, notes). The models start
        as for a run, MODFLOW 6 reads its first day (the packages have their entries then) and the links are built:
        the SWAT+ settings, the parameters, the recharge map and its areas, the wells', head outputs' and rules'
        cells, the routes, users and rules, the exchange and its maps (SFR's checks included), the rivers,
        groundwater into the soil and transport. A part with a problem is listed and left out, and the rest go on;
        a problem that leaves a model unable to start ends the check there. Notes are the start-up's warnings.
        ``takeover``, if given, gets what Basin Ledger takes over in MODFLOW 6 and what it leaves, one
        ``existing.PackageRole`` per package (coupling/existing.py); ``units``, if given, MODFLOW 6's units and
        the factors converting them at every link (``Modflow6.units_summary``)."""
        problems: list[str] = []
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                b = cls(cfg, run=run, users=users, days=1, _problems=problems)
            except Exception as e:  # noqa: BLE001
                problems.append(f"the start-up stopped: {e}")
                b = None
            if b is not None:
                try:
                    if b.mf6 is not None and b._part("MODFLOW 6's first day", lambda: b.mf6.prepare_day() or True):
                        if cfg.exchange:
                            b._part("[exchange]", b._build_exchange)
                        if cfg.rivers:
                            b._part("[rivers]", b._build_rivers)
                        for name, (pkg, node) in b.router.wells.items():
                            if not isinstance(node, tuple):  # a multi-aquifer well is checked by its number
                                b._part(f"well {name}", lambda pkg=pkg, node=node: b.mf6.entry(pkg, node))
                        if b.gw_soil is not None:
                            b._part("[groundwater_to_soil]", lambda: b.gw_soil.attach(b.mf6))
                        if b.transport:
                            b._part("[transport]", b._build_mass_links)
                        b._part("[modflow6]", b._start_takeover)
                        if takeover is not None:
                            takeover.extend(getattr(b, "takeover", []))
                    if b.mf6 is not None and units is not None:
                        units.append(b.mf6.units_summary())
                finally:
                    b.finalize()
        return problems, list(dict.fromkeys(str(w.message) for w in caught))

    @classmethod
    def spin_up(cls, cfg: C.ExampleConfig, days: int = 365, cycles: int = 5, tol: float = 0.01, run: str = "spinup",
                users: bool = False, prepare=None, heads: np.ndarray | None = None) -> tuple[np.ndarray, list[dict]]:
        """MODFLOW 6's starting heads under the coupling: the first ``days`` run again and again, SWAT+ from its own
        start each time and MODFLOW 6 from the heads the last cycle ended with, until the cycle's net release from
        storage is within ``tol`` of the recharge it received (or ``cycles`` run). A start that isn't consistent with
        the coupling drains or fills storage over the run (a study basin's first coupled run drained 49 million m3,
        15% of its recharge). Returns the heads and each cycle's storage release and recharge, m3."""
        report = []
        for i in range(cycles):
            b = cls(cfg, run=f"{run}-{i}", users=users, days=days, prepare=prepare, heads=heads)
            try:
                b.run()
                t = b.ledger.totals()
                heads = b.mf6.head()
            finally:
                b.finalize()
            rel, rch = t.get("groundwater:storage_release", 0.0), t.get("recharge:received", 0.0)
            report.append({"cycle": i, "storage_release_m3": rel, "recharge_m3": rch,
                           "share": rel / rch if rch else float("nan")})
            if rch and abs(rel / rch) < tol:
                break
        return heads, report

    @classmethod
    def from_config(cls, path, run: str = "base", users: bool = True, days: int | None = None,
                    restart: str | Path | None = None, prepare=None) -> "Basin":
        """``users=False`` runs the example without its water users (the baseline); ``days`` shortens the run;
        ``restart`` carries on from a folder written by ``save``."""
        return cls(C.load(path), run, users, days, restart, prepare)
