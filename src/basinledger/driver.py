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
from pathlib import Path

import numpy as np
from swatplus_bmi import SwatPlusBmi

from . import config as C
from . import rules as R
from .chain import UtilityChain
from .coupling import swatplus as S
from .coupling.exchange import ChannelMap, StreamExchange, requests as exchange_requests
from .coupling.modflow6 import Modflow6, remaining_simulation
from .coupling.recharge import RechargeLink
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
                 wrap_surface=None):
        """``prepare(model, folder)``, if given, is called on the run's copy of each model's files ("swatplus",
        "modflow6") before that model starts: where a calibration writes a parameter set. ``heads``: MODFLOW 6's
        starting heads (per reduced node), in place of its input's, as from ``spin_up``. ``wrap_surface(model)``,
        if given, returns what the driver steps in place of the surface model once it has started: a study's way to
        hold or replace one of its values through BMI each day (a process substitution), with everything else
        passed through."""
        self.cfg = cfg
        self.status = cfg.calibration_status  # 'screening' unless the example's record says 'calibrated'
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
            if prepare:
                prepare("swatplus", swat_dir)
            problems = S.coupling_problems(swat_dir, two_way=bool(ex))
            if problems:
                raise ValueError("SWAT+ isn't set up to be coupled:\n- " + "\n- ".join(problems))
            self.swat = SwatPlusBmi()
            self.swat.initialize(str(swat_dir))
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
            elif saved:
                remaining_simulation(cfg.mf6_workspace, ws, self.day)
            else:
                shutil.copytree(cfg.mf6_workspace, ws)
            if prepare:
                prepare("modflow6", ws)
            self.mf6 = Modflow6(ws, cfg.mf6_model)
            if saved:
                self.mf6.set_heads(np.load(Path(restart) / "heads.npy"))
            elif heads is not None:
                self.mf6.set_heads(np.asarray(heads, float))
            if self.mf6.end < self.days - self.day - 1e-9:
                raise ValueError(f"MODFLOW 6 simulation is {self.mf6.end} days; the run needs {self.days - self.day}")
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
            raise ValueError("[rivers] stage must be fixed or surface_depth")
        if cfg.rivers and cfg.rivers.get("stage") in ("surface_depth", "nextgen_depth") and not self.sw.has(S.DEPTH):
            raise ValueError("[rivers] stage = surface_depth needs the surface model's channel depth "
                             "(NextGen's t-route, or swatplus-bmi 0.2.0)")

        self.ledger = Ledger()
        L = self.ledger
        # recharge: SWAT+ aquifers -> MODFLOW 6 cells
        if self.mf6 is not None:
            aquifers = self.sw.names("soil_water_table__recharge_rate")
            self.map = self._area_map(cfg.recharge_map, aquifers)
            areas = S.aquifer_areas(swat_dir) if swat_dir else self.swat.areas
            source = np.array([areas[a] for a in aquifers])
            self.map.check(self.mf6.cell_area, source, cfg.area_rtol)
            self.recharge = RechargeLink(self.swat, self.mf6, self.map, cfg.recharge_package, source)
            L.pair("recharge:sent", "recharge:received")
            # what the surface model made, and what the map added to it (+) or couldn't place (-): water from outside
            L.pair("recharge:sent", ["recharge:made", "recharge:by_the_map"])
        wells = {w.name: (w.package, self.mf6.node(*w.cell)) for w in cfg.wells}
        self.router = Router(self.sw, self.mf6, wells, L)

        # stream-aquifer exchange (built on the first day, once MODFLOW 6's river package has its entries)
        self.exchange = None
        self.due = {c: 0.0 for c in ex_channels}  # yesterday's exchange, m3 per channel, applied in SWAT+ today
        if saved and saved.get("exchange_due"):
            self.due = dict(saved["exchange_due"])
        if ex:
            for c in ex_channels:
                if not self.direct_exchange:
                    self.router.check(f"swat:point_source:gw_{c}")
                    self.router.check(f"swat:place_of_use:seep_{c}")
                k = f"exchange:{c}:"
                L.pair(k + "due", [k + "applied", k + "shortfall"], rtol=FLOAT32)
            L.zero("swatplus:aquifer_outflow", atol=1e-9)

        # utilities, users, rules
        self.utilities = {u["name"]: UtilityChain.from_config(u) for u in cfg.utilities}
        for u in self.utilities.values():
            u.declare(L)
            for s in u.sources:
                self.router.check(s.route)
            self.router.check(u.outfall)
        self.users = []
        if users:
            for u in cfg.users:
                comp = load_component(u.component)
                wu = WaterUser(u.name, comp, str(u.config) if u.config else None,
                               variables=u.variables or None, constituents=u.constituents,
                               report=u.report or None)
                wu.inputs = self._forcing(u, comp)  # after initialize, when a component's inputs are known
                self.users.append((u, wu))
                for route in (u.withdrawal, u.returns):
                    k, name = kind(route)
                    if k == "utility":
                        if name not in self.utilities:
                            raise KeyError(f"user {u.name}: no utility named {name!r}")
                    else:
                        self.router.check(route)
                k = f"user:{u.name}:"
                L.pair(k + "withdrawal", [k + "consumption", k + "return", k + "storage_less_gain"])
        # intakes from the aquifer: each user on a configured well, for each species its component takes
        self.intakes = []  # (user, WaterUser, input variable, species, node)
        for u, wu in self.users:
            k, name = kind(u.withdrawal)
            if k != "well":
                continue
            takes = set(wu.bmi.get_input_var_names())
            for sp in self.transport:
                var = f"water_user_intake_water_{sp}__mass_concentration"
                if var not in takes:
                    continue
                if var in (u.inputs or {}):
                    raise ValueError(f"user {u.name}: {var} is set from the aquifer ([transport.{sp}]); "
                                     "remove its constant input")
                self.intakes.append((u.name, wu, var, sp, wells[name][1]))
        self._set_intakes()
        users_by_name = {u.name: wu for u, wu in self.users}
        self.rules = []
        for r in cfg.rules:
            settings = {k: v for k, v in r.items() if k not in ("name", "kind")}
            fn = R.BUILT_IN.get(r["kind"]) or _example_function(cfg, r["kind"], "rule")
            if users and settings.get("user") not in users_by_name:
                raise KeyError(f"rule {r.get('name', r['kind'])}: no user {settings.get('user')!r}")
            need = {"curtailment": "limit", "drawdown_limit": "limit", "dry_mode": "dry_mode"}.get(r["kind"])
            if users and need and need not in users_by_name[settings["user"]].control_vars:
                raise KeyError(f"rule {r.get('name')}: user {settings['user']!r} takes no {need!r} input")
            self.rules.append(R.Rule(r.get("name", r["kind"]), fn, settings))
        self.rule_memory = {}

        # outputs
        out = cfg.outputs
        self.out_channels = out.get("channels", channels)
        self._ch_index = {c: channels.index(c) for c in self.out_channels}
        self.out_heads = [(h["name"], self.mf6.node(*h["cell"])) for h in out.get("heads", [])]
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
        pkg = rv.get("package", "riv")
        if self.rivers is None:
            cmap = ChannelMap.from_spec(rv.get("map", ""), self.cfg.folder, self.mf6, pkg, self.sw.names(S.QOUT))
            self.rivers = (pkg, cmap, self.mf6.get("RBOT", self.mf6.model, pkg.upper()).astype(float),
                           self.mf6.get("STAGE", self.mf6.model, pkg.upper()).astype(float))
        pkg, cmap, rbot, stage0 = self.rivers
        if rv.get("stage") in ("surface_depth", "nextgen_depth"):
            depth = dict(zip(self.sw.names(S.DEPTH), self.sw.get(S.DEPTH)))  # yesterday's (SWAT+) or today's (replay)
            e = cmap.entry_channel
            stage = stage0.copy()
            for i, c in enumerate(cmap.channels):
                stage[e == i] = rbot[e == i] + float(rv.get("bed_thickness", 1.0)) + depth[c]
            self.mf6.set_stage(pkg, stage)

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
                series = swatplus_inputs(self.folder / "swatplus", u.forcing["station"], start)
            else:
                raise NotImplementedError(f"user {u.name}: forcing from {src!r} with a {self.cfg.surface} surface "
                                          "model (nextgen with NextGen, swatplus with SWAT+; AORC is planned)")
            out.update({k: v for k, v in series.items() if k in takes})
        for k, v in (u.inputs or {}).items():
            if k not in takes:
                raise KeyError(f"user {u.name}: {k!r} is not an input of {u.component}")
            out[k] = (lambda value: (lambda t: value))(float(v))
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

    def _transport_day(self):
        """The day's inflow concentrations, after MODFLOW 6 prepares the day [A]. The couplings are found on the first
        day, once the packages have their entries."""
        if not hasattr(self, "_mass_links"):
            rp = self.cfg.recharge_package
            wells = {w.name: (w.package, self.mf6.entry(w.package, self.mf6.node(*w.cell))) for w in self.cfg.wells}
            streams = []
            if self.exchange is not None:
                streams.append((self.exchange.package, self.exchange.map.entry_channel, self.exchange.map.channels))
                if self.exchange.drains is not None:
                    dp, dmap = self.exchange.drains
                    streams.append((dp, dmap.entry_channel, dmap.channels))
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
        if self.direct_exchange:
            names = self.sw.names(S.GW_EXCHANGE)
            self.sw.set(S.GW_EXCHANGE, np.array([self.due.get(c, 0.0) / DAY for c in names]))
        self.swat.update()
        rt.read_swat()
        if self.direct_exchange:
            got = dict(zip(self.sw.names(S.GW_EXCHANGE), self.sw.get(S.GW_EXCHANGE) * DAY))
            self._exchanged = {c: got[c] for c in self.due}

        # 5. MODFLOW 6
        if self.mf6 is not None:
            self._modflow_day(L, rt)
        rt.close_day()

        # 6. outputs and checks
        q = self.sw.get(S.QOUT)
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
            ex = self.cfg.exchange
            names = self.sw.names(S.QOUT)
            cmap = ChannelMap.from_spec(ex.get("map", ""), self.cfg.folder, self.mf6, ex.get("package", "riv"),
                                        names, channels=list(self.due))
            if cmap.channels != list(self.due):
                raise ValueError(f"exchange map channels {cmap.channels} differ from {list(self.due)}")
            drains = None
            if ex.get("drains"):
                dp = ex["drains"].get("package", "drn")
                drains = (dp, ChannelMap.from_spec(ex["drains"]["map"], self.cfg.folder, self.mf6, dp, names,
                                                   channels=list(self.due)))
            self.exchange = StreamExchange(self.mf6, ex.get("package", "riv"), cmap, drains)
        if self.transport:
            self._transport_day()
        sent = self.recharge.transfer()
        L.enter("recharge:sent", sent)
        L.enter("recharge:made", self.recharge.made)
        L.enter("recharge:by_the_map", sent - self.recharge.made)
        rt.apply_wells()
        if self.cfg.rivers:
            self._rivers_day()
        limit = None
        if self.exchange is not None and self.cfg.exchange.get("limit", True):
            # a losing reach takes no more than its channel carries today, within the solve (exchange.limit)
            # what each channel has without today's exchange, which won't be there tomorrow when this one is taken
            flow = dict(zip(self.sw.names(S.QOUT), self.sw.get(S.QOUT) * DAY))
            if self.direct_exchange:
                applied = self._exchanged
            else:
                applied = {c: L.today(f"route:swat:point_source:gw_{c}:delivered")
                           - L.today(f"route:swat:place_of_use:seep_{c}:delivered") for c in self.due}
            limit = self.exchange.limit({c: max(flow[c] - applied[c], 0.0) for c in self.due})
        self.mf6.solve_day(limit)
        if limit is not None:
            for c, d in self.exchange.lowered.items():
                L.enter(f"exchange:{c}:stage_lowered_m", d)
            self.exchange.end_day()
        L.enter("recharge:received", self.recharge.received())
        rt.read_wells()
        if self.cfg.rivers:
            self._rivers_out()
        if self.cfg.exchange:
            self._enter_exchange(L)
            L.enter("swatplus:aquifer_outflow", float(self.sw.get(S.AQU_OUT).sum()) * DAY)
        self._groundwater_accounts(L)
        if self.transport:
            self._transport_accounts(L)

    def _groundwater_accounts(self, L):
        """Every flow in or out of the groundwater model that no link carries, by package (water from or to outside
        the coupled models: fixed or general heads, drains and rivers not coupled, recharge packages Basin Ledger
        doesn't set), the release from storage, and the model's closure, all in m3 for the day. Positive is into
        the aquifer. The closure is checked every day, against the day's gross flow."""
        mf = self.mf6
        if not hasattr(self, "_gw_packages"):
            prefix = mf.model.upper() + "/"
            names = set(mf.api.get_input_var_names()) | set(mf.api.get_output_var_names())
            pkgs = sorted({n.split("/")[1] for n in names if n.startswith(prefix) and n.endswith("/SIMVALS")})
            coupled = {}
            for p in pkgs:
                n = mf.nbound(p)
                m = np.zeros(n, bool)
                if p == self.cfg.recharge_package.upper():
                    m[:] = True
                if self.exchange is not None and p == self.exchange.package.upper():
                    m |= self.exchange.map.entry_channel[:n] >= 0
                if self.exchange is not None and self.exchange.drains and p == self.exchange.drains[0].upper():
                    m |= self.exchange.drains[1].entry_channel[:n] >= 0
                for w in self.cfg.wells:
                    if w.package.upper() == p:
                        m[mf.entry(p, mf.node(*w.cell))] = True
                coupled[p] = m
            self._gw_packages = coupled
            L.zero("groundwater:closure_share", atol=1e-6)
        net = gross = 0.0
        for p, m in self._gw_packages.items():
            q = mf.simvals(p)
            L.enter(f"groundwater:boundary:{p.lower()}", float(q[~m[: len(q)]].sum()))
            net += float(q.sum())
            gross += float(np.abs(q).sum())
        sto = 0.0
        for v in ("STRGSS", "STRGSY"):
            try:
                x = mf.get(v, mf.model, "STO")
            except Exception:
                continue
            sto += float(np.sum(x))
        L.enter("groundwater:storage_release", sto)  # MODFLOW 6's STRGSS + STRGSY: positive is released into the flow
        # against the day's gross flow, with a floor of 1 m3 so a day with almost no flow isn't noise over noise
        L.enter("groundwater:closure_share", (net + sto) / max(gross + abs(sto), 1.0))

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
        (out / "summary.json").write_text(json.dumps(summary, indent=2))
        return out

    def finalize(self):
        for _, wu in self.users:
            wu.finalize()
        if self.mf6 is not None:
            self.mf6.finalize()
        self.swat.finalize()

    @classmethod
    def spin_up(cls, cfg: C.ExampleConfig, days: int = 365, cycles: int = 5, tol: float = 0.01, run: str = "spinup",
                users: bool = False, prepare=None, heads: np.ndarray | None = None) -> tuple[np.ndarray, list[dict]]:
        """MODFLOW 6's starting heads under the coupling: the first ``days`` run again and again, SWAT+ from its own
        start each time and MODFLOW 6 from the heads the last cycle ended with, until the cycle's net release from
        storage is within ``tol`` of the recharge it received (or ``cycles`` run). A start that isn't consistent with
        the coupling drains or fills storage over the run (the Clam's first coupled run drained 49 million m3,
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
