"""Stream temperature: a reach heat balance after a coupled run (Andrew, Nov 6, 2026: "yes I like that").

Each reach, each day, upstream first:

1. **Mixing.** The water entering the reach is mixed by volume (m3), each source at its temperature [A]:
   upstream reaches at the temperature this module gave them; lakes at ``5 + 0.75 x`` the day's air temperature
   (what SWAT+ gives any inflow that isn't a channel); water users' returns through point sources at their own
   temperature (the ledger's ``route:<point source>:heat_m3C`` over its volume); surface runoff, lateral flow and
   snowmelt at
   ``ch_temp.f90``'s lagged air temperatures (0.75 x the 5-, 10- and 2-day means, melt capped at 2.5 C); and
   **groundwater from MODFLOW 6 at a stated temperature**, by default the mean annual air temperature at the reach's
   station over the run [A] (shallow groundwater tracks the mean annual air temperature, typically within 1 to 2 C;
   Kurylyk and others 2014, to check). SWAT+'s own ``ch_temp`` can't do this step: no aquifer routine fills the
   groundwater term it reads (``hdsep%flo_gwsw``), and it mixes HRU water yields in mm with flows in m3/s.
2. **Heat exchange with the air** over the reach's travel time, by ``ch_temp``'s equilibrium-temperature method
   (Edinger and others 1974) and coefficients [A], integrated exactly (an exponential approach to the equilibrium
   temperature) where ``ch_temp`` takes one explicit step that overshoots it on long reaches: the exchange coefficient from wind and temperature, the
   equilibrium temperature from solar and long-wave radiation and the dew point, the flow's depth from the channel's
   geometry and ``hex_coef2``, travel time from the reach length and a Manning velocity at the Fuller and Peck peak
   rate (SWAT+ takes it from its rating curve, which swatplus-bmi doesn't expose) [A]. The dew point is
   ``tave - (100 - RH%) / 5`` (Lawrence 2005); SWAT+'s ``tave - (1 - rhum) / 5`` takes RH as a fraction and sits
   about 100 times too close to the air temperature (``swat_dew_point`` uses it).
3. **Energy.** ``reach_heat_balance`` returns the heat each source brings, the heat exchanged and the heat leaving,
   in J relative to 0 C; they close by construction and a test holds them to it. Water at 0 C stays at least
   0.001 C (as ``ch_temp``); the heat that floor adds is returned too (ice isn't modelled). A loss to the aquifer
   leaves at the mixed temperature and doesn't change it.

It reads the run's SWAT+ outputs and the ledger, and needs SWAT+ to print ``lsunit_wb`` daily (``print_daily``).
"""

from __future__ import annotations

import csv
import datetime as dt
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class Coefficients:
    """``ch_temp``'s defaults (channel_data_module.f90, water_temperature_data)."""
    gw_temperature: float | None = None  # C; None: the station's mean annual air temperature over the run [A]
    sno_lag: int = 2
    surf_lag: int = 5
    lat_lag: int = 10
    surf_lag_coef: float = 0.75   # also snowmelt
    lat_lag_coef: float = 0.75
    hex_coef1: float = 0.75       # dew point
    hex_coef2: float = 1.5        # channel geometry
    ssff: float = 0.5             # shade
    pk_rto: float = 1.0           # SWAT+'s peak rate ratio (sd_ch%pk_rto) [A]
    swat_dew_point: bool = False
    heat_exchange: bool = True    # off: the mixed temperature only, for checks
    humidity_if_generated: float = 0.7  # relative humidity (fraction) where SWAT+ generates it [A]
    seepage_as_runoff: bool = True  # drains' seepage reaches the stream over land, at the runoff's temperature [A]


def print_daily(project: Path, objects: list[str]):
    """Switch SWAT+'s daily print on for the given objects (``print.prt`` rows, e.g. ``lsunit_wb``), from the first
    day: ``nyskip`` (years not printed) is set to 0, since a balance or a heat balance needs every day, and the rows
    are read by their labels."""
    path = Path(project) / "print.prt"
    lines = path.read_text().splitlines()
    f = lines[2].split()
    if f and f[0] != "0":
        lines[2] = "  ".join(["0"] + f[1:])
    # SWAT+ 62.0.1 reads the object rows by position unless use_obj_labels (the csvout line's second flag) is y
    # (basin_print_codes_read.f90); by label, a row edited by its name is the row SWAT+ applies
    for i, ln in enumerate(lines[:-1]):
        if ln.split()[:1] == ["csvout"]:
            v = lines[i + 1].split()
            if len(v) >= 2 and v[1] != "y":
                v[1] = "y"
                lines[i + 1] = "  ".join(v)
            break
    for i, ln in enumerate(lines):
        f = ln.split()
        if f and f[0] in objects and len(f) == 5:
            lines[i] = f"{f[0]:<20s}{'y':>14s}{f[2]:>14s}{f[3]:>14s}{f[4]:>14s}"
    path.write_text("\n".join(lines) + "\n")


def _rows(path: Path, skip: int = 1):
    lines = Path(path).read_text().splitlines()
    head = lines[skip].split()
    rows = [ln.split() for ln in lines[skip + 1:] if ln.split()]
    if "jday" in head:  # SWAT+'s output files carry a units row under the header
        rows = [r for r in rows if r[0].lstrip("-").isdigit()]
    return head, rows


def _weather(swat: Path, name: str, humidity_if_generated: float) -> dict:
    """A station's daily series, keyed by date: tmax, tmin, solar (MJ/m2), relative humidity (fraction), wind (m/s),
    from the files ``weather-sta.cli`` names, and its latitude. A variable SWAT+ generates (``sim``) has no series:
    humidity then takes ``humidity_if_generated`` [A]; the others are needed."""
    lines = (swat / "weather-sta.cli").read_text().splitlines()
    head = lines[1].split()
    files = dict(zip(head, next(ln.split() for ln in lines[2:] if ln.split() and ln.split()[0] == name)))
    out, lat = {}, None
    for kind, cols in (("tmp", ("tmax", "tmin")), ("slr", ("slr",)), ("hmd", ("hmd",)), ("wnd", ("wnd",))):
        fname = files.get(kind, "sim")
        if fname in ("sim", "null"):
            if kind != "hmd":
                raise ValueError(f"station {name!r}: SWAT+ generates its {kind}; the heat balance needs a series")
            continue
        flines = (swat / fname).read_text().splitlines()
        lat = float(flines[2].split()[2])
        for ln in flines[3:]:
            f = ln.split()
            if not f:
                continue
            d = dt.date(int(f[0]), 1, 1) + dt.timedelta(days=int(f[1]) - 1)
            rec = out.setdefault(d, {})
            for c, v in zip(cols, f[2:]):
                rec[c] = float(v)
    for rec in out.values():
        rec.setdefault("hmd", humidity_if_generated)
    return {"lat": lat, "days": out}


def _solradmx(lat_deg: float, jday: int) -> float:
    """SWAT+'s maximum possible solar radiation, MJ/m2 (cli_clgen.f90)."""
    latr = math.radians(lat_deg)
    sd = math.asin(0.4 * math.sin((jday - 82.0) / 58.09))
    dd = 1.0 + 0.033 * math.cos(jday / 58.09)
    sdlat = -math.sin(latr) * math.tan(sd) / math.cos(latr)
    h = 0.0 if sdlat > 1 else math.acos(sdlat) if sdlat >= -1 else 3.1416
    ys, yc = math.sin(latr) * math.sin(sd), math.cos(latr) * math.cos(sd)
    return 30.0 * dd * (h * ys + yc * math.sin(h))


def _routing(swat: Path) -> tuple[dict[str, str | None], list[str]]:
    """Each channel's and lake's downstream object, and the channels upstream first."""
    down = {}
    for f in ("chandeg.con", "reservoir.con"):
        p = swat / f
        if not p.exists():
            continue
        head, rows = _rows(p)
        for r in rows:
            typ, i = r[-4], r[-3]
            down[r[1]] = {"sdc": f"cha{int(i):03d}", "res": f"res{int(i):03d}"}.get(typ)
    order, seen = [], set()

    def visit(o):
        if o in seen:
            return
        seen.add(o)
        for u, dn in down.items():
            if dn == o:
                visit(u)
        order.append(o)

    for o in down:
        visit(o)
    return down, [o for o in order if o.startswith("cha")]


RHO_C = 1000.0 * 4186.0  # J m-3 K-1, water (as ch_temp)


def reach_heat_balance(parts, tave: float, rhum: float, wind: float, solar: float, solar_max: float, geo: dict,
                       q_in: float, q_out: float, area_ha: float, coef: Coefficients) -> dict:
    """One reach, one day. ``parts``: [(m3, C), ...], the water entering and its temperatures; the weather (air
    temperature C, relative humidity as a fraction, wind m/s, solar and its clear-sky maximum MJ/m2); the reach's
    geometry (``wd``, ``dp``, ``slp``, ``len`` km, ``mann``); its inflow and outflow (m3/s) and drainage area (ha).

    Returns the mixed temperature ``t_init``, the equilibrium temperature ``t_eq``, the outflow temperature
    ``t_out``, and the energy relative to 0 C (J): ``in_J`` brought by the sources, ``exchanged_J`` with the air
    (the surface flux integrated over the parcel's surface, V / depth, and its travel time), ``floor_J`` added by the
    0.001 C floor, and ``out_J`` leaving with the mixed volume. in + exchanged + floor = out."""
    vol = sum(v for v, _ in parts)
    in_j = RHO_C * sum(v * t for v, t in parts)
    t_init = in_j / (RHO_C * vol) if vol > 0 else 5.0 + 0.75 * tave
    dewpt = tave - (1.0 - rhum) / 5.0 if coef.swat_dew_point else tave - (100.0 - 100.0 * rhum) / 5.0
    tdx = min(dewpt * coef.hex_coef1, tave)
    t_md = (t_init + tdx) / 2.0
    ke_beta = 0.35 + 0.015 * t_md + 0.0012 * t_md ** 2
    f_wind = 9.2 + 0.46 * wind ** 2
    k_e = 4.48 + 0.05 * t_init + (ke_beta + 0.47) * f_wind                  # W m-2 K-1
    h_sr = 0.97 * (solar * 11.57) * (1.0 - coef.ssff)                       # MJ m-2 d-1 to W m-2
    e_a = 6.1275 * math.exp(17.62 * tave / (237.3 + tave)) * rhum
    cloud = 0.0 if solar_max < 1e-4 else 0.9 * (solar / solar_max) + 0.1
    e_atm = 0.74 + 0.0065 * e_a * (1.0 + 0.17 * cloud ** 2)
    h_atm = 0.96 * e_atm * 5.67e-8 * (tave + 273.15) ** 4
    t_eq = tdx + h_sr / k_e + (h_atm - 305.5 - 4.48 * tdx) / k_e
    pk = coef.pk_rto * (1.0 + 2.66 * (area_ha / 100.0) ** -0.3) * q_in     # Fuller and Peck peak rate, m3/s
    wd, dp, slp = geo["wd"], geo["dp"], geo["slp"]
    dep = (max(pk, 1e-6) * geo["mann"] / (wd * math.sqrt(slp))) ** 0.6 if slp > 0 else 0.1
    vc = max(0.01, max(pk, 1e-6) / (wd * max(dep, 1e-3)))                  # Manning velocity, m/s [A]
    wid_chan = wd if wd >= 0.1 else wd * 0.25
    wid_flow = math.sqrt(((q_out / vc) * 2) / ((dp / coef.hex_coef2) / wid_chan)) if q_out > 0 else 0.0
    dep_flow = min(max(0.1, ((dp / coef.hex_coef2) / wid_chan) * wid_flow), dp * 2)
    tau = geo["len"] * 1000.0 / vc                                          # travel time, s
    # the exchange integrated over the travel time: dT/dt = k_e (T_eq - T) / (rho c d) gives an exponential approach
    # to T_eq. ch_temp takes one explicit step, k_e (T_eq - T) tau / (rho c d), which overshoots T_eq once tau is long
    # against rho c d / k_e (Osu_1hru's channel reached 90 C); the exact solution can't
    rate = k_e * tau / (RHO_C * dep_flow)
    dt_ex = (t_eq - t_init) * -math.expm1(-rate) if coef.heat_exchange else 0.0
    exchanged = RHO_C * vol * dt_ex                                         # J through the parcel's surface
    t = t_init + dt_ex
    floor = 0.0
    if t < 0.001:
        floor = RHO_C * vol * (0.001 - t)
        t = 0.001
    out = RHO_C * vol * t
    return {"t_init": t_init, "t_eq": t_eq, "t_out": t, "in_J": in_j, "exchanged_J": exchanged, "floor_J": floor,
            "out_J": out, "volume_m3": vol, "travel_s": tau, "depth_m": dep_flow}


def stream_temperature(run: str | Path, coef: Coefficients | None = None) -> dict:
    """Each channel's outflow temperature each day of a finished run (``run`` holds ``swatplus/`` and
    ``outputs/ledger.csv``). Returns {"dates": [...], "channels": [...], "temperature": array days x channels, C,
    "swat_ch_temp": SWAT+'s own ch_temp result for comparison, "sources": volumes by source, m3}."""
    coef = coef or Coefficients()
    run = Path(run)
    swat = run / "swatplus"

    # channels: station, geometry, routing unit draining to each (rtu k -> cha k, as SWAT+'s builders write them)
    head, rows = _rows(swat / "chandeg.con")
    chans = [r[1] for r in rows]
    wst = {r[1]: r[head.index("wst")] for r in rows}
    # each channel's geometry: chandeg.con's lcha row in channel-lte.cha names its hyd-sed-lte.cha row (cha_hyd)
    h2, r2 = _rows(swat / "hyd-sed-lte.cha")
    hyd = {r[0]: {k: float(r[h2.index(k)]) for k in ("wd", "dp", "slp", "len", "mann")} for r in r2}
    h3, r3 = _rows(swat / "channel-lte.cha")
    lte = {r[h3.index("id")]: r[h3.index("cha_hyd")] for r in r3}
    geo = {r[1]: hyd[lte[r[head.index("lcha")]]] for r in rows}
    area_ha = {r[1]: float(r[head.index("area")]) for r in rows}  # the channel's drainage area, for the peak rate
    lsu_area = {}
    hl, rl = _rows(swat / "ls_unit.def", skip=2)
    for r in rl:
        lsu_area[r[1]] = float(r[2]) * 1e4  # m2
    down, order = _routing(swat)

    # SWAT+ outputs by day
    def daily(path, key_col, cols):
        h, rr = _rows(path)
        idx = {c: h.index(c) for c in cols}
        out = {}
        for r in rr:
            d = dt.date(int(r[3]), 1, 1) + dt.timedelta(days=int(r[0]) - 1)
            out.setdefault(d, {})[r[h.index(key_col)]] = {c: float(r[i]) for c, i in idx.items()}
        return out

    cha = daily(swat / "channel_sd_day.txt", "name", ["flo_in", "flo_out"])
    hc, rc = _rows(swat / "channel_sd_day.txt")
    swat_t = {}
    for r in rc:  # ch_temp's own result is the unlabeled column before water_temp (ch_out_d%temp)
        d = dt.date(int(r[3]), 1, 1) + dt.timedelta(days=int(r[0]) - 1)
        swat_t.setdefault(d, {})[r[hc.index("name")]] = float(r[-2])
    lsu = daily(swat / "lsunit_wb_day.txt", "name", ["precip", "snofall", "snomlt", "surq_gen", "latq"])
    ledger, seep_made = {}, {}
    with open(run / "outputs" / "ledger.csv", newline="") as f:
        for r in csv.DictReader(f):
            d = dt.date.fromisoformat(r["date"])
            ledger[d] = {c: float(r[f"exchange:{c}:applied"]) for c in chans if f"exchange:{c}:applied" in r}
            seep_made[d] = {c: float(r[f"exchange:{c}:seepage"]) for c in chans if f"exchange:{c}:seepage" in r}
    dates = sorted(set(cha) & set(lsu) & set(ledger))
    # water users' returns into channels through point sources: volume and heat (m3 C) per channel and day
    returns = {}
    if (swat / "recall.con").exists():
        hr, rr = _rows(swat / "recall.con")
        to_ch = {r[1]: f"cha{int(r[-3]):03d}" for r in rr if r[-4] == "sdc"}
        with open(run / "outputs" / "ledger.csv", newline="") as f:
            for r in csv.DictReader(f):
                d = dt.date.fromisoformat(r["date"])
                for name, ch in to_ch.items():
                    v = r.get(f"route:swat:point_source:{name}:delivered")
                    if v in (None, "", "nan") or float(v) <= 0:
                        continue
                    h = r.get(f"route:swat:point_source:{name}:heat_m3C")
                    vol, heat = returns.get((d, ch), (0.0, 0.0))
                    returns[(d, ch)] = (vol + float(v), heat + (float(h) if h not in (None, "", "nan") else np.nan))
    stations = {s: _weather(swat, s, coef.humidity_if_generated) for s in set(wst.values())}

    def lagged(st, d, n):
        vals = []
        for k in range(n):
            w = st["days"].get(d - dt.timedelta(days=k))
            if w is None:
                break
            vals.append(0.5 * (w["tmax"] + w["tmin"]))
        return float(np.mean(vals))

    T = np.zeros((len(dates), len(chans)))
    S = np.zeros_like(T)
    names = ("upstream", "lakes", "surface", "lateral", "snowmelt", "groundwater", "returns")
    src = {k: np.zeros_like(T) for k in names}
    heat = {k: np.zeros_like(T) for k in ("in_J", "exchanged_J", "out_J", "floor_J")}
    col = {c: j for j, c in enumerate(chans)}
    t_gw_of = {}
    for name, st in stations.items():  # groundwater: the station's mean annual air temperature over the run [A]
        t_gw_of[name] = coef.gw_temperature if coef.gw_temperature is not None else float(
            np.mean([0.5 * (w["tmax"] + w["tmin"]) for w in st["days"].values()]))
    for i, d in enumerate(dates):
        temp_today = {}
        for c in order:
            st = stations[wst[c]]
            w = st["days"][d]
            tave = 0.5 * (w["tmax"] + w["tmin"])
            t_surf = max(0.01, lagged(st, d, coef.surf_lag)) * coef.surf_lag_coef
            t_lat = max(0.01, lagged(st, d, coef.lat_lag)) * coef.lat_lag_coef
            t_sno = min(2.5, max(0.01, lagged(st, d, coef.sno_lag))) * coef.surf_lag_coef
            t_lake = 5.0 + 0.75 * tave
            rtu = "rtu" + c[3:]
            u = lsu[d].get(rtu, {"precip": 0, "snofall": 0, "snomlt": 0, "surq_gen": 0, "latq": 0})
            a = lsu_area.get(rtu, 0.0)
            v_surf, v_lat = u["surq_gen"] / 1000 * a, u["latq"] / 1000 * a
            rain = max(u["precip"] - u["snofall"], 0.0)
            melt = u["snomlt"] / (u["snomlt"] + rain) if (u["snomlt"] + rain) > 0 else 0.0
            v_gw = max(ledger[d].get(c, 0.0), 0.0)
            # yesterday's seepage arrives today with the rest of the exchange (the one-day lag)
            v_seep = min(max(seep_made.get(d - dt.timedelta(days=1), {}).get(c, 0.0), 0.0), v_gw) \
                if coef.seepage_as_runoff else 0.0
            v_gw -= v_seep
            ups = [x for x, dn in down.items() if dn == c and x.startswith("cha")]
            v_up = sum(cha[d][x]["flo_out"] * 86400 for x in ups)
            t_up = (sum(cha[d][x]["flo_out"] * 86400 * temp_today[x] for x in ups) / v_up) if v_up > 0 else 0.0
            v_in = cha[d][c]["flo_in"] * 86400
            v_ret, heat_ret = returns.get((d, c), (0.0, 0.0))
            t_ret = heat_ret / v_ret if v_ret > 0 and np.isfinite(heat_ret) else t_surf  # [A] if no temperature given
            v_lake = max(v_in - v_up - v_surf - v_lat - v_gw - v_seep - v_ret, 0.0) if any(
                dn == c and x.startswith("res") for x, dn in down.items()) else 0.0
            parts = [(v_up, t_up), (v_lake, t_lake), ((1 - melt) * v_surf + v_seep, t_surf),
                     ((1 - melt) * v_lat, t_lat), (melt * (v_surf + v_lat), t_sno), (v_gw, t_gw_of[wst[c]]),
                     (v_ret, t_ret)]
            jd = d.timetuple().tm_yday
            r = reach_heat_balance(parts, tave, w["hmd"], w["wnd"], w["slr"], _solradmx(st["lat"], jd), geo[c],
                                   cha[d][c]["flo_in"], cha[d][c]["flo_out"], area_ha[c], coef)
            temp_today[c] = r["t_out"]
            j = col[c]
            T[i, j] = r["t_out"]
            S[i, j] = swat_t.get(d, {}).get(c, np.nan)
            for k, (v, _) in zip(names, parts):
                src[k][i, j] = v
            for k in heat:
                heat[k][i, j] = r[k]
    return {"dates": [d.isoformat() for d in dates], "channels": chans, "temperature": T, "swat_ch_temp": S,
            "sources": src, "heat": heat, "groundwater_temperature": t_gw_of}


def write(result: dict, path: str | Path):
    """temperature.csv: date, then one column per channel, C."""
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["date", *result["channels"]])
        for d, row in zip(result["dates"], result["temperature"]):
            w.writerow([d, *(f"{v:.4f}" for v in row)])
