"""What each scenario of water users does to the basin, and whether users' effects add up (docs/design.md,
"Scenarios of water users").

Every case of a set (``scenarios.load_set``) is the baseline with its changes, run from one warm-up through the same
weather. Over the compared period each case is set against the baseline: the flow at the channels compared (daily,
its mean, its flow duration curve and its low flows: the lowest n-day mean, each climatic year's lowest, and the
nQT by log-Pearson type III), the heads at the cells compared, the capture of its extra pumping
(``analysis.capture``), and the water of each user and system it changes. For a combination, the interaction is its
effect less the sum of its scenarios' effects alone.

``write_tables(set)`` writes the CSVs and ``tables.md`` to the set's ``tables/``; the functions it uses work on plain
arrays and are tested on their own.
"""

from __future__ import annotations

import csv
import datetime as dt
import json
import math
from pathlib import Path
from statistics import NormalDist

import numpy as np

from .analysis import capture, read_ledger

BASELINE = "baseline"         # the case that is the example itself
EXCEEDANCE = (0.10, 0.50, 0.90, 0.95, 0.99)
CFS = 0.028316846592          # m3 s-1 per ft3 s-1
FT = 0.3048                   # m per ft
MG = 3785.411784              # m3 per million gallons
UNITS = {  # name: (flow, head, volume, rate) labels and factors from SI
    "SI": {"flow": ("m3/s", 1.0), "head": ("m", 1.0), "volume": ("m3", 1.0), "rate": ("m3/s", 1.0)},
    "US": {"flow": ("cfs", 1 / CFS), "head": ("ft", 1 / FT), "volume": ("MG", 1 / MG), "rate": ("MGD", 1 / MG)},
}


# ---------- statistics on a daily series ----------

def rolling_mean(q: np.ndarray, days: int) -> np.ndarray:
    """The mean of each ``days``-day window ending on each day; NaN until the first full window."""
    q = np.asarray(q, float)
    out = np.full(len(q), np.nan)
    if len(q) >= days:
        c = np.cumsum(np.insert(q, 0, 0.0))
        out[days - 1:] = (c[days:] - c[:-days]) / days
    return out


def climatic_years(dates: list[dt.date], year_start: tuple[int, int] = (4, 1)) -> dict[int, np.ndarray]:
    """Each complete climatic year in ``dates`` (named by the calendar year it starts in) -> its days' indices."""
    def year(d):
        return d.year if (d.month, d.day) >= year_start else d.year - 1
    ys: dict[int, list[int]] = {}
    for i, d in enumerate(dates):
        ys.setdefault(year(d), []).append(i)
    out = {}
    for y, idx in ys.items():
        first = dt.date(y, *year_start)
        n = (dt.date(y + 1, *year_start) - first).days
        if len(idx) == n and dates[idx[0]] == first:
            out[y] = np.array(idx)
    return out


def lp3_low_flow(minima: np.ndarray, return_period: float) -> float:
    """The low flow not exceeded once in ``return_period`` years on average, by log-Pearson type III on the annual
    minima (all > 0): the logs' mean, standard deviation and sample skew, the frequency factor by Wilson and
    Hilferty's approximation (as USGS practice)."""
    x = np.log10(np.asarray(minima, float))
    n = len(x)
    m, s = x.mean(), x.std(ddof=1)
    if s == 0:
        return float(10 ** m)
    g = n / ((n - 1) * (n - 2)) * np.sum(((x - m) / s) ** 3)
    z = NormalDist().inv_cdf(1.0 / return_period)
    k = z if abs(g) < 1e-9 else (2.0 / g) * ((1.0 + g * z / 6.0 - g * g / 36.0) ** 3 - 1.0)
    return float(10 ** (m + k * s))


def empirical_low_flow(minima: np.ndarray, return_period: float) -> float:
    """The same quantile from the minima's Weibull plotting positions (i / (n + 1)), interpolated; NaN beyond them."""
    x = np.sort(np.asarray(minima, float))
    p = np.arange(1, len(x) + 1) / (len(x) + 1)
    target = 1.0 / return_period
    return float(np.interp(target, p, x)) if p[0] <= target <= p[-1] else float("nan")


def low_flows(q: np.ndarray, dates: list[dt.date], days: int = 7, return_period: float = 10,
              year_start: tuple[int, int] = (4, 1)) -> dict:
    """The lowest ``days``-day mean in the period; each complete climatic year's lowest and their mean; and the
    ``days``Q``return_period`` (7Q10 by default): log-Pearson type III, or empirical where a year's minimum is zero,
    NaN with fewer than 3 years."""
    r = rolling_mean(q, days)
    years = climatic_years(dates, year_start)
    minima = np.array([np.nanmin(r[idx]) for idx in years.values() if np.isfinite(r[idx]).any()])
    out = {"lowest": float(np.nanmin(r)) if np.isfinite(r).any() else float("nan"),
           "annual_minima": minima.tolist(), "years": len(minima),
           "mean_annual_minimum": float(minima.mean()) if len(minima) else float("nan")}
    if len(minima) < 3:
        out |= {"nQT": float("nan"), "method": f"fewer than 3 climatic years ({len(minima)})"}
    elif (minima <= 0).any():
        out |= {"nQT": empirical_low_flow(minima, return_period), "method": "empirical (a year's minimum is zero)"}
    else:
        out |= {"nQT": lp3_low_flow(minima, return_period), "method": "log-Pearson III"}
    return out


def flow_duration(q: np.ndarray, exceedance) -> np.ndarray:
    """The flow exceeded on each share of days (0.9: the flow exceeded 90% of the time)."""
    return np.quantile(np.asarray(q, float), 1.0 - np.asarray(exceedance, float))


def interaction(together: float, alone: list[float]) -> dict:
    """A combination's effect against the sum of its parts' effects alone."""
    s = float(np.sum(alone))
    i = together - s
    return {"sum_alone": s, "together": together, "interaction": i,
            "interaction_share": i / together if together else float("nan")}


def daily_interaction(together: np.ndarray, alone: list[np.ndarray]) -> dict:
    """Day by day: the share of the combination's change that isn't the sum of the parts' (sum |I| / sum |E|) and
    the largest |I| on a day."""
    i = np.asarray(together) - np.sum(alone, axis=0)
    e = np.nansum(np.abs(together))
    return {"interaction_share": float(np.nansum(np.abs(i)) / e) if e else float("nan"),
            "largest": float(np.nanmax(np.abs(i))) if len(i) else float("nan")}


# ---------- a set's cases ----------

def read_case(folder: Path) -> dict:
    """A finished case's outputs: dates, ledger and series."""
    out = Path(folder) / "outputs"
    with open(out / "series.csv", newline="") as f:
        dates = [dt.date.fromisoformat(r["date"]) for r in csv.DictReader(f)]
    return {"dates": dates, "ledger": read_ledger(out / "ledger.csv"), "series": read_ledger(out / "series.csv"),
            "summary": json.loads((out / "summary.json").read_text())}


def _period(s, dates: list[dt.date]) -> slice:
    """The compared days: from the warm-up's end (or [compare] from) to the run's end (or [compare] to)."""
    done = s.folder / "warm_up" / "done.json"
    day = json.loads(done.read_text())["day"] if done.exists() and not s.warm_up_state else \
        json.loads((s.state / "basin.json").read_text())["day"]
    i0, i1 = day, len(dates)
    if s.compare.get("from"):
        i0 = max(i0, next((i for i, d in enumerate(dates) if d >= dt.date.fromisoformat(str(s.compare["from"]))), i1))
    if s.compare.get("to"):
        last = dt.date.fromisoformat(str(s.compare["to"]))
        i1 = min(i1, sum(1 for d in dates if d <= last))
    if i1 <= i0:
        raise ValueError(f"{s.path}: the compared period is empty")
    return slice(i0, i1)


def _changed_users(s, case: str) -> list[str]:
    out = []
    for m in s.members(case):
        out += [u.name for u in m.users] + [x["target"][1] for x in m.scale if x["target"][0] == "user"]
        out += [r.split(":", 1)[1] for r in m.remove if r.startswith("user:")]
    return list(dict.fromkeys(out))


def _total(t: dict, k: str, p: slice) -> float | None:
    return float(np.nansum(t[k][p])) if k in t else None


def _user_use(t: dict, u: str, p: slice) -> dict:
    k = f"user:{u}:"
    def first(*names):
        for n in names:
            v = _total(t, k + n, p)
            if v is not None:
                return v
        return 0.0
    out = {"withdrawal": first("withdrawal"), "delivered": first("delivered", "withdrawal"), "unmet": first("unmet"),
           "consumption": first("consumption_met", "consumption"), "return": first("return_due", "return")}
    if k + "MWh" in t:
        out["MWh"] = first("MWh")
    return out


def _systems(t: dict, p: slice) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    for k in t:
        if not k.startswith("utility:"):
            continue
        _, name, acct = k.split(":", 2)
        if acct in ("demand", "produced", "delivered", "shortfall", "outfall", "reuse") or \
                (acct.startswith("source:") and acct.endswith(":delivered")):
            out.setdefault(name, {})[acct] = _total(t, k, p)
    return out


def compare(s) -> dict:
    """Every case against the baseline over the compared period: the rows of each table, in SI units."""
    cases = {c: read_case(s.folder / "cases" / c) for c in s.cases
             if (s.folder / "cases" / c / "done.json").exists()}
    if BASELINE not in cases:
        raise FileNotFoundError(f"{s.folder}: the baseline hasn't finished; nothing to compare")
    base = cases[BASELINE]
    p = _period(s, base["dates"])
    dates = base["dates"][p]
    lf = s.compare.get("low_flow", {})
    days, T = int(lf.get("days", 7)), float(lf.get("return_period", 10))
    ys = tuple(int(x) for x in str(lf.get("year_start", "04-01")).split("-"))
    exc = [float(x) for x in s.compare.get("exceedance", EXCEEDANCE)]
    nqt = f"{days}Q{T:g}"
    channels = s.compare.get("channels") or [k.split(":", 1)[1] for k in base["series"] if k.startswith("flow:")]
    heads = [h["name"] for h in s.compare.get("heads", [])] or \
        [k.split(":", 1)[1] for k in base["series"] if k.startswith("head:")]

    def flow_stats(q):
        lows = low_flows(q, dates, days, T, ys)
        st = {"mean": float(np.mean(q))}
        st |= {f"Q{round(100 * e):g}": float(v) for e, v in zip(exc, flow_duration(q, exc))}
        st |= {f"lowest {days}-day": lows["lowest"], f"mean annual {days}-day minimum": lows["mean_annual_minimum"],
               nqt: lows["nQT"]}
        return st, lows

    def head_stats(h):
        return {"mean": float(np.mean(h)), "lowest": float(np.min(h)), "last day": float(h[-1])}

    metrics: dict[str, dict] = {}   # case -> {(kind, point, statistic): value}, for the interaction
    daily: dict[str, dict] = {}     # case -> {(kind, point): daily values}
    rows = {k: [] for k in ("flow", "heads", "capture", "use", "systems", "flow_duration")}
    notes = {}
    for c, r in cases.items():
        if r["dates"][p] != dates:
            raise ValueError(f"case {c} covers other days than the baseline")
        mt, dl = {}, {}
        for ch in channels:
            q = r["series"][f"flow:{ch}"][p]
            st, lows = flow_stats(q)
            notes[(c, ch)] = f"{lows['years']} climatic years; {lows['method']}"
            for k, v in st.items():
                mt[("flow", ch, k)] = v
            dl[("flow", ch)] = q
        for h in heads:
            x = r["series"][f"head:{h}"][p]
            for k, v in head_stats(x).items():
                mt[("head", h, k)] = v
            dl[("head", h)] = x
        for u in {u for cc in s.cases if cc in cases for u in _changed_users(s, cc)}:
            for k, v in _user_use(r["ledger"], u, p).items():
                mt[("use", u, k)] = v
        for name, accts in _systems(r["ledger"], p).items():
            for k, v in accts.items():
                mt[("system", name, k)] = v
        metrics[c], daily[c] = mt, dl
    bm, bd = metrics[BASELINE], daily[BASELINE]

    def eff(c, key):
        v, b = metrics[c].get(key), bm.get(key, 0.0)
        return None if v is None else v - b

    for c in cases:
        for ch in channels:
            for k in [k for (kind, pt, k) in bm if kind == "flow" and pt == ch]:
                b, v = bm[("flow", ch, k)], metrics[c][("flow", ch, k)]
                rows["flow"].append({"case": c, "channel": ch, "statistic": k, "baseline": b, "value": v,
                                     "change": v - b, "change_percent": 100 * (v - b) / b if b else float("nan"),
                                     "note": notes[(c, ch)] if k == nqt else ""})
        for h in heads:
            d = daily[c][("head", h)] - bd[("head", h)]
            for k in ("mean", "lowest", "last day"):
                b, v = bm[("head", h, k)], metrics[c][("head", h, k)]
                rows["heads"].append({"case": c, "head": h, "statistic": k, "baseline": b, "value": v,
                                      "change": v - b, "note": ""})
            i = int(np.argmin(d))
            rows["heads"].append({"case": c, "head": h, "statistic": "largest drawdown", "baseline": float("nan"),
                                  "value": float("nan"), "change": float(d[i]), "note": dates[i].isoformat()})
        if c != BASELINE:
            cap = capture({k: v[p] for k, v in base["ledger"].items()}, {k: v[p] for k, v in cases[c]["ledger"].items()})
            for part, v in cap["total"].items():
                metrics[c][("capture", part, "total")] = v
                rows["capture"].append({"case": c, "part": part, "m3": v,
                                        "share": cap["share"].get(part, float("nan"))})
        for u in _changed_users(s, c) if c != BASELINE else []:
            for k, v in _user_use(cases[c]["ledger"], u, p).items():
                b = bm.get(("use", u, k), 0.0)
                rows["use"].append({"case": c, "user": u, "account": k, "baseline": b, "value": v, "change": v - b,
                                    "days": len(dates)})
        if c != BASELINE:
            for (kind, name, k), v in metrics[c].items():
                if kind != "system":
                    continue
                b = bm.get(("system", name, k), 0.0)
                if v != b:
                    rows["systems"].append({"case": c, "system": name, "account": k, "baseline": b, "value": v,
                                            "change": v - b, "days": len(dates)})
    for ch in channels:
        grid = np.arange(1, 100) / 100
        for e in grid:
            row = {"channel": ch, "exceedance": float(e)}
            for c in cases:
                row[c] = float(flow_duration(daily[c][("flow", ch)], [e])[0])
            rows["flow_duration"].append(row)

    inter = []
    for combo, members in s.combinations.items():
        if combo not in cases or any(m not in cases for m in members):
            continue
        for key in metrics[combo]:
            if key[0] == "use" or key not in bm and key[0] != "capture":
                continue
            alone = [eff(m, key) if key[0] != "capture" else metrics[m].get(key, 0.0) for m in members]
            tog = eff(combo, key) if key[0] != "capture" else metrics[combo][key]
            if tog is None or any(a is None for a in alone) or not np.isfinite([tog, *alone]).all() or \
                    key[1] == "residual" or (tog == 0 and not any(alone)):
                continue
            inter.append({"combination": combo, "kind": key[0], "point": key[1], "statistic": key[2],
                          **{f"alone:{m}": a for m, a in zip(members, alone)}, **interaction(tog, alone)})
        for key in daily[combo]:
            d = daily[combo][key] - bd[key]
            parts = [daily[m][key] - bd[key] for m in members]
            di = daily_interaction(d, parts)
            inter.append({"combination": combo, "kind": key[0], "point": key[1], "statistic": "daily",
                          "sum_alone": float("nan"), "together": float("nan"), "interaction": di["largest"],
                          "interaction_share": di["interaction_share"]})
        for u in _changed_users(s, combo):  # the users' own water: it adds exactly when nothing feeds back
            for k in ("withdrawal", "consumption", "return"):
                key = ("use", u, k)
                alone = [eff(m, key) if u in _changed_users(s, m) else 0.0 for m in members]
                inter.append({"combination": combo, "kind": "use", "point": u, "statistic": k,
                              **{f"alone:{m}": a for m, a in zip(members, alone)},
                              **interaction(eff(combo, key), alone)})
    rows["interaction"] = inter
    rows["daily"] = {"dates": dates, "series": {c: daily[c] for c in cases}}
    rows["status"] = base["summary"].get("status"), base["summary"].get("status_note")
    return rows


# ---------- writing ----------

def _write_csv(path: Path, rows: list[dict]):
    cols = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: (repr(v) if isinstance(v, float) else v) for k, v in r.items()})


def _g(x, digits=4) -> str:
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "–"
    return f"{x:.{digits}g}" if abs(x) >= 1e-12 or x == 0 else "0"


def _pct(x) -> str:
    return "–" if x is None or not math.isfinite(x) else f"{x:+.1f}%"


def _md(header: list[str], rows: list[list[str]]) -> str:
    return "\n".join(["| " + " | ".join(header) + " |", "|" + "---|" * len(header)] +
                     ["| " + " | ".join(r) + " |" for r in rows]) + "\n"


def write_tables(s) -> Path:
    """The set's tables (module docstring) to ``<folder>/tables``; returns the folder."""
    out = s.folder / "tables"
    out.mkdir(parents=True, exist_ok=True)
    r = compare(s)
    u = UNITS[s.compare.get("units", "SI")]
    fl, fq = u["flow"]
    hl, hq = u["head"]
    vl, vq = u["volume"]
    rl, rq = u["rate"]

    def conv(rows, cols, f):
        return [{k: (v * f if k in cols and isinstance(v, float) else v) for k, v in x.items()} for x in rows]

    _write_csv(out / "flow.csv", [{**x, "units": fl} for x in conv(r["flow"], ("baseline", "value", "change"), fq)])
    _write_csv(out / "flow_duration.csv", [{**x, "units": fl} for x in
                                           conv(r["flow_duration"], set(s.cases), fq)])
    _write_csv(out / "heads.csv", [{**x, "units": hl} for x in conv(r["heads"], ("baseline", "value", "change"), hq)])
    _write_csv(out / "capture.csv", [{"case": x["case"], "part": x["part"], vl: x["m3"] * vq, "share": x["share"]}
                                     for x in r["capture"]])
    for name in ("use", "systems"):
        rows = []
        for x in r[name]:
            y = conv([x], ("baseline", "value", "change"), 1.0 if x["account"] == "MWh" else vq)[0]
            y["units"] = "MWh" if x["account"] == "MWh" else vl
            if x["account"] != "MWh":  # the mean rate over the period
                y[f"change_mean_rate_{rl}"] = x["change"] / x["days"] * (1.0 if rl == "MGD" else 1 / 86400.0) * \
                    (1 / MG if rl == "MGD" else 1.0)
            rows.append(y)
        _write_csv(out / f"{name}.csv", rows)
    _write_csv(out / "interaction.csv", r["interaction"])
    dates, series = r["daily"]["dates"], r["daily"]["series"]
    for kind, label, f, fname in (("flow", fl, fq, "daily_flow.csv"), ("head", hl, hq, "daily_heads.csv")):
        keys = [k for k in series[BASELINE] if k[0] == kind]
        if not keys:
            continue
        with open(out / fname, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["date"] + [f"baseline:{k[1]} ({label})" for k in keys] +
                       [f"{c} change:{k[1]} ({label})" for c in series if c != BASELINE for k in keys])
            for i, d in enumerate(dates):
                w.writerow([d.isoformat()] + [repr(series[BASELINE][k][i] * f) for k in keys] +
                           [repr((series[c][k][i] - series[BASELINE][k][i]) * f)
                            for c in series if c != BASELINE for k in keys])
    (out / "tables.md").write_text(_markdown(s, r, u, dates))
    return out


def _markdown(s, r, u, dates) -> str:
    fl, fq = u["flow"]
    hl, hq = u["head"]
    vl, vq = u["volume"]
    status, note = r["status"]
    md = [f"# {s.name}: what each scenario does\n"]
    if s.description:
        md.append(s.description.strip() + "\n")
    if status != "calibrated":
        md.append(f"**Screening run.** {note or ''}\n")
    md.append(f"Compared: {dates[0]} to {dates[-1]} ({len(dates):,} days), each case against the baseline "
              f"(the example's existing system), from one warm-up and the same weather.\n")
    wu = s.state / "warm_up.json"
    if wu.exists():
        w = json.loads(wu.read_text())
        if not w.get("settled", True):
            md.append(f"**The warm-up's soils hadn't settled:** soil water changed {w['soil_water_change_mm']:+.0f} mm "
                      f"over its last {w['days']} days ({w['share']:.0%} of their precipitation). Every case carries "
                      "the trend.\n")
    cases = [c for c in s.cases if any(x["case"] == c for x in r["flow"] + r["heads"])]
    flow = {(x["case"], x["channel"], x["statistic"]): x for x in r["flow"]}
    stats = list(dict.fromkeys(x["statistic"] for x in r["flow"]))
    pick = [k for k in stats if k == "mean" or k in ("Q90", "Q95") or k.startswith("lowest") or "Q" in k[1:]
            and k[0].isdigit()]
    for ch in dict.fromkeys(x["channel"] for x in r["flow"]):
        md.append(f"## Flow at {ch} ({fl})\n")
        rows = [["baseline"] + [_g(flow[(BASELINE, ch, k)]["value"] * fq) for k in pick]]
        for c in cases:
            if c == BASELINE:
                continue
            rows.append([c] + [f"{_g(flow[(c, ch, k)]['change'] * fq)} ({_pct(flow[(c, ch, k)]['change_percent'])})"
                               for k in pick])
        md.append(_md(["case"] + pick, rows))
        md.append("The baseline's row is its flow; each case's is its change, and the change in percent.\n")
        n = flow[(BASELINE, ch, pick[-1])]["note"] if pick else ""
        if n:
            md.append(f"Low flows: {n}. Flow exceeded the share of days given (Q90: 90% of days).\n")
    heads = {(x["case"], x["head"], x["statistic"]): x for x in r["heads"]}
    for h in dict.fromkeys(x["head"] for x in r["heads"]):
        md.append(f"## Heads at {h} ({hl})\n")
        rows = [[c, _g(heads[(c, h, "mean")]["change"] * hq), _g(heads[(c, h, "largest drawdown")]["change"] * hq),
                 heads[(c, h, "largest drawdown")]["note"], _g(heads[(c, h, "last day")]["change"] * hq)]
                for c in cases if c != BASELINE]
        md.append(_md(["case", "mean change", "largest drawdown", "on", "last day change"], rows))
    if r["capture"]:
        cap = {(x["case"], x["part"]): x for x in r["capture"]}
        parts = [p for p in dict.fromkeys(x["part"] for x in r["capture"]) if p != "pumping"]
        md.append(f"## Capture of the extra pumping\n")
        rows = [[c, _g(cap[(c, "pumping")]["m3"] * vq)] + [_pct(100 * cap[(c, p)]["share"]).lstrip("+") for p in parts]
                for c in cases if (c, "pumping") in cap]
        md.append(_md(["case", f"extra pumping ({vl})"] + [p.replace("_", " ") for p in parts], rows))
    if r["use"]:
        md.append(f"## Water of the users added or changed (totals, {vl}; energy MWh)\n")
        use = {}
        for x in r["use"]:
            use.setdefault((x["case"], x["user"]), {})[x["account"]] = x
        accts = list(dict.fromkeys(x["account"] for x in r["use"]))
        rows = [[c, user] + [_g(v[a]["change"] * (1.0 if a == "MWh" else vq)) if a in v else "–" for a in accts]
                for (c, user), v in use.items()]
        md.append(_md(["case", "user"] + [a + " change" for a in accts], rows))
    if r["systems"]:
        md.append(f"## Systems: the change in their water ({vl})\n")
        rows = [[x["case"], x["system"], x["account"], _g(x["baseline"] * vq), _g(x["change"] * vq)]
                for x in r["systems"]]
        md.append(_md(["case", "system", "account", "baseline", "change"], rows))
    if r["interaction"]:
        md.append("## Do the effects add up?\n")
        md.append("Interaction = the effect together less the sum of the effects alone. Daily: the share of the "
                  "combination's day-by-day change that isn't the sum of its parts', and the largest on a day.\n")
        rows = []
        for x in r["interaction"]:
            f = {"flow": fq, "head": hq, "capture": vq, "use": vq, "system": vq}[x["kind"]]
            if x["statistic"] == "daily":
                rows.append([x["combination"], x["kind"], x["point"], "daily", "–", "–", _g(x["interaction"] * f),
                             _pct(100 * x["interaction_share"]).lstrip("+")])
            else:
                rows.append([x["combination"], x["kind"], x["point"], x["statistic"], _g(x["sum_alone"] * f),
                             _g(x["together"] * f), _g(x["interaction"] * f), _pct(100 * x["interaction_share"])])
        md.append(_md(["combination", "kind", "point", "statistic", "sum alone", "together", "interaction",
                       "share"], rows))
    return "\n".join(md)
