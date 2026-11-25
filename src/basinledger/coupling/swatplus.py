"""SWAT+ through swatplus-bmi: project edits an example needs, and the daily reads and sets the driver uses.

The edits add objects to a copy of a SWAT+ project before it starts: water allocation places of use (withdrawals from
a channel, set daily through ``place_of_use_water__volume_demand_rate``) and daily point sources (recall objects,
set daily through ``point_source_water__volume_flow_rate``). They follow the input readers of SWAT+ 62.0.1
(water_allocation_read.f90, recall_read.f90), as swatplus-bmi's tests do. Places of use are of type "out": SWAT+
delivers their water nowhere, so it leaves SWAT+ and the driver places it.
"""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import numpy as np

DEMAND = "place_of_use_water__volume_demand_rate"        # m3 s-1, set for the next day
WITHDRAWN = "place_of_use_water__volume_withdrawal_rate"  # m3 s-1, what SWAT+ took
POINT = "point_source_water__volume_flow_rate"           # m3 s-1, set for the next day
QOUT = "channel_water__volume_outflow_rate"
AQU_OUT = "aquifer_water__volume_outflow_rate"
DEPTH = "channel_water__depth"                                   # m, swatplus-bmi 0.2.0
GW_EXCHANGE = "channel_water__groundwater_volume_exchange_rate"  # m3 s-1, positive into the channel, 0.2.0


def _swatplus_bmi_version() -> tuple[int, int]:
    from importlib.metadata import version
    major, minor = (int(x) for x in version("swatplus-bmi").split(".")[:2])
    return major, minor


def has_groundwater_exchange() -> bool:
    """Whether the installed swatplus-bmi sets a channel's groundwater exchange directly (0.2.0 and later). Earlier
    versions route the exchange through a point source and a water allocation place of use per channel."""
    return _swatplus_bmi_version() >= (0, 2)


def places_of_use_fixed() -> bool:
    """Whether the installed swatplus-bmi patches SWAT+'s place of use (0.14.0 and later). SWAT+ 62.0.1's
    wallo_withdraw.f90 sizes a withdrawal from ``pou(ipou)%pod(ipod)%duty``, with ``ipou`` the place's position in
    its point of diversion's list and ``ipod`` the point of diversion's number. So a place of use on any channel but
    the first takes nothing, or a stray value past the end of an array, and places sharing a channel take each other's
    demands. Found on the Rappahannock study's intakes (Nov 24, 2026)."""
    return _swatplus_bmi_version() >= (0, 14)


def places_of_use_problems(uses: list[tuple[str, str]]) -> list[str]:
    """What an older SWAT+ would get wrong in places of use [(name, channel), ...]: with swatplus-bmi before 0.14.0,
    only places of use on one channel take their own demands (``places_of_use_fixed``)."""
    channels = list(dict.fromkeys(ch for _, ch in uses))
    if len(channels) <= 1 or places_of_use_fixed():
        return []
    major, minor = _swatplus_bmi_version()
    return [f"places of use on {len(channels)} channels ({', '.join(channels[:5])}{', ...' if len(channels) > 5 else ''}) "
            f"need swatplus-bmi 0.14.0 or later (installed: {major}.{minor}): SWAT+ 62.0.1 gives a place of use on any "
            "channel but the first another's demand or none (wallo_withdraw.f90), which 0.14.0 patches"]


def add_places_of_use(project: Path, uses: list[tuple[str, str]], rate_max: float = 100.0):
    """Water allocation places of use, each taking from one channel: [(name, channel), ...].

    Several places of use on one channel share its point of diversion, which lists them in the order given here, and
    SWAT+ serves them in that order: the first takes its demand from the day's outflow, the next what it asks of what
    is left, and so on, so the first listed is the senior (wallo_control.f90). Each gets its own demand only with
    swatplus-bmi 0.14.0 or later (``places_of_use_fixed``)."""
    project = Path(project)
    cio = (project / "file.cio").read_text().splitlines()
    cio = ["water_rights      place_of_use.wro  point_of_diver.wro  null" if l.startswith("water_rights") else l
           for l in cio]
    (project / "file.cio").write_text("\n".join(cio) + "\n")
    chans = _channels(project)
    pou = [f"place_of_use.wro: {len(uses)} places of use (Basin Ledger)", str(len(uses)),
           "num name typ typ_num pods pors dtbl_mx rate_max dtbl_pod_fr dtbl_por_fr"]
    # SWAT+ keeps one point of diversion per channel (sd_ch%wallo_pod), so the places of use on a channel share it
    order = list(dict.fromkeys(ch for _, ch in uses))
    for i, (name, ch) in enumerate(uses, 1):
        pou.append(f"{i} {name} out {i} 1 0 null {rate_max} null null")
        pou.append(f"1 {order.index(ch) + 1} {ch} cha {chans.index(ch) + 1} null 0 null 0.0 null 1000000. 1.0 n")
    pod = [f"point_of_diver.wro: one per channel (Basin Ledger)", str(len(order)),
           "num name typ typ_num pous pou_num pou_name pod_num typ typ_num right"]
    for k, ch in enumerate(order, 1):
        mine = [(i, name) for i, (name, c) in enumerate(uses, 1) if c == ch]
        pod.append(f"{k} {ch} cha {chans.index(ch) + 1} {len(mine)} "
                   + " ".join(f"{i} {name} 1 out {i} senior" for i, name in mine))
    (project / "place_of_use.wro").write_text("\n".join(pou) + "\n")
    (project / "point_of_diver.wro").write_text("\n".join(pod) + "\n")


def add_point_sources(project: Path, sources: list[tuple[str, str]]):
    """Daily point sources (recall objects) into channels, zero flow in the files: [(name, channel), ...]."""
    project = Path(project)
    start, end = _period(project)
    if not (project / "recall.con").exists():
        _no_recalls_yet(project, sources[0][1])
    con = (project / "recall.con").read_text().splitlines()
    hdr = con[:2]
    template = con[2].split()
    rows = []
    for i, (name, ch) in enumerate(sources, 1):
        f = list(template)
        c = _channels(project).index(ch) + 1
        f[0], f[1], f[2], f[7] = str(i), name, str(i), str(i)
        f[13], f[14] = "sdc", str(c)
        rows.append("  ".join(f))
    (project / "recall.con").write_text("\n".join(hdr + rows) + "\n")
    cnt = (project / "object.cnt").read_text().splitlines()
    f = cnt[2].split()  # name ls_area tot_area obj hru lhru rtu mfl aqu cha res rec ...
    old = int(f[11])
    f[3] = str(int(f[3]) + len(sources) - old)
    f[11] = str(len(sources))
    cnt[2] = "  ".join(f)
    (project / "object.cnt").write_text("\n".join(cnt) + "\n")
    rec = ["recall.rec: daily point sources (Basin Ledger)", "NUM     NAME    TYPE     FNAME"]
    for i, (name, ch) in enumerate(sources, 1):
        rec.append(f"   {i}    {name}      1    {name}.rec")
        lines = [f"{name}.rec: daily point source, zero in the file; set daily through BMI", "14",
                 "jday mo day_mo yr ob_typ ob_name flo sed orgn sedp no3 solp chla nh3 no2 "
                 "cbod dox san sil cla sag lag grv temp"]
        d = start
        while d <= end:
            lines.append(f"{d.timetuple().tm_yday} {d.month} {d.day} {d.year} sdc {ch} 0 " + " ".join(["0"] * 17))
            d += dt.timedelta(days=1)
        (project / f"{name}.rec").write_text("\n".join(lines) + "\n")
    (project / "recall.rec").write_text("\n".join(rec) + "\n")


def _no_recalls_yet(project: Path, channel: str):
    """A project with no point sources: write recall.con with one template row (at the channel's position and
    station) and name it in file.cio (the recall line, and recall.con in the connect line's ninth place, as SWAT+
    62.0.1's readcio_read.f90 reads them)."""
    lines = (project / "chandeg.con").read_text().splitlines()
    head = lines[1].split()
    row = next(ln.split() for ln in lines[2:] if ln.split() and ln.split()[1] == channel)
    get = dict(zip(head, row))
    (project / "recall.con").write_text(
        "recall.con: point sources (Basin Ledger)\n"
        "      NUMB    NAME    GIS_ID    AREA_HA      LAT     LONG       ELEV    RECALL     WST   CONSTIT  OVERFLOW  "
        "RULESET OUT_TOT OBTYP_OUT1 HTYPNO_OUT1  HTYP_OUT1  FRAC_OUT1\n"
        f"         1  template      1        0.000    {get['lat']}   {get['lon']}     {get['elev']}      1    "
        f"{get['wst']}         0         0        0       1         sdc          {get['id']}        tot     1.0000\n")
    cio = (project / "file.cio").read_text().splitlines()
    for i, ln in enumerate(cio):
        f = ln.split()
        if f[:1] == ["connect"]:
            while len(f) < 10:
                f.append("null")
            f[9] = "recall.con"
            cio[i] = "  ".join(f)
        elif f[:1] == ["recall"]:
            cio[i] = "recall            recall.rec"
    (project / "file.cio").write_text("\n".join(cio) + "\n")


def aquifer_pass_through(project: Path):
    """SWAT+'s aquifers keep the recharge they get and return none of it (no return flow, revap or deep recharge),
    for two-way coupling, where MODFLOW 6 handles groundwater. SWAT+ 62.0.1 (aqu_1d_control.f90) returns flow when
    the depth to water is at or above flo_min (m) and allows revap above revap_min (m), so both are set to a depth the
    water table can't reach (-1e6 m); rchg_dp is zero.

    Their soluble P too: SWAT+ 62.0.1 sends an aquifer's channels soluble P of the day's recharge volume times the
    aquifer's concentration (``aqu_d%minp = hin%flo * aqu_dat%minp / 1000``, ``hd(1)%solp``), whether or not any
    water flows out, so with return flow off a project with ``sol_p`` above zero (aquifer.aqu, ``aqu_dat%minp``) put
    phosphorus in its channels with no water, and in two-way mode MODFLOW 6's base flow is where groundwater P comes
    from. ``sol_p`` is set to zero. Its nitrate leaves only with return flow and seepage (``conc_no3 * flo``), so with
    both off its recharge nitrate stays in the aquifer's store, as its water does."""
    path = Path(project) / "aquifer.aqu"
    lines = path.read_text().splitlines()
    cols = lines[1].split()
    i_flo, i_rev, i_dp = cols.index("flo_min"), cols.index("revap_min"), cols.index("rchg_dp")
    i_p = AQUIFER_SOLP_COLUMN if len(cols) > AQUIFER_SOLP_COLUMN else None
    out = lines[:2]
    for l in lines[2:]:
        f = l.split()
        if not f:
            continue
        f[i_flo] = f[i_rev] = "-1000000.0"
        f[i_dp] = "0.0"
        if i_p is not None:
            f[i_p] = "0.0"
        out.append("  ".join(f))
    path.write_text("\n".join(out) + "\n")


# aquifer.aqu's soluble P column, by position as SWAT+ 62.0.1 reads it (aqu_read.f90 into aquifer_database: id,
# name, init, flo, dep_bot, dep_wt, no3, minp...), whatever the header calls it (the editor writes sol_p)
AQUIFER_SOLP_COLUMN = 7

# calibration.cal parameters that would undo aquifer_pass_through once SWAT+ applies them
PASS_THROUGH_UNDONE = ("flo_min", "revap_min", "rchg_dp", "chk", "minp_init")


def coupling_problems(project: Path, two_way: bool) -> list[str]:
    """What in a SWAT+ project (the run's copy, after its parameters are written) would make it fight MODFLOW 6.

    - a daily time step (``time.sim`` step 0), since Basin Ledger steps one day at a time;
    - two-way: the aquifers pass recharge through, and nothing in ``calibration.cal`` (applied after aquifer.aqu)
      turns their return flow, revap or deep recharge back on;
    - two-way: no channel transmission loss (``k`` in ``hyd-sed-lte.cha``): the stream-aquifer exchange is MODFLOW
      6's, and SWAT+'s own seepage would go to an aquifer that passes nothing on;
    - two-way: nothing drawn from SWAT+'s aquifers (``aquifer_withdrawals``)."""
    project = Path(project)
    out = []
    step = (project / "time.sim").read_text().splitlines()[2].split()[4]
    if int(step) != 0:
        out.append(f"time.sim: step {step}; Basin Ledger needs a daily time step (0)")
    if not two_way:
        return out
    lines = (project / "aquifer.aqu").read_text().splitlines()
    cols = lines[1].split()
    for ln in lines[2:]:
        f = ln.split()
        if f and (float(f[cols.index("flo_min")]) > -1e5 or float(f[cols.index("revap_min")]) > -1e5
                  or float(f[cols.index("rchg_dp")]) != 0.0):
            out.append(f"aquifer.aqu: {f[1]} doesn't pass recharge through (flo_min, revap_min, rchg_dp)")
            break
        if f and len(f) > AQUIFER_SOLP_COLUMN and float(f[AQUIFER_SOLP_COLUMN]) != 0.0:
            out.append(f"aquifer.aqu: {f[1]} has soluble P ({cols[AQUIFER_SOLP_COLUMN]} "
                       f"{f[AQUIFER_SOLP_COLUMN]}), which SWAT+ sends its channels with the recharge's volume and no "
                       "water; in two-way mode groundwater P is MODFLOW 6's, so set it to 0")
            break
    cio = (project / "file.cio").read_text().split()
    if "calibration.cal" in cio and (project / "calibration.cal").exists():
        names = {ln.split()[0] for ln in (project / "calibration.cal").read_text().splitlines()[3:] if ln.split()}
        bad = sorted(names & set(PASS_THROUGH_UNDONE))
        if bad:
            out.append(f"calibration.cal sets {', '.join(bad)}, which turns SWAT+'s own groundwater or channel loss "
                       "back on; MODFLOW 6 owns both in two-way mode, so leave them out")
    out += aquifer_withdrawals(project)
    cha = project / "hyd-sed-lte.cha"
    if cha.exists():
        lines = cha.read_text().splitlines()
        cols = lines[1].split()
        k = [float(ln.split()[cols.index("k")]) for ln in lines[2:] if ln.split()]
        if any(v != 0.0 for v in k):
            out.append(f"hyd-sed-lte.cha: {sum(v != 0 for v in k)} channels have transmission loss k > 0; in "
                       "two-way mode the river exchange is MODFLOW 6's, so set k to 0")
    return out


class StartWarning(UserWarning):
    """A setting or a warm-up that doesn't stop a run but can leave its start drifting."""


def start_warnings(project: Path) -> list[str]:
    """What in a SWAT+ project can leave the run's first years drifting: warnings, not refusals.

    - soils that start empty (``sw_frac`` 0 in ``soil_plant.ini``, for an initialization the HRUs use): SWAT+ then
      fills the profile from the first day, which can take years where soils are deep or drain slowly. A study
      basin's took more than 15 years, so a one-year warm-up left soil water, percolation and base flow rising
      through a 20-year run (Nov 7, 2026). Warm up long enough, or start the soils wetter."""
    project = Path(project)
    ini, hru = project / "soil_plant.ini", project / "hru-data.hru"
    if not ini.exists() or not hru.exists():
        return []
    lines = ini.read_text().splitlines()
    i = lines[1].lower().split().index("sw_frac")  # SWAT+ reads by position; some projects' headers are capitals
    frac = {f[0]: float(f[i]) for f in (ln.split() for ln in lines[2:]) if f}
    lines = hru.read_text().splitlines()
    j = lines[1].lower().split().index("soil_plant_init")
    used = [f[j] for f in (ln.split() for ln in lines[2:]) if len(f) > j]
    empty = sum(frac.get(u, 1.0) == 0.0 for u in used)
    if not empty:
        return []
    return [f"soil_plant.ini: {empty} of {len(used)} HRUs start with empty soils (sw_frac 0); SWAT+ fills them from "
            "the first day, which can take years, so check the warm-up's soil water (scenarios.warm_up reports it) or "
            "start the soils wetter"]


def records_every_year(project: Path):
    """SWAT+ makes its daily output records from the first year: ``nyskip`` (print.prt's years not summarized) set to
    0. swatplus-bmi's reads of those records (``hru_pw:percn``, ``hru_ls:lchlabp``, ``aquifer:minp``) are NaN in the
    years SWAT+ skips (command.f90 summarizes only after nyskip), and the transport's links read them every day. Only
    what SWAT+ prints and sums changes (its average-annual files cover the skipped years too), not what it simulates."""
    path = Path(project) / "print.prt"
    lines = path.read_text().splitlines()
    f = lines[2].split()
    if f and f[0] != "0":
        lines[2] = "  ".join(["0"] + f[1:])
        path.write_text("\n".join(lines) + "\n")


def no_transmission_loss(project: Path):
    """No channel transmission loss (``k`` = 0 in hyd-sed-lte.cha), for two-way coupling. SWAT+ 62.0.1
    (ch_rtmusk.f90) takes the loss out of the channel, records it as ``seep`` and sends it nowhere, so in two-way mode
    it would leave the basin unbooked; the stream-aquifer exchange is MODFLOW 6's."""
    path = Path(project) / "hyd-sed-lte.cha"
    if not path.exists():
        return
    lines = path.read_text().splitlines()
    i = lines[1].split().index("k")
    out = lines[:2]
    for ln in lines[2:]:
        f = ln.split()
        if f:
            f[i] = "0.00000"
            out.append("  ".join(f))
    path.write_text("\n".join(out) + "\n")


def aquifer_areas(project: Path) -> dict[str, float]:
    """Each aquifer's area in m2, from aquifer.con (hectares)."""
    lines = (Path(project) / "aquifer.con").read_text().splitlines()
    cols = lines[1].split()
    i, a = cols.index("name"), cols.index("area")
    return {f[i]: float(f[a]) * 1e4 for f in (l.split() for l in lines[2:]) if f}


def hru_areas(project: Path) -> dict[str, float]:
    """Each HRU's area in m2, from hru.con (hectares), as SWAT+ reads it into ``ob%area_ha``."""
    lines = (Path(project) / "hru.con").read_text().splitlines()
    cols = lines[1].split()
    i, a = cols.index("name"), cols.index("area")
    return {f[i]: float(f[a]) * 1e4 for f in (l.split() for l in lines[2:]) if f}


def aquifer_withdrawals(project: Path) -> list[str]:
    """What draws water from SWAT+'s own aquifers: an irrigate action whose source is an aquifer (``irrigate ...
    aqu <n>`` in a decision table, actions.f90) or a water allocation point of diversion of type ``aqu``
    (point_of_diver.wro, wallo_withdraw.f90). In two-way mode an aquifer only holds the recharge it has passed to
    MODFLOW 6, so water drawn from it would be pumped twice: groundwater is MODFLOW 6's, and its wells are where it is
    pumped (docs/design.md, "Large users and irrigation, MODFLOW first")."""
    project = Path(project)
    out = []
    used = set()  # decision tables the project's management schedules and land uses name
    for name in ("management.sch", "landuse.lum"):
        if (project / name).exists():
            used |= set((project / name).read_text(errors="replace").split())
    for f in sorted(project.glob("*.dtl")):
        lines = f.read_text(errors="replace").splitlines()
        table = None
        for i, ln in enumerate(lines):
            w = ln.split()
            if len(w) >= 2 and w[0].lower() == "name" and w[1].lower() == "conds":  # a table's header; its name next
                nxt = next((x.split() for x in lines[i + 1:] if x.split()), [None])
                table = nxt[0]
            elif len(w) >= 3 and w[0] == "irrigate" and w[1] == "aqu" and table in used:
                out.append(f"{f.name}: decision table {table} irrigates from SWAT+'s aquifer {w[2]}; in two-way mode "
                           "groundwater is MODFLOW 6's, so irrigation from it is pumped from MODFLOW 6's wells")
    cio = (project / "file.cio").read_text().splitlines()
    names = next((ln.split()[1:] for ln in cio if ln.startswith("water_rights")), [])
    for name in names:
        f = project / name
        if name == "null" or not f.exists():
            continue
        for ln in f.read_text(errors="replace").splitlines()[3:]:
            w = ln.split()
            if (len(w) >= 3 and w[2] == "aqu") or (len(w) >= 4 and w[3] == "aqu"):
                out.append(f"{name}: water allocation draws from a SWAT+ aquifer; in two-way mode groundwater is "
                           "MODFLOW 6's, so that withdrawal is pumped from MODFLOW 6's wells")
                break
    return out


def _half_ulp(text: str) -> float:
    """Half a unit in the last digit of a number as written ("0.0000806" -> 5e-8; "2.61659" -> 5e-6)."""
    t = text.lower().split("e")[0]
    d = len(t.split(".")[1]) if "." in t else 0
    e = int(text.lower().split("e")[1]) if "e" in text.lower() else 0
    return 0.5 * 10.0 ** (e - d)


def hru_area_precision(project: Path) -> dict[str, float]:
    """How far each HRU's area as SWAT+ uses it can be from the area Basin Ledger reads, m2, from the precision the
    project's files are written to: half the last digit of its area in hru.con, plus, where rout_unit.ele lists it,
    half the last digit of its share of the routing unit times that unit's area (its area over its share). A small
    HRU in a large unit, its share written to 7 decimals, can be off by 6e-4 of its area. Used to bound the daily pair
    of percolation sent and SWAT+'s own aquifer recharge."""
    lines = (Path(project) / "hru.con").read_text().splitlines()
    cols = lines[1].split()
    i, a, k = cols.index("name"), cols.index("area"), cols.index("id")
    rows = [f for f in (ln.split() for ln in lines[2:]) if f]
    out = {f[i]: _half_ulp(f[a]) * 1e4 for f in rows}
    area = {f[i]: float(f[a]) * 1e4 for f in rows}
    name_of = {int(f[k]): f[i] for f in rows}
    ele = Path(project) / "rout_unit.ele"
    if ele.exists():
        lines = ele.read_text().splitlines()
        cols = lines[1].split()
        t, n, fr = cols.index("obj_typ"), cols.index("obj_id"), cols.index("frac")
        for f in (ln.split() for ln in lines[2:]):
            if f and f[t] == "hru" and int(f[n]) in name_of and float(f[fr]) > 0:
                h = name_of[int(f[n])]
                out[h] += _half_ulp(f[fr]) * area[h] / float(f[fr])
    return out


def _channels(project: Path) -> list[str]:
    if not (Path(project) / "chandeg.con").exists():  # a project of HRUs only (Ames_sub1)
        return []
    lines = (Path(project) / "chandeg.con").read_text().splitlines()
    return [l.split()[1] for l in lines[2:] if l.split()]


def _period(project: Path) -> tuple[dt.date, dt.date]:
    f = (Path(project) / "time.sim").read_text().splitlines()[2].split()
    d0, y0, d1, y1 = (int(x) for x in f[:4])
    start = dt.date(y0, 1, 1) + dt.timedelta(days=max(d0, 1) - 1)
    end = dt.date(y1, 12, 31) if d1 == 0 else dt.date(y1, 1, 1) + dt.timedelta(days=d1 - 1)
    return start, end


class Swat:
    """Index lookups and float32 buffers for the variables the driver reads and sets."""

    def __init__(self, bmi):
        self.bmi = bmi
        self._buf = {}

    def has(self, var):
        return var in set(self.bmi.get_input_var_names()) | set(self.bmi.get_output_var_names())

    def names(self, var):
        return list(self.bmi.get_grid_object_names(self.bmi.get_var_grid(var)))

    def get(self, var) -> np.ndarray:
        if var not in self._buf:
            self._buf[var] = np.empty(self.bmi.get_grid_size(self.bmi.get_var_grid(var)), np.float32)
        self.bmi.get_value(var, self._buf[var])
        return self._buf[var].astype(np.float64)

    def set(self, var, values: np.ndarray):
        self.bmi.set_value(var, np.asarray(values, np.float32))
