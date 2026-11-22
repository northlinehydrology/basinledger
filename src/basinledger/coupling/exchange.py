"""Two-way stream-aquifer exchange: MODFLOW 6 river leakage <-> SWAT+ channels, explicit with a one-day lag.

MODFLOW 6's river package (RIV) computes the exchange between each river cell and the aquifer from the river stage
and the head. Each day, after the solve, the exchange is summed per SWAT+ channel through a channel map (which RIV
entries belong to which channel). The next day it goes to SWAT+:

- **swatplus-bmi 0.2.0 and later:** as the channel's groundwater exchange
  (``channel_water__groundwater_volume_exchange_rate``, positive into the channel). A gain joins the channel's inflow
  before routing. A loss comes out of the day's outflow, no more than is there, and counts as SWAT+'s channel
  seepage. Reading the variable after the day gives what was actually exchanged, which is what the ledger enters.
- **Earlier versions:** a net discharge to the stream enters the channel as a SWAT+ point source, ``gw_<channel>``.
  A net leakage from the stream is taken as a SWAT+ water allocation demand, ``seep_<channel>``. A negative point
  source would not do: SWAT+ 62.0.1 records a negative point source in full even when the channel can't supply it.

SWAT+'s own aquifers are set to pass recharge through (``swatplus.aquifer_pass_through``: no return flow, revap or
deep recharge), so groundwater reaches the streams once, through MODFLOW 6. The ledger holds SWAT+'s aquifer outflow
at zero every day.

**Limited by supply** (``[exchange] limit``, on by default): a losing reach can't take more than its channel carries
that day. The limit works inside MODFLOW 6's solve, on the heads as they are: after each converged outer iteration,
any coupled channel losing more than its flow that day has its stages lowered by one depth, no entry below its
bottom, until the loss fits (``limit_to_supply``), and the solve goes on until heads and stages agree. Lowering the
stage is the smooth form of the river drying: the reach loses less, and gains more where the head is above it.
Today's loss is taken from tomorrow's flow (the one-day lag), so a small shortfall can remain; it stays in the
ledger. The limit works from each day's own stages, and its search starts where yesterday's ended (``[exchange]
limit_warm_start``, on by default): each channel's stages are first lowered by yesterday's depth, so a reach that was
losing more than it carried yesterday doesn't have to be found again. The answer is the same either way, since the
limit lowers from the day's own stages and puts them back up where the start was too low; on a study basin (1,157
channels, Jan-Sep 1991) summer days took 7-13 outer iterations instead of 35-54. A limit set before the solve, at the day's starting heads, wasn't
enough: wells and drains beside a river pull its heads down during the day. Without the limit, a study basin's headwater
channels, left nearly dry once SWAT+'s aquifers pass recharge through, lost 84 million m3 to MODFLOW 6 that they didn't
have.

**Through SFR** (``[exchange] package`` a streamflow routing package): MODFLOW 6 computes each reach's stage from
its flow (Manning's equation) and the reach's leakage at that stage, inside the solve, and a reach can't lose more
than it carries. Each coupled channel's water for the day (its flow without today's exchange, as for the limit
above) goes in as the inflow of the channel's first reach; the reaches below carry it, less what leaks and plus what
the aquifer gives. The channel's net exchange goes to SWAT+ the next day, as with RIV; SFR's own outflow from the
channel's last reach isn't entered anywhere, since SWAT+ routes the channel's water on. So SFR stands in for RIV's
fixed stage and for the stage-lowering limit, and SWAT+ stays the only model routing water between channels. Refused
at start-up (``sfr_problems``): a connection between channels anywhere but from one's last reach into another's
top reach, and a coupled channel flowing into a reach not coupled (SFR would route between channels too, counting
their water twice); a channel whose reaches aren't one chain from one top reach;
diversions, the mover where it doesn't line up (below), and SFR's storage option (kinematic-wave routing: stored
water a restart would need); and a
Manning's constant for other units than meters and days (UNIT_CONVERSION 86400). MODFLOW 6's default, 1, is for
meters and seconds, and makes a stream centimeters deep come out meters deep.

**SFR networks Basin Ledger didn't build** (an existing model's, docs/design.md): connections between channels are
allowed where they line up with SWAT+'s network, from one channel's last reach into another's top reach (or from the
model's own reaches into a channel's top reach), and are cut each day by setting the top reach's upstream fraction to
0. A confluence inside a channel is refused by name (``sfr_problems``); the study puts a SWAT+ outlet there
(``gridmaps.sfr_confluences`` lists them). A package built from SWAT+'s channel lines has no such connections and runs
as before.

**The mover** (coupling/mover.py, read from its files): a move from one of the package's reaches to another is a
connection to the rule, and so is a move into a channel's reach from another model's SFR or another SFR package (as
the model's own reaches). Lined up, it is cut each day by its value (it moves nothing), and the water stays with the
provider: down its own connections, or out of SFR as its external outflow. Refused by name: a move within a channel,
from a coupled reach to a reach no channel has or to any other package or model, and into a channel from a package
that isn't SFR (a drain's, a lake's, a well's: water SWAT+'s channel doesn't see). Moves among the model's own
entries stay the model's.

**Lakes at a fixed stage** (``[exchange] lakes``, a RIV package of their own): a lake's cells exchange with the aquifer
at the lake's surface, set by the study in the package's input. A lake couples through a channel, as the drains do (built before
swatplus-bmi 0.11.0 had a reservoir exchange; a LAK lake coupled to its reservoir is coupling/lakes.py): the study maps each lake's cells to the channel its reservoir
flows into, one of the exchange's coupled channels, and the lakes' net exchange joins that channel's the next day.
A lake can't take more than its channel has: with a RIV exchange the supply limit lowers a channel's river and lake
stages by one depth together; with SFR the lakes have a limit of their own, which doesn't count SFR's reach losses
from the same channel (a combined excess is a shortfall in the ledger). Basin Ledger never changes a lake's stage
otherwise, and puts it back after the day.

Other limits of this version:
- MODFLOW 6's river stage is fixed in its input unless the example's ``[rivers]`` section sets it from the surface
  model's channel depth (``stage = "surface_depth"``; SWAT+'s depth needs swatplus-bmi 0.2.0);
- the lag moves water by a day but keeps its volume; the last day's exchange is reported as in transit.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..mapping import cell_node
from . import mover as MV
from .modflow6 import Modflow6


class ChannelMap:
    """Which entries of a MODFLOW 6 river, drain or SFR package belong to which SWAT+ channel. A map's rows give
    each entry's cell (``mapping.cell_node``), or, for an SFR package, its 0-based ``reach`` number."""

    def __init__(self, channels: list[str], entry_channel: np.ndarray):
        self.channels = channels
        self.entry_channel = entry_channel  # per package entry: index into channels, or -1 (not coupled)

    @classmethod
    def from_spec(cls, spec: str, folder: Path, mf6: Modflow6, package: str, swat_channels: list[str],
                  channels: list[str] | None = None, first_model: str | None = None):
        """``channels``, if given, fixes the channel order (and so the indices), as when rivers and drains share
        one list of coupled channels. ``first_model``: with several groundwater models coupled, the first's name;
        then only the rows of ``mf6``'s model are read (a row without a ``model`` is the first's), and ``channels``
        keeps only the channels those rows give."""
        nodes = mf6.package_nodes(package)
        if spec.startswith("all:"):
            ch = spec[4:]
            if ch not in swat_channels:
                raise KeyError(f"exchange map: {ch!r} is not a SWAT+ channel")
            chans = list(channels) if channels else [ch]
            return cls(chans, np.full(len(nodes), chans.index(ch), dtype=int))
        entry = np.full(len(nodes), -1)
        chans: list[str] = list(channels) if channels else []
        sfr = package_type(mf6, package) == "sfr"
        with open(folder / spec, newline="") as f:
            for r in csv.DictReader(f):
                if first_model is not None and \
                        ((r.get("model") or "").strip() or first_model).upper() != mf6.model.upper():
                    continue
                ch = r["object"]
                if ch not in swat_channels:
                    raise KeyError(f"{spec}: {ch!r} is not a SWAT+ channel")
                if ch not in chans:
                    chans.append(ch)
                if (r.get("reach") or "").strip():
                    if not sfr:
                        raise ValueError(f"{spec}: reach numbers are for an SFR package; {package!r} isn't one")
                    k = int(r["reach"])
                    if not 0 <= k < len(nodes):
                        raise ValueError(f"{spec}: reach {k} isn't in {package!r} ({len(nodes)} reaches, from 0)")
                    hit = np.array([k])
                else:
                    n = cell_node(r, mf6.node, spec)
                    hit = np.flatnonzero(nodes == n)
                    if len(hit) == 0:
                        raise ValueError(f"{spec}: cell {n} has no {package} entry")
                entry[hit] = chans.index(ch)
        if first_model is not None:  # this model's channels only, in the order given
            used = [c for i, c in enumerate(chans) if (entry == i).any()]
            entry = np.array([used.index(chans[i]) if i >= 0 else -1 for i in entry], dtype=int)
            chans = used
        return cls(chans, entry)


def part_map(spec: str, folder: Path, mf6: Modflow6, package: str, swat_channels: list[str], channels: list[str],
             first_model: str) -> tuple[ChannelMap, list[str]]:
    """With several groundwater models coupled, a drains' or lakes' map in one model (``mf6``, a view): its rows
    only, onto that model's coupled channels (``channels``, in their order). Returns the map and the channels its
    rows name that aren't among them (a drain or lake whose channel's reaches are in another model)."""
    m = ChannelMap.from_spec(spec, folder, mf6, package, swat_channels, first_model=first_model)
    extra = [c for c in m.channels if c not in channels]
    idx = np.array([channels.index(m.channels[i]) if i >= 0 and m.channels[i] in channels else -1
                    for i in m.entry_channel], dtype=int)
    return ChannelMap(list(channels), idx), extra


def package_type(mf6: Modflow6, package: str) -> str:
    """A MODFLOW 6 package's type, lower case: "riv", "drn", "sfr", ..."""
    return str(mf6.get("PACKAGE_TYPE", mf6.model, package.upper())[0]).strip().lower()


DAY_S = 86400.0  # Manning's constant for meters and days (MODFLOW 6's UNIT_CONVERSION)


def sfr_reach_links(mf6: Modflow6, package: str) -> list[tuple[int, int]]:
    """Each (upstream, downstream) pair of connected reaches in an SFR package, 0-based."""
    P = package.upper()
    ia = mf6.get("IA", mf6.model, P).astype(int) - 1
    ja = mf6.get("JA", mf6.model, P).astype(int) - 1
    idir = mf6.get("IDIR", mf6.model, P).astype(int)  # per connection: 1 upstream of the reach, -1 downstream
    return [(i, int(ja[k])) for i in range(len(ia) - 1) for k in range(ia[i], ia[i + 1]) if idir[k] == -1]


SFR_MODES = ("route", "route_or_stage", "stage")  # [exchange] sfr: A everywhere, A where a channel lines up, B everywhere


@dataclass
class SfrLayout:
    """How an SFR package lines up with SWAT+'s channels (``sfr_problems``)."""
    problems: list[str]                    # refused whatever the mode
    channel_problems: dict[str, list[str]]  # what keeps a channel from routing (A); it can run at a fixed stage (B)
    top: dict[str, int]                    # each channel's top reach, where it lines up
    last: dict[str, int]                   # each channel's last reach, where it lines up
    cut: list[int]                         # reaches whose upstream fraction is set to 0 each day
    moves_cut: list = field(default_factory=list)  # the mover's moves set to move nothing each day (coupling/mover.py)
    movers: list = field(default_factory=list)     # the simulation's movers, read from their files
    drains_in: list = field(default_factory=list)  # (a drain entry, the channel its move fed): seepage to that channel


def sfr_problems(mf6: Modflow6, package: str, cmap: ChannelMap,
                 others: dict[tuple[str, str], ChannelMap] | None = None) -> SfrLayout:
    """How an SFR package lines up with SWAT+'s channels: what is refused, what keeps each channel from routing its
    own water (A), each channel's top and last reach, and the reaches to cut from what flows into them at run time
    (their upstream fraction set to 0).

    The line-up rule, the same for a package Basin Ledger built from SWAT+'s channel lines and for one a study
    brought (docs/design.md, "Extending an existing MODFLOW 6 model"): each reach maps to one channel, or none and
    stays the model's; each channel's reaches form one connected path from one top reach; and every confluence falls
    on a channel boundary: a connection into a coupled channel goes into its top reach, from another channel's last
    reach or from a reach no channel has. Such a connection is cut each day, so SWAT+ alone carries water from one
    channel to the next. A channel that breaks the rule can't route (A), but can run at a fixed stage (B,
    ``[exchange] sfr``); its incoming connections are cut wherever they are. Refused whatever the mode: the package's
    settings below, and a coupled channel flowing into a reach no channel has (SFR would carry SWAT+'s water on).

    ``others``: with several groundwater models coupled, the other models' SFR maps, by (model, package), lower case.
    A coupled reach giving water through the mover to another model's coupled reach is then lined up as a
    connection between channels is (from its channel's last reach; the other side checks it reaches a top reach)
    and cut; into the same channel there, it is a move within a channel, refused."""
    P = package.upper()
    g = lambda v: mf6.get(v, mf6.model, P)  # noqa: E731
    out = []
    # Manning's constant in the model's own units: 86400 for meters and days, 1.486 x 86400 for feet and days (the
    # 1.486 written for feet is 1.4859 to five figures, so a thousandth is allowed)
    unitconv = float(mf6.get_model_units("UNITCONV", mf6.model, P)[0])
    want = DAY_S * mf6.days / mf6.meters ** (1.0 / 3.0)
    if abs(unitconv / want - 1.0) > 1e-3:
        out.append(f"{package}: Manning's constant is {unitconv:g}, which isn't for {mf6.length_units} and "
                   f"{mf6.time_units}; set UNIT_CONVERSION {want:g} in its options")
    if int(g("ISTORAGE")[0]):
        out.append(f"{package}: STORAGE (kinematic-wave routing) keeps water in the reaches from day to day, which a "
                   "restart would need saved; not built")
    if int(g("IDIVERSIONS")[0]) or np.any(g("NDIV")):
        out.append(f"{package}: diversions take water out of a channel that SWAT+ doesn't see")
    e = cmap.entry_channel
    name = lambda i: cmap.channels[e[i]] if e[i] >= 0 else "no channel"  # noqa: E731
    links = sfr_reach_links(mf6, package)
    movers, moved, moves_cut = [], [], []  # moved: (up, down, move) between this package's reaches, cut each day
    drains_in = []  # (drain entry, channel): a drain's move into a coupled reach, cut, its water the channel's
    leaving = []  # (up, move): from a coupled reach into another coupled model's channel
    if int(g("IMOVER")[0]):
        movers = MV.movers(mf6.ws)
        me = (mf6.model.lower(), package.lower())
        for m in sorted({m for mv in movers for m in mv.moves}, key=lambda m: (m.mover, m.period, m.index)):
            p_mine = (m.provider.model, m.provider.package) == me
            r_mine = (m.receiver.model, m.receiver.package) == me
            u = m.provider.id if p_mine else None
            d = m.receiver.id if r_mine else None
            via = f"through the mover ({m.mover.split('/')[0]}, period {m.period}'s move {m.index + 1})"
            if u is not None and d is not None:
                if e[u] < 0 and e[d] < 0:
                    continue  # the model's own reaches
                if e[u] == e[d]:
                    out.append(f"{package}: reach {u} gives water to reach {d} in channel {name(u)} {via}; a "
                               "mover within a channel isn't built")
                elif e[d] < 0:
                    out.append(f"{package}: reach {u} ({name(u)}) gives water to reach {d} (no channel) {via}; SWAT+ "
                               "routes water between channels, so SFR would carry the channel's water on")
                else:
                    moved.append((u, d, m))
            elif u is not None and e[u] >= 0:
                other = (others or {}).get((m.receiver.model, m.receiver.package))
                there = other.entry_channel[m.receiver.id] if other is not None else -1
                if there >= 0 and other.channels[there] == name(u):
                    out.append(f"{package}: reach {u} gives water to {m.receiver} in channel {name(u)} {via}; a mover "
                               "within a channel isn't built")
                elif there >= 0:
                    leaving.append((u, m))
                elif other is not None:
                    out.append(f"{package}: reach {u} ({name(u)}) gives water to {m.receiver} (no channel) {via}; "
                               "SWAT+ routes water between channels, so SFR would carry the channel's water on")
                else:
                    out.append(f"{package}: reach {u} ({name(u)}) gives water to {m.receiver} {via}; SWAT+ routes the "
                               "channel's water, so the mover can't take it elsewhere")
            elif d is not None and e[d] >= 0:
                kind = _type_of(mf6, m.provider)
                if kind == "sfr":  # another model's reaches, or another SFR package's, as the model's own reaches
                    moved.append((None, d, m))
                elif kind == "drn":  # seepage at land surface: cut, and the drain's water the channel's next day
                    moves_cut.append(m)
                    drains_in.append((m.provider, name(d)))
                else:
                    out.append(f"{package}: {m.provider} ({kind.upper() or 'a package'}) gives water to reach {d} "
                               f"({name(d)}) {via}: water SWAT+'s channel doesn't see; not built")
    within = [(u, d) for u, d in links if e[u] == e[d] and e[u] >= 0]
    fed = {d for _, d in within}
    bad: dict[str, list[str]] = {ch: [] for ch in cmap.channels}
    top, last = {}, {}
    for i, ch in enumerate(cmap.channels):
        mine = [int(r) for r in np.flatnonzero(e == i)]
        tops = [r for r in mine if r not in fed]
        if len(tops) != 1:
            bad[ch].append(f"{package}: channel {ch} has {len(tops)} top reaches ({tops}); its reaches must be one "
                           "chain from one top reach, where the channel's water goes in")
            continue
        down = {}
        for u, d in within:
            if e[u] == i:
                down.setdefault(u, []).append(d)
        split = [u for u, ds in down.items() if len(ds) > 1]
        if split:
            bad[ch].append(f"{package}: channel {ch} branches below reach {split[0]} (into {down[split[0]]}); its "
                           "reaches must be one chain")
            continue
        top[ch] = tops[0]
        last[ch] = [r for r in mine if r not in down][0]
    cut = set()
    joins = [(up, down, None) for up, down in links] + moved
    for up, down, m in joins:
        if m is None and e[up] == e[down]:
            continue
        if m is None and e[down] < 0:
            out.append(f"{package}: reach {up} ({name(up)}) flows into reach {down} ({name(down)}); SWAT+ routes "
                       "water between channels, so SFR would carry the channel's water on (cut that connection, or "
                       "give the reaches below a channel)")
            continue
        if m is None:
            cut.add(int(down))
        else:
            moves_cut.append(m)
        how = "flows into" if m is None else "gives water through the mover to"
        src = f"reach {up} ({name(up)})" if up is not None else str(m.provider)
        if top.get(name(down)) != down:
            bad[name(down)].append(
                f"{package}: {src} {how} reach {down} ({name(down)}), inside channel "
                f"{name(down)}; every confluence must fall on a channel boundary (a SWAT+ outlet at reach {down}'s top)")
        if up is not None and e[up] >= 0 and last.get(name(up)) != up:
            bad[name(up)].append(f"{package}: {src} {how} reach {down} ({name(down)}) from "
                                 f"inside channel {name(up)}; a channel leaves SFR only from its last reach")
    for up, m in leaving:
        moves_cut.append(m)
        if last.get(name(up)) != up:
            bad[name(up)].append(f"{package}: reach {up} ({name(up)}) gives water to {m.receiver} through the mover "
                                 f"from inside channel {name(up)}; a channel leaves SFR only from its last reach")
    bad = {ch: v for ch, v in bad.items() if v}
    for ch in bad:
        top.pop(ch, None)
        last.pop(ch, None)
    return SfrLayout(out, bad, top, last, sorted(cut), sorted(moves_cut, key=lambda m: (m.mover, m.period, m.index)),
                     movers, drains_in)


def _type_of(mf6: Modflow6, entry) -> str:
    """The type of a mover's provider or receiver package, from libmf6 ("" where it isn't found)."""
    try:
        return str(mf6.api.get_value(f"{entry.model.upper()}/{entry.package.upper()}/PACKAGE_TYPE")[0]).strip().lower()
    except Exception:  # noqa: BLE001
        return ""


def stage_entries(mf6: Modflow6, package: str, cmap: ChannelMap, channel: str, cut: list[int]) -> dict[int, float]:
    """For a channel at a fixed stage (B): where its water goes in, and each place's share of it. The places are the
    channel's reaches nothing in the channel flows into, and its reaches cut from water coming from outside it; each
    takes the share of the channel's reach length from it down to the next such place (SFR's reach lengths)."""
    i = cmap.channels.index(channel)
    e = cmap.entry_channel
    mine = [int(r) for r in np.flatnonzero(e == i)]
    within = [(u, d) for u, d in sfr_reach_links(mf6, package) if e[u] == i and e[d] == i]
    fed = {d for _, d in within}
    entries = [r for r in mine if r not in fed or r in cut]
    down: dict[int, list[int]] = {}
    for u, d in within:
        down.setdefault(u, []).append(d)
    length = mf6.get("LENGTH", mf6.model, package.upper()).astype(float)
    owner: dict[int, int] = {}
    for r0 in entries:
        stack = [r0]
        while stack:
            r = stack.pop()
            if r in owner or (r != r0 and r in entries):
                continue
            owner[r] = r0
            stack.extend(down.get(r, []))
    total = sum(length[r] for r in owner)
    share = {r0: sum(length[r] for r, o in owner.items() if o == r0) / total for r0 in entries}
    return share


STAGE_TOLERANCE_M = 1e-3  # the supply limit's stages agree with the heads when no stage moves more than this


def limit_to_supply(stage: np.ndarray, rbot: np.ndarray, cond: np.ndarray, head: np.ndarray, entry_channel: np.ndarray,
                    available: np.ndarray, iterations: int = 60) -> tuple[np.ndarray, np.ndarray]:
    """Stages lowered so that no channel's net exchange takes more than it has, at the given heads.

    Per entry, MODFLOW 6's river package exchanges ``cond * (stage - max(head, rbot))`` with the aquifer: into it
    when positive (a loss), out of it when negative (a gain). A channel's net loss is the sum over its entries
    (``entry_channel``, -1 for entries not coupled). Where it is more than ``available`` (m3 for the day; a channel's
    other gains, such as seepage, can be added to it), the channel's stages are lowered by one depth d,
    ``max(stage - d, rbot)``, with d found by bisection: the losses shrink and the gains grow as d grows, and at
    d = max(stage - rbot) the channel loses nothing, so a d always exists. Returns the new stages and each
    channel's d (m, 0 where nothing was lowered). The channels are bisected together, each on its own interval."""
    stage = np.asarray(stage, float).copy()
    rbot, cond = np.asarray(rbot, float), np.asarray(cond, float)
    floor = np.maximum(head, rbot)
    head = np.asarray(head, float)
    e = np.asarray(entry_channel)
    n = len(available)
    have = np.maximum(np.asarray(available, float), 0.0)

    def net_loss(st, idx, ch):  # each channel's net loss over the entries idx (their channels ch)
        q = cond[idx] * (np.clip(st - floor[idx], 0.0, None) - np.clip(head[idx] - st, 0.0, None))
        return np.bincount(ch, weights=q, minlength=n)

    coupled = np.flatnonzero(e >= 0)
    over = net_loss(stage[coupled], coupled, e[coupled]) > have
    lowered = np.zeros(n)
    if not over.any():
        return stage, lowered
    idx = coupled[over[e[coupled]]]
    ch = e[idx]
    lo, hi = np.zeros(n), np.full(n, -np.inf)
    np.maximum.at(hi, ch, stage[idx] - rbot[idx])
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        ok = net_loss(np.maximum(stage[idx] - mid[ch], rbot[idx]), idx, ch) <= have
        lo, hi = np.where(ok, lo, mid), np.where(ok, mid, hi)
    stage[idx] = np.maximum(stage[idx] - hi[ch], rbot[idx])
    lowered[over] = hi[over]
    return stage, lowered


def requests(due: dict[str, float]) -> dict[str, float]:
    """Yesterday's exchange (m3 per channel, positive to the stream) as today's SWAT+ routes."""
    out = {}
    for ch, v in due.items():
        out[f"swat:point_source:gw_{ch}"] = max(v, 0.0)
        out[f"swat:place_of_use:seep_{ch}"] = max(-v, 0.0)
    return out


def check_lake_package(mf6: Modflow6, lakes: str, package: str, drains: str | None = None):
    """The lakes' package is a river (RIV) package of their own: not the exchange's, nor the drains'."""
    if lakes.lower() == package.lower():
        raise ValueError(f"[exchange] lakes: {lakes!r} is the exchange's own package; lakes need a river (RIV) "
                         "package of their own")
    if drains is not None and lakes.lower() == drains.lower():
        raise ValueError(f"[exchange] lakes: {lakes!r} is the drains' package; lakes need a river (RIV) package of "
                         "their own")
    if package_type(mf6, lakes) != "riv":
        raise ValueError(f"[exchange] lakes: {lakes!r} is {package_type(mf6, lakes).upper()}; lakes at a fixed stage "
                         "are a river (RIV) package")


def lake_map(spec: str, folder: Path, mf6: Modflow6, package: str, swat_channels: list[str],
             coupled: list[str], first_model: str | None = None) -> ChannelMap:
    """The lakes' map, onto the exchange's coupled channels and in their order: a lake couples through its outlet
    channel, which must be one of them. A lake may map to some of them only. ``first_model``: with several
    groundwater models coupled, the rows of ``mf6``'s model only, onto its own channels (``coupled``)."""
    if spec.startswith("all:") and spec[4:] in swat_channels and spec[4:] not in coupled:
        extra = [spec[4:]]
    elif first_model is not None:
        lmap, extra = part_map(spec, folder, mf6, package, swat_channels, coupled, first_model)
    else:
        lmap = ChannelMap.from_spec(spec, folder, mf6, package, swat_channels, channels=coupled)
        extra = lmap.channels[len(coupled):]
    if extra:
        raise ValueError(f"[exchange] lakes: {', '.join(extra)} isn't a channel the exchange couples ({coupled}); a "
                         "lake couples through its outlet channel, which must be one of the rivers' or drains'")
    return lmap


class StreamExchange:
    """The river (RIV) or streamflow routing (SFR) package's exchange, plus any drains coupled to the same channels
    (``drains``: package, map) and any lakes at a fixed stage (``lakes``: a RIV package, map). Drains only take
    water from the aquifer (seepage at land surface) and give it to their channel the next day, as the rivers'
    gains do. Lakes gain or lose, and their net exchange joins their outlet channel's the next day."""

    def __init__(self, mf6: Modflow6, package: str, cmap: ChannelMap, drains: tuple[str, ChannelMap] | None = None,
                 *, sfr: str = "route", others: dict[tuple[str, str], ChannelMap] | None = None,
                 lakes: tuple[str, ChannelMap] | None = None, stage_tolerance: float = STAGE_TOLERANCE_M,
                 warm_start: bool = True):
        self.mf6, self.package, self.map = mf6, package, cmap
        self.drains = drains
        self.lakes = lakes
        if not stage_tolerance > 0:
            raise ValueError(f"[exchange] limit_tolerance must be a depth above 0 m, not {stage_tolerance!r}")
        self.stage_tolerance = float(stage_tolerance)
        self.warm_start = bool(warm_start)
        self.last_lowered = None  # each channel's depth (m) at the end of the last day's limit, where tomorrow's starts
        if drains is not None and drains[1].channels != cmap.channels:
            raise ValueError("the drains' channels must be the rivers' channels, in the same order")
        if lakes is not None:
            lp, lmap = lakes
            check_lake_package(mf6, lp, package, drains[0] if drains is not None else None)
            if lmap.channels != cmap.channels:
                raise ValueError("the lakes' channels must be the rivers' channels, in the same order")
            L = lp.upper()
            self._lake_cond = mf6.get("COND", mf6.model, L).astype(float)
            self._lake_rbot = mf6.get("RBOT", mf6.model, L).astype(float)
        P = package.upper()
        self.kind = package_type(mf6, package)
        self.modes: dict[str, tuple[str, list[str]]] = {}  # SFR: channel -> ("route" or "stage", why)
        if self.kind == "sfr":
            if sfr not in SFR_MODES:
                raise ValueError(f"[exchange] sfr = {sfr!r}: one of {SFR_MODES}")
            lay = sfr_problems(mf6, package, cmap, others)
            problems = list(lay.problems)
            if sfr == "route":
                problems += [m for ms in lay.channel_problems.values() for m in ms]
            if problems:
                raise ValueError("the SFR package can't be coupled:\n- " + "\n- ".join(problems))
            self._cut = lay.cut
            self._movers, self._moves_cut = lay.movers, lay.moves_cut
            self._drains_in = self._drains_by_mover(lay.drains_in)
            for ch in cmap.channels:
                if sfr == "stage" or ch in lay.channel_problems:
                    why = lay.channel_problems.get(ch, ["[exchange] sfr = \"stage\""])
                    self.modes[ch] = ("stage", why)
                else:
                    self.modes[ch] = ("route", [])
            self._top = {ch: r for ch, r in lay.top.items() if self.modes[ch][0] == "route"}
            self._stage = {ch: stage_entries(mf6, package, cmap, ch, lay.cut) for ch, (m, _) in self.modes.items()
                           if m == "stage"}
            e = cmap.entry_channel
            self._stage_reaches = np.flatnonzero(np.isin(e, [cmap.channels.index(c) for c in self._stage]))
            self._strtop = mf6.get("STRTOP", mf6.model, P).astype(float)
        elif self.kind == "riv":
            self._cond = mf6.get("COND", mf6.model, P).astype(float)
            self._rbot = mf6.get("RBOT", mf6.model, P).astype(float)
        else:
            raise ValueError(f"[exchange] package {package!r} is {self.kind.upper()}; the exchange takes a river (RIV) "
                             "or a streamflow routing (SFR) package")

    def _drains_by_mover(self, drains_in) -> dict[str, list]:
        """Each channel's drain entries whose moves into its reaches are cut: (model, package, entry), upper case.
        An entry the exchange's drains already map to the same channel is counted there; to another, refused."""
        out: dict[str, list] = {}
        for entry, ch in drains_in:
            key = (entry.model.upper(), entry.package.upper(), int(entry.id))
            if self.drains is not None and key[:2] == (self.mf6.model.upper(), self.drains[0].upper()):
                mapped = int(self.drains[1].entry_channel[key[2]])
                if mapped >= 0 and self.drains[1].channels[mapped] == ch:
                    continue
                if mapped >= 0:
                    raise ValueError(f"{entry} seeps to channel {self.drains[1].channels[mapped]} by the drains' map "
                                     f"and gives water to channel {ch} through the mover; one channel takes it")
            out.setdefault(ch, []).append(key)
        return out

    def set_inflow(self, available: dict[str, float], depth: dict[str, float] | None = None):
        """SFR: each channel's water for the day (m3) into its top reach, for today's solve; the channel's other
        reaches get none from outside. Reaches not coupled keep their own inflow. A reach fed by another channel's
        last reach, or by the model's own reaches, has its upstream fraction set to 0: that water leaves SFR there
        (booked as SFR's external outflow), and SWAT+ carries it on.

        A channel at a fixed stage (B, ``modes``): its reaches take status simple, each one's stage its streambed top
        plus SWAT+'s channel depth for the day (``depth``, m), and its water goes in at its entries
        (``stage_entries``), each by its share. A reach at a fixed stage still routes what reaches it and loses no
        more than that."""
        P = self.package.upper()
        inflow = self.mf6.get("INFLOW", self.mf6.model, P).astype(float)
        inflow[self.map.entry_channel[: len(inflow)] >= 0] = 0.0
        for ch, r in self._top.items():
            inflow[r] = max(float(available[ch]), 0.0)
        for ch, shares in self._stage.items():
            for r, f in shares.items():
                inflow[r] = max(float(available[ch]), 0.0) * f
        self.mf6.set("INFLOW", P, inflow)
        if self._cut:  # what flows into a channel from another, or from the model's own reaches, leaves SFR there
            ustrf = self.mf6.get("USTRF", self.mf6.model, P).astype(float)
            ustrf[self._cut] = 0.0
            self.mf6.set("USTRF", P, ustrf)
        if self._moves_cut:  # and what the mover would carry into a channel, or out of one into another
            MV.cut(self.mf6, self._movers, self._moves_cut)
        if self._stage:
            if depth is None:
                raise ValueError("SFR channels at a fixed stage need SWAT+'s channel depth for the day")
            e = self.map.entry_channel
            status = self.mf6.get("IBOUNDPAK", self.mf6.model, P)
            stage = self.mf6.get("SSTAGE", self.mf6.model, P).astype(float)
            r = self._stage_reaches
            status[r] = -1
            stage[r] = self._strtop[r] + np.array([max(float(depth[self.map.channels[e[k]]]), 0.0) for k in r])
            self.mf6.set("IBOUNDPAK", P, status)
            self.mf6.set("SSTAGE", P, stage)

    def limit(self, available: dict[str, float]):
        """The supply limit for today's solve, as a callback for ``Modflow6.solve_day``: after each converged outer
        iteration it lowers, from the day's own stages, the stages of any channel whose net exchange (rivers, drains
        and lakes together, at the heads as they are then) takes more than ``available``, and says whether anything
        changed, so the solve goes on until heads and stages agree. ``available`` (m3 for the day) is the channel's
        flow without today's exchange: the exchange is taken tomorrow, when today's won't be there. Counting it let
        a channel lose on alternate days what it gained the day before (on a study basin). A channel's river and lake
        stages are lowered by the same depth. With an SFR exchange only the lakes' stages are limited here (SFR
        limits its own reaches), and SFR's losses from the same channel aren't counted. ``end_day`` puts the day's
        stages back after the solve."""
        have = np.array([available[c] for c in self.map.channels])
        parts = []  # (package, rbot, cond, map): the rivers' (RIV only), then the lakes'
        if self.kind == "riv":
            parts.append((self.package, self._rbot, self._cond, self.map))
        if self.lakes is not None:
            parts.append((self.lakes[0], self._lake_rbot, self._lake_cond, self.lakes[1]))
        self._stage0 = {}
        day = []  # (package, entries, nodes, rbot, cond, entry_channel), each package's entries in turn
        for pkg, rbot, cond, cmap in parts:
            n = self.mf6.nbound(pkg)
            self._stage0[pkg] = self.mf6.get("STAGE", self.mf6.model, pkg.upper()).astype(float)
            day.append((pkg, n, self.mf6.package_nodes(pkg), rbot[:n], cond[:n], cmap.entry_channel[:n]))
        # the packages' entries together, so a channel's rivers and lakes are lowered by one depth
        stage0 = np.concatenate([self._stage0[p][:n] for p, n, *_ in day])
        nodes, rbot, cond, e = (np.concatenate([x[i] for x in day]) for i in (2, 3, 4, 5))
        lake_channels = (np.isin(np.arange(len(have)), self.lakes[1].entry_channel) if self.lakes is not None
                         else np.zeros(len(have), bool))
        self.lowered = dict.fromkeys(self.map.channels, 0.0)
        self.lake_lowered = dict.fromkeys(self.map.channels, 0.0)
        last = self.last_lowered
        if self.warm_start and last is not None and last.any():
            # start the day's search where yesterday's ended: each channel's stages lowered by yesterday's depth.
            # The limit still works from the day's own stages, so where this is too low it puts them back up.
            k = e >= 0
            start = stage0.copy()
            start[k] = np.maximum(stage0[k] - last[e[k]], rbot[k])
            at = 0
            for pkg, n, *_ in day:
                st = self.mf6.get("STAGE", self.mf6.model, pkg.upper()).astype(float)
                st[:n] = start[at: at + n]
                at += n
                self.mf6.set_stage(pkg, st)
        self.last_lowered = np.zeros(len(have))

        drains = None
        if self.drains is not None:
            dp, dmap = self.drains
            nd = self.mf6.nbound(dp)
            D = dp.upper()
            drains = (self.mf6.package_nodes(dp), self.mf6.get("ELEV", self.mf6.model, D).astype(float)[:nd],
                      self.mf6.get("COND", self.mf6.model, D).astype(float)[:nd], dmap.entry_channel[:nd])

        def callback() -> bool:
            stages = [self.mf6.get("STAGE", self.mf6.model, p.upper()).astype(float) for p, *_ in day]
            h = self.mf6.head()
            extra = have.copy()
            if drains is not None:  # the channel's seepage at these heads
                dn, elev, dc, de = drains
                q = dc * np.clip(h[dn] - elev, 0.0, None)
                extra += np.bincount(de[de >= 0], weights=q[de >= 0], minlength=len(have))
            new, d = limit_to_supply(stage0, rbot, cond, h[nodes], e, extra)
            self.last_lowered = d
            if self.kind == "riv":
                self.lowered = dict(zip(self.map.channels, d.tolist()))
            self.lake_lowered = dict(zip(self.map.channels, np.where(lake_channels, d, 0.0).tolist()))
            now = np.concatenate([st[:n] for st, (_, n, *_) in zip(stages, day)])
            if np.max(np.abs(new - now), initial=0.0) <= self.stage_tolerance:
                return False
            at = 0
            for st, (pkg, n, *_) in zip(stages, day):
                st[:n] = new[at: at + n]
                at += n
                self.mf6.set_stage(pkg, st)
            return True

        return callback

    def limit_state(self) -> list[float] | None:
        """What a restart needs of the limit: each channel's depth at the end of the last day (where the next day's
        search starts)."""
        return None if self.last_lowered is None else [float(x) for x in self.last_lowered]

    def set_limit_state(self, state: list[float] | None):
        self.last_lowered = None if state is None else np.asarray(state, float)

    def end_day(self):
        """The day's own stages back, so tomorrow starts from them (a fixed stage would otherwise stay lowered)."""
        for pkg, stage in (getattr(self, "_stage0", None) or {}).items():
            self.mf6.set_stage(pkg, stage)
        self._stage0 = None

    def after_solve(self) -> dict[str, float]:
        """Today's exchange per channel (m3, positive to the stream), due in SWAT+ tomorrow: the rivers' net
        exchange plus the drains' seepage and the lakes' net exchange."""
        q = self.mf6.simvals(self.package)  # m3 d-1, positive into the aquifer
        e = self.map.entry_channel[: len(q)]
        out = {ch: -float(q[e == i].sum()) for i, ch in enumerate(self.map.channels)}
        self.seepage = {ch: 0.0 for ch in out}
        if self.drains is not None:
            pkg, dmap = self.drains
            qd = self.mf6.simvals(pkg)
            ed = dmap.entry_channel[: len(qd)]
            for i, ch in enumerate(dmap.channels):
                self.seepage[ch] = -float(qd[ed == i].sum())
                out[ch] += self.seepage[ch]
        self.lake_flow = {ch: 0.0 for ch in out}
        if self.lakes is not None:
            pkg, lmap = self.lakes
            ql = self.mf6.simvals(pkg)
            el = lmap.entry_channel[: len(ql)]
            for i, ch in enumerate(lmap.channels):
                self.lake_flow[ch] = -float(ql[el == i].sum())
                out[ch] += self.lake_flow[ch]
        self.drains_moved = {ch: 0.0 for ch in out}  # drains whose moves into the channel's reaches are cut
        for ch, entries in getattr(self, "_drains_in", {}).items():
            for model, package, i in entries:
                q = -float(self.mf6.get("SIMVALS", model, package)[i])  # out of the aquifer, positive
                self.drains_moved[ch] += q
            out[ch] += self.drains_moved[ch]
        return out


class StreamExchanges:
    """The stream exchange across several coupled groundwater models (``[exchange]`` with a map whose rows name
    their model): one ``StreamExchange`` per model, each with its own channels; a channel's reaches are in one model.
    What the driver asks of one exchange is asked of each, and the channels' results put together."""

    def __init__(self, parts: list[StreamExchange]):
        kinds = {p.kind for p in parts}
        if len(kinds) != 1:
            raise ValueError(f"[exchange]: the models' packages are {sorted(kinds)}; one kind, RIV or SFR, for all")
        self.parts = parts
        self.kind = parts[0].kind
        self.package = parts[0].package
        # drains and lakes at a fixed stage: each model's own, on its own channels (any part's, for the driver's
        # questions whether there are any; the values are put together from every part's)
        self.drains = next((p.drains for p in parts if p.drains is not None), None)
        self.lakes = next((p.lakes for p in parts if p.lakes is not None), None)
        self.modes = {c: m for p in parts for c, m in p.modes.items()}
        self._stage = {c: x for p in parts for c, x in getattr(p, "_stage", {}).items()}
        self.channels = [c for p in parts for c in p.map.channels]

    def set_inflow(self, available: dict[str, float], depth: dict[str, float] | None = None):
        for p in self.parts:
            p.set_inflow(available, depth)

    def limit_state(self) -> list | None:
        return [p.limit_state() for p in self.parts]

    def set_limit_state(self, state: list | None):
        for p, x in zip(self.parts, state or [None] * len(self.parts)):
            p.set_limit_state(x)

    def limit(self, available: dict[str, float]):
        callbacks = [p.limit(available) for p in self.parts]

        def callback() -> bool:
            changed = [cb() for cb in callbacks]  # every model's, not stopping at the first that changed
            return any(changed)

        return callback

    @property
    def lowered(self) -> dict[str, float]:
        return {c: d for p in self.parts for c, d in p.lowered.items()}

    @property
    def lake_lowered(self) -> dict[str, float]:
        return {c: getattr(p, "lake_lowered", {}).get(c, 0.0) for p in self.parts for c in p.map.channels}

    def end_day(self):
        for p in self.parts:
            p.end_day()

    def after_solve(self) -> dict[str, float]:
        out, self.seepage, self.lake_flow, self.drains_moved = {}, {}, {}, {}
        for p in self.parts:
            out.update(p.after_solve())
            self.seepage.update(getattr(p, "seepage", None) or dict.fromkeys(p.map.channels, 0.0))
            self.lake_flow.update(getattr(p, "lake_flow", None) or dict.fromkeys(p.map.channels, 0.0))
            self.drains_moved.update(getattr(p, "drains_moved", None) or dict.fromkeys(p.map.channels, 0.0))
        return out
