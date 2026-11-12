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
ledger. Each day starts from its own stages. A limit set before the solve, at the day's starting heads, wasn't
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
at start-up (``sfr_problems``): a reach connected to a reach of another channel, or to one not coupled (SFR would
route between channels too, counting their water twice); a channel whose reaches have more than one top reach;
diversions, the mover and SFR's storage option (kinematic-wave routing: stored water a restart would need); and a
Manning's constant for other units than meters and days (UNIT_CONVERSION 86400). MODFLOW 6's default, 1, is for
meters and seconds, and makes a stream centimeters deep come out meters deep.

Other limits of this version:
- MODFLOW 6's river stage is fixed in its input unless the example's ``[rivers]`` section sets it from the surface
  model's channel depth (``stage = "surface_depth"``; SWAT+'s depth needs swatplus-bmi 0.2.0);
- the lag moves water by a day but keeps its volume; the last day's exchange is reported as in transit.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from ..mapping import cell_node
from .modflow6 import Modflow6


class ChannelMap:
    """Which entries of a MODFLOW 6 river, drain or SFR package belong to which SWAT+ channel. A map's rows give
    each entry's cell (``mapping.cell_node``), or, for an SFR package, its 0-based ``reach`` number."""

    def __init__(self, channels: list[str], entry_channel: np.ndarray):
        self.channels = channels
        self.entry_channel = entry_channel  # per package entry: index into channels, or -1 (not coupled)

    @classmethod
    def from_spec(cls, spec: str, folder: Path, mf6: Modflow6, package: str, swat_channels: list[str],
                  channels: list[str] | None = None):
        """``channels``, if given, fixes the channel order (and so the indices), as when rivers and drains share
        one list of coupled channels."""
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
        return cls(chans, entry)


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


def sfr_problems(mf6: Modflow6, package: str, cmap: ChannelMap) -> tuple[list[str], dict[str, int]]:
    """What stops an SFR package from being coupled (see the module's docstring), and each channel's top reach."""
    P = package.upper()
    g = lambda v: mf6.get(v, mf6.model, P)  # noqa: E731
    out = []
    unitconv = float(g("UNITCONV")[0])
    if abs(unitconv / DAY_S - 1.0) > 1e-9:
        out.append(f"{package}: Manning's constant is {unitconv:g}, which isn't for meters and days; set "
                   f"UNIT_CONVERSION {DAY_S:g} (or TIME_CONVERSION {DAY_S:g}) in its options")
    if int(g("ISTORAGE")[0]):
        out.append(f"{package}: STORAGE (kinematic-wave routing) keeps water in the reaches from day to day, which a "
                   "restart would need saved; not built")
    if int(g("IDIVERSIONS")[0]) or np.any(g("NDIV")):
        out.append(f"{package}: diversions take water out of a channel that SWAT+ doesn't see")
    if int(g("IMOVER")[0]):
        out.append(f"{package}: the mover moves water SWAT+ doesn't see")
    e = cmap.entry_channel
    name = lambda i: cmap.channels[e[i]] if e[i] >= 0 else "no channel"  # noqa: E731
    links = sfr_reach_links(mf6, package)
    for up, down in links:
        if e[up] != e[down]:
            out.append(f"{package}: reach {up} ({name(up)}) flows into reach {down} ({name(down)}); SWAT+ routes "
                       "water between channels, so SFR's reaches connect only within a channel (cut that connection)")
    fed = {down for _, down in links}
    top = {}
    for i, ch in enumerate(cmap.channels):
        tops = [int(r) for r in np.flatnonzero(e == i) if r not in fed]
        if len(tops) != 1:
            out.append(f"{package}: channel {ch} has {len(tops)} top reaches ({tops}); its reaches must be one chain "
                       "from one top reach, where the channel's water goes in")
        else:
            top[ch] = tops[0]
    return out, top


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
    channel's d (m, 0 where nothing was lowered)."""
    stage = np.asarray(stage, float).copy()
    floor = np.maximum(head, rbot)
    head = np.asarray(head, float)

    def net_loss(st, sel):
        lose = np.clip(st[sel] - floor[sel], 0.0, None)
        gain = np.clip(head[sel] - st[sel], 0.0, None)
        return float(np.sum(cond[sel] * (lose - gain)))

    lowered = np.zeros(len(available))
    for i, have in enumerate(np.maximum(available, 0.0)):
        sel = np.flatnonzero(entry_channel == i)
        if len(sel) == 0 or net_loss(stage, sel) <= have:
            continue
        lo, hi = 0.0, float(np.max(stage[sel] - rbot[sel]))
        s0 = stage.copy()
        for _ in range(iterations):
            mid = 0.5 * (lo + hi)
            s0[sel] = np.maximum(stage[sel] - mid, rbot[sel])
            lo, hi = (lo, mid) if net_loss(s0, sel) <= have else (mid, hi)
        stage[sel] = np.maximum(stage[sel] - hi, rbot[sel])
        lowered[i] = hi
    return stage, lowered


def requests(due: dict[str, float]) -> dict[str, float]:
    """Yesterday's exchange (m3 per channel, positive to the stream) as today's SWAT+ routes."""
    out = {}
    for ch, v in due.items():
        out[f"swat:point_source:gw_{ch}"] = max(v, 0.0)
        out[f"swat:place_of_use:seep_{ch}"] = max(-v, 0.0)
    return out


class StreamExchange:
    """The river (RIV) or streamflow routing (SFR) package's exchange, plus any drains coupled to the same channels
    (``drains``: package, map). Drains only take water from the aquifer (seepage at land surface) and give it to
    their channel the next day, as the rivers' gains do."""

    def __init__(self, mf6: Modflow6, package: str, cmap: ChannelMap, drains: tuple[str, ChannelMap] | None = None):
        self.mf6, self.package, self.map = mf6, package, cmap
        self.drains = drains
        if drains is not None and drains[1].channels != cmap.channels:
            raise ValueError("the drains' channels must be the rivers' channels, in the same order")
        P = package.upper()
        self.kind = package_type(mf6, package)
        if self.kind == "sfr":
            problems, self._top = sfr_problems(mf6, package, cmap)
            if problems:
                raise ValueError("the SFR package can't be coupled:\n- " + "\n- ".join(problems))
        elif self.kind == "riv":
            self._cond = mf6.get("COND", mf6.model, P).astype(float)
            self._rbot = mf6.get("RBOT", mf6.model, P).astype(float)
        else:
            raise ValueError(f"[exchange] package {package!r} is {self.kind.upper()}; the exchange takes a river (RIV) "
                             "or a streamflow routing (SFR) package")

    def set_inflow(self, available: dict[str, float]):
        """SFR: each channel's water for the day (m3) into its top reach, for today's solve; the channel's other
        reaches get none from outside. Reaches not coupled keep their own inflow."""
        P = self.package.upper()
        inflow = self.mf6.get("INFLOW", self.mf6.model, P).astype(float)
        inflow[self.map.entry_channel[: len(inflow)] >= 0] = 0.0
        for ch, r in self._top.items():
            inflow[r] = max(float(available[ch]), 0.0)
        self.mf6.set("INFLOW", P, inflow)

    def limit(self, available: dict[str, float]):
        """The supply limit for today's solve, as a callback for ``Modflow6.solve_day``: after each converged outer
        iteration it lowers, from the day's own stages, the stages of any channel whose net exchange (rivers and
        drains together, at the heads as they are then) takes more than ``available``, and says whether anything
        changed, so the solve goes on until heads and stages agree. ``available`` (m3 for the day) is the channel's
        flow without today's exchange: the exchange is taken tomorrow, when today's won't be there. Counting it let
        a channel lose on alternate days what it gained the day before (on a study basin). ``end_day`` puts the day's
        stages back after the solve."""
        P = self.package.upper()
        n = self.mf6.nbound(self.package)
        self._stage0 = self.mf6.get("STAGE", self.mf6.model, P).astype(float)
        nodes = self.mf6.package_nodes(self.package)
        e = self.map.entry_channel[:n]
        have = np.array([available[c] for c in self.map.channels])
        self.lowered = dict.fromkeys(self.map.channels, 0.0)

        drains = None
        if self.drains is not None:
            dp, dmap = self.drains
            nd = self.mf6.nbound(dp)
            D = dp.upper()
            drains = (self.mf6.package_nodes(dp), self.mf6.get("ELEV", self.mf6.model, D).astype(float)[:nd],
                      self.mf6.get("COND", self.mf6.model, D).astype(float)[:nd], dmap.entry_channel[:nd])

        def callback() -> bool:
            stage = self.mf6.get("STAGE", self.mf6.model, P).astype(float)
            h = self.mf6.head()
            extra = have.copy()
            if drains is not None:  # the channel's seepage at these heads
                dn, elev, dc, de = drains
                q = dc * np.clip(h[dn] - elev, 0.0, None)
                extra += np.bincount(de[de >= 0], weights=q[de >= 0], minlength=len(have))
            new, d = limit_to_supply(self._stage0[:n], self._rbot[:n], self._cond[:n], h[nodes], e, extra)
            self.lowered = dict(zip(self.map.channels, d.tolist()))
            if np.max(np.abs(new - stage[:n]), initial=0.0) <= STAGE_TOLERANCE_M:
                return False
            stage[:n] = new
            self.mf6.set_stage(self.package, stage)
            return True

        return callback

    def end_day(self):
        """The day's own stages back, so tomorrow starts from them (a fixed stage would otherwise stay lowered)."""
        if getattr(self, "_stage0", None) is not None:
            self.mf6.set_stage(self.package, self._stage0)
            self._stage0 = None

    def after_solve(self) -> dict[str, float]:
        """Today's exchange per channel (m3, positive to the stream), due in SWAT+ tomorrow: the rivers' net
        exchange plus the drains' seepage."""
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
        return out
