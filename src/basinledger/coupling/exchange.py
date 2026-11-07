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
enough: wells and drains beside a river pull its heads down during the day. Without the limit, the Clam's headwater channels, left
nearly dry once SWAT+'s aquifers pass recharge through, lost 84 million m3 to MODFLOW 6 that they didn't have.

Other limits of this version:
- MODFLOW 6's river stage is fixed in its input unless the example's ``[rivers]`` section sets it from the surface
  model's channel depth (``stage = "surface_depth"``; SWAT+'s depth needs swatplus-bmi 0.2.0);
- the lag moves water by a day but keeps its volume; the last day's exchange is reported as in transit.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .modflow6 import Modflow6


class ChannelMap:
    """Which entries of a MODFLOW 6 river package belong to which SWAT+ channel."""

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
        with open(folder / spec, newline="") as f:
            for r in csv.DictReader(f):
                ch = r["object"]
                if ch not in swat_channels:
                    raise KeyError(f"{spec}: {ch!r} is not a SWAT+ channel")
                if ch not in chans:
                    chans.append(ch)
                n = int(r["node"]) if r.get("node") else mf6.node(int(r["layer"]), int(r["row"]), int(r["col"]))
                hit = np.flatnonzero(nodes == n)
                if len(hit) == 0:
                    raise ValueError(f"{spec}: cell {n} has no {package} entry")
                entry[hit] = chans.index(ch)
        return cls(chans, entry)


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
    """The river package's exchange, plus any drains coupled to the same channels (``drains``: package, map).
    Drains only take water from the aquifer (seepage at land surface) and give it to their channel the next day,
    as the rivers' gains do."""

    def __init__(self, mf6: Modflow6, package: str, cmap: ChannelMap, drains: tuple[str, ChannelMap] | None = None):
        self.mf6, self.package, self.map = mf6, package, cmap
        self.drains = drains
        if drains is not None and drains[1].channels != cmap.channels:
            raise ValueError("the drains' channels must be the rivers' channels, in the same order")
        P = package.upper()
        self._cond = mf6.get("COND", mf6.model, P).astype(float)
        self._rbot = mf6.get("RBOT", mf6.model, P).astype(float)

    def limit(self, available: dict[str, float]):
        """The supply limit for today's solve, as a callback for ``Modflow6.solve_day``: after each converged outer
        iteration it lowers, from the day's own stages, the stages of any channel whose net exchange (rivers and
        drains together, at the heads as they are then) takes more than ``available``, and says whether anything
        changed, so the solve goes on until heads and stages agree. ``available`` (m3 for the day) is the channel's
        flow without today's exchange: the exchange is taken tomorrow, when today's won't be there. Counting it let
        a channel lose on alternate days what it gained the day before (the Clam). ``end_day`` puts the day's stages
        back after the solve."""
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
