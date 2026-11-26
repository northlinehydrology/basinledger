"""Feedback rules: decisions made at the pause between days, from the models' state, applied to the coming day.

A rule reads the state at the end of yesterday (SWAT+ channel flow, MODFLOW 6 heads) and returns controls for a user's
coming day: ``limit`` (a cap on withdrawal, m3 s-1) and ``dry_mode`` (1 = run without evaporative cooling). The user
component must take that control as an input (``users.CONTROLS``); otherwise the rule is refused when the run starts.

Built-in rules (``[[rules]]`` in an example's config):

    kind = "curtailment"      user, channel, threshold (m3 s-1), limit (m3 s-1, default 0)
                              cap the user's withdrawal the day after its channel's flow was below the threshold
    kind = "dry_mode"         user, channel, threshold (m3 s-1)
                              run the user dry the day after the flow was below the threshold
    kind = "drawdown_limit"   user, cell [layer, row, col] or [layer, icell2d] (DISV), head_min (m),
                              limit (m3 s-1, default 0)
                              cap the user's withdrawal the day after the head in that cell was below head_min

Operating rules (docs/design.md, "Operating rules for systems and reservoirs"); a transfer's rate is 0 on a day no
rule sets it, and a reservoir's pool is its normal volume (``pool`` m3, or SWAT+'s principal spillway volume):

    kind = "pumped_storage"   transfer, channel, above (m3 s-1), reservoir, below (share of the pool, default 1),
                              pool (m3, optional)
                              pump the river's flow above ``above`` into the reservoir while it is below
                              ``below`` x pool: the least of the transfer's capacity, the flow above the threshold and
                              the room to the target over a day
    kind = "release"          transfer, channel, below (m3 s-1), reservoir and keep (share of the pool, optional),
                              pool (optional)
                              release the channel's deficit below ``below``: the least of the capacity, the deficit and
                              what the reservoir holds above ``keep`` x pool over a day
    kind = "restriction"      utility, stages = [{when = {channel, below} or {reservoir, below (share)}, cut =
                              {class = share}}, ...]
                              cut the system's customer classes by the last stage whose trigger held yesterday

A study can add its own: ``kind = "rules:my_rule"`` names a function in the example folder,
``my_rule(state, **settings) -> {user: {control: value}}``.

When several rules give the same user the same control, the strictest wins: the smallest limit, and dry if any says
dry; the smallest rate for a transfer (``transfer:<name>``), the largest cut for a system's class
(``utility:<name>``, control ``cut:<class>``). Every rule's decision is entered in the ledger (``rule:<name>:active``, 1 or 0) so a run shows when it acted.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field


@dataclass
class State:
    """What a rule may read at the pause before day ``day``: yesterday's results."""

    day: int
    flows: dict[str, float]                 # channel -> outflow, m3 s-1, yesterday
    head_at: object                         # cell (layer, row, col or layer, icell2d) -> head, m, after yesterday
    memory: dict = field(default_factory=dict)  # a rule's own memory between days (saved with restarts)
    volumes: dict[str, float] = field(default_factory=dict)  # reservoir -> volume, m3, at the end of yesterday
    pools: dict[str, float] = field(default_factory=dict)    # reservoir -> normal pool, m3 (principal spillway)
    capacities: dict[str, float] = field(default_factory=dict)  # transfer -> its capacity, m3 s-1

    def flow(self, channel: str, rule: str) -> float:
        if channel not in self.flows:
            raise KeyError(f"{rule}: channel {channel!r} not in SWAT+")
        return self.flows[channel]

    def volume(self, reservoir: str, rule: str) -> float:
        if reservoir not in self.volumes:
            raise KeyError(f"{rule}: reservoir {reservoir!r} not in SWAT+")
        return self.volumes[reservoir]

    def pool(self, reservoir: str, given: float | None) -> float:
        return float(given) if given is not None else self.pools[reservoir]


def curtailment(state: State, *, user: str, channel: str, threshold: float, limit: float = 0.0):
    q = state.flows.get(channel)
    if q is None:
        raise KeyError(f"curtailment: channel {channel!r} not in SWAT+")
    return {user: {"limit": float(limit)}} if (state.day > 0 and q < threshold) else {}


def dry_mode(state: State, *, user: str, channel: str, threshold: float):
    q = state.flows.get(channel)
    if q is None:
        raise KeyError(f"dry_mode: channel {channel!r} not in SWAT+")
    return {user: {"dry_mode": 1.0}} if (state.day > 0 and q < threshold) else {}


def drawdown_limit(state: State, *, user: str, cell, head_min: float, limit: float = 0.0):
    if state.day == 0:
        return {}
    h = state.head_at(tuple(cell))
    return {user: {"limit": float(limit)}} if h < head_min else {}


DAY = 86400.0


def pumped_storage(state: State, *, transfer: str, channel: str, above: float, reservoir: str,
                   below: float = 1.0, pool: float | None = None):
    if state.day == 0:
        return {}
    q = state.flow(channel, "pumped_storage")
    room = below * state.pool(reservoir, pool) - state.volume(reservoir, "pumped_storage")
    rate = min(state.capacities[transfer], q - above, room / DAY)
    return {f"transfer:{transfer}": {"rate": rate if rate > 0 else 0.0}}


def release(state: State, *, transfer: str, channel: str, below: float, reservoir: str | None = None,
            keep: float = 0.0, pool: float | None = None):
    if state.day == 0:
        return {}
    rate = min(state.capacities[transfer], below - state.flow(channel, "release"))
    if reservoir is not None:
        rate = min(rate, (state.volume(reservoir, "release") - keep * state.pool(reservoir, pool)) / DAY)
    return {f"transfer:{transfer}": {"rate": rate if rate > 0 else 0.0}}


def _triggered(state: State, when: dict, where: str) -> bool:
    if set(when) == {"channel", "below"}:
        return state.flow(when["channel"], where) < when["below"]
    if set(when) - {"reservoir", "below", "pool"} or not {"reservoir", "below"} <= set(when):
        raise ValueError(f"{where}: a stage's when is {{channel, below}} (m3 s-1) or {{reservoir, below}} (a share of "
                         f"its pool, which pool = m3 may give); got {when!r}")
    return state.volume(when["reservoir"], where) < when["below"] * state.pool(when["reservoir"], when.get("pool"))


def restriction(state: State, *, utility: str, stages: list):
    if state.day == 0:
        return {}
    cut = None
    for i, s in enumerate(stages, 1):
        if _triggered(state, s["when"], f"restriction on {utility}, stage {i}"):
            cut = s["cut"]
    return {f"utility:{utility}": {f"cut:{c}": float(v) for c, v in cut.items()}} if cut else {}


BUILT_IN = {"curtailment": curtailment, "dry_mode": dry_mode, "drawdown_limit": drawdown_limit,
            "pumped_storage": pumped_storage, "release": release, "restriction": restriction}
TARGET = {"pumped_storage": "transfer", "release": "transfer", "restriction": "utility"}  # what else a rule acts on


def check(kind: str, settings: dict, where: str):
    """An operating rule's settings, at the start."""
    if kind == "restriction":
        stages = settings.get("stages")
        if not isinstance(stages, list) or not stages:
            raise ValueError(f"{where}: stages is a list of {{when, cut}}")
        for s in stages:
            if not isinstance(s, dict) or set(s) != {"when", "cut"} or not isinstance(s["cut"], dict) or \
                    not all(isinstance(v, (int, float)) and 0 <= v <= 1 for v in s["cut"].values()):
                raise ValueError(f"{where}: each stage is {{when = {{...}}, cut = {{class = share 0 to 1}}}}; got {s!r}")
            w = s["when"]
            if not isinstance(w, dict) or not (set(w) == {"channel", "below"} or
                                               {"reservoir", "below"} <= set(w) <= {"reservoir", "below", "pool"}):
                raise ValueError(f"{where}: a stage's when is {{channel, below}} (m3 s-1) or {{reservoir, below}} (a "
                                 f"share of its pool, which pool = m3 may give); got {w!r}")


@dataclass
class Rule:
    name: str
    fn: object
    settings: dict

    def __call__(self, state: State) -> dict:
        return self.fn(state, **self.settings) or {}


def combine(decisions: list[dict]) -> dict[str, dict[str, float]]:
    """The strictest control per user: the smallest limit, and dry if any rule says dry."""
    out: dict[str, dict[str, float]] = {}
    for d in decisions:
        for user, controls in d.items():
            c = out.setdefault(user, {})
            for k, v in controls.items():
                if k == "limit":
                    c[k] = min(c.get(k, math.inf), v)
                elif k == "dry_mode" or k.startswith("cut:"):
                    c[k] = max(c.get(k, 0.0), v)
                elif k == "rate":
                    c[k] = min(c.get(k, math.inf), v)
                else:
                    raise KeyError(f"unknown control {k!r}")
    return out
