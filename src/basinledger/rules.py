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

A study can add its own: ``kind = "rules:my_rule"`` names a function in the example folder,
``my_rule(state, **settings) -> {user: {control: value}}``.

When several rules give the same user the same control, the strictest wins: the smallest limit, and dry if any says
dry. Every rule's decision is entered in the ledger (``rule:<name>:active``, 1 or 0) so a run shows when it acted.
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


BUILT_IN = {"curtailment": curtailment, "dry_mode": dry_mode, "drawdown_limit": drawdown_limit}


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
                elif k == "dry_mode":
                    c[k] = max(c.get(k, 0.0), v)
                else:
                    raise KeyError(f"unknown control {k!r}")
    return out
