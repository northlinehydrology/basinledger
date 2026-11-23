"""The stream side of groundwater transport (stage 2): one species' mass across every link between SWAT+ and MODFLOW
6's transport model, as the water crosses them (docs/design.md, Water quality across links). Nitrate (``no3``) and
soluble phosphorus (``solp``), the species SWAT+ carries by default, through swatplus-bmi.

**Recharge, HRU by HRU** (``[recharge] source = "hru"``): each HRU's mass leached below its soil profile on the day
just done, SWAT+'s own records read through swatplus-bmi's generated reads: nitrate ``hru_pw:percn`` (``percn``, the
nitrate in the bottom layer's percolate, nut_nlch.f90, which SWAT+ also puts in the HRU's recharge hydrograph,
hru_hyds.f90) and soluble P ``hru_ls:lchlabp`` (the labile P leached past the bottom layer, nut_solp.f90, which SWAT+
62.0.1 leaves out of the recharge hydrograph, so in SWAT+ alone it leaves the model). Both kg ha-1. The mass is
spread over the cells by the same map, scaled to each HRU's area, as the water (``AreaMap.conserving``), so each
cell's recharge takes the concentration of the mass and water its HRUs sent it. A cell given mass and no water is
booked (``no_water``): a recharge package can't carry mass without water. With ``[recharge] source = "aquifer"`` the
mass is each SWAT+ aquifer's recharge mass, spread over all its cells, a warning says so; soluble P has no aquifer
path in SWAT+ and is refused.

**The aquifer to the streams.** What the transport model's coupled rivers and drains discharge to a channel today
(gains and drains, at their cells' concentrations) is due in SWAT+ tomorrow, with the water. Where the day's net
exchange for the channel is a gain, the gain's mass (``channel_water_<sp>__groundwater_mass_exchange_rate``) is what
the gains and drains carried, less what the channel's losing entries took the same day at SWAT+'s outflow
concentration for the day just run (none below zero): SWAT+ receives the net, as it receives the net water.

**The streams to the aquifer.** A loss's water enters the aquifer today at no concentration, and its mass a day later,
when SWAT+ has given it: tomorrow, once SWAT+ has applied the exchange, what it reports it received (negative where
its loss took mass, at its own outflow concentration) is set against what the aquifer sent it, and the difference
goes back into the aquifer through the mass sources on the channel's cells (``src_<package>``, add_transport's
stream_mass), on the entries that lost water the day before, by the volume each lost (on the entries that gained,
by their mass, if none lost). So the link closes exactly, every day, whatever the channel's concentration did
overnight: sent = received in SWAT+ + put back. A loss's mass reaches its cells a day after its water [A], the
exchange's own lag.

**Into the soil** (``[groundwater_to_soil]``): the soil link's water leaves each cell at the cell's concentration and
is split among the HRUs' layers as the water is (``GroundwaterToSoil.mass_split``); SWAT+ takes it the next day with
the water (``soil_layer_water_<sp>__groundwater_inflow_mass_rate``).

Masses kg, concentrations kg m-3."""

from __future__ import annotations

import warnings

import numpy as np

from ..units import DAY

SPECIES = ("no3", "solp")
LEACHED = {"no3": "hru_pw:percn", "solp": "hru_ls:lchlabp"}   # kg ha-1, the day just done, per HRU
AQUIFER_RECHARGE = "aquifer_water_{x}__recharge_mass_rate"      # kg s-1, per SWAT+ aquifer
GAIN = "channel_water_{x}__groundwater_mass_exchange_rate"      # kg s-1, set with a gain, read back as in or out
CONCENTRATION = "channel_water_{x}__mass_concentration"         # kg m-3, the day's outflow
SOIL = "soil_layer_water_{x}__groundwater_inflow_mass_rate"     # kg s-1 per soil layer, set with the water


class AquiferRechargeMass(UserWarning):
    """Recharge mass mapped per SWAT+ aquifer rather than HRU by HRU."""


def spread(amount: float, *weights: np.ndarray) -> np.ndarray:
    """``amount`` over entries by the first set of weights with a positive sum (each an array over the same entries);
    by equal shares of the last's entries if none has."""
    for w in weights:
        w = np.clip(np.asarray(w, float), 0.0, None)
        if w.sum() > 0:
            return amount * w / w.sum()
    w = np.ones(len(weights[-1]))
    return amount * w / w.sum() if len(w) else w


class SurfaceMass:
    """One species' mass across the links, for the driver. ``x``: SWAT+'s species (``no3`` or ``solp``).
    ``recharge``: recharge mass from SWAT+; ``stream``: the stream side; ``soil``: the soil link's mass."""

    def __init__(self, sw, transport, x: str, recharge: bool, stream: bool, soil: bool, source: str = "hru",
                 source_area: np.ndarray | None = None):
        if x not in SPECIES:
            raise ValueError(f"surface species {x!r}: SWAT+ carries {SPECIES} through swatplus-bmi")
        self.sw, self.t, self.x = sw, transport, x
        self.recharge, self.stream, self.soil = recharge, stream, soil
        self.source = source
        self.source_area = None if source_area is None else np.asarray(source_area, float)
        need = []
        if recharge:
            if source == "hru":
                need.append(LEACHED[x])
            else:
                if x == "solp":
                    raise ValueError("soluble P has no recharge path through SWAT+'s aquifers (SWAT+ 62.0.1 leaves "
                                     "leached P out of the recharge hydrograph): its recharge mass is read HRU by "
                                     "HRU, so it needs [recharge] source = \"hru\"")
                need.append(AQUIFER_RECHARGE.format(x=x))
        if stream:
            need += [GAIN.format(x=x), CONCENTRATION.format(x=x)]
        if soil:
            need.append(SOIL.format(x=x))
        missing = [v for v in need if not sw.has(v)]
        if missing:
            raise ValueError(f"[transport] {x}: swatplus-bmi has no {missing} (0.8.0 or later)")
        if recharge and source == "aquifer":
            warnings.warn(f"[transport] {x}: recharge mass is SWAT+'s per aquifer ({AQUIFER_RECHARGE.format(x=x)}), "
                          "spread evenly over each aquifer's cells, NOT HRU by HRU: every cell under an aquifer "
                          "takes the same concentration. [recharge] source = \"hru\" maps it as the water is",
                          AquiferRechargeMass, stacklevel=3)
        if recharge and source == "hru" and sw.names(LEACHED[x]) != sw.names("soil_water__percolation_rate"):
            raise ValueError(f"{LEACHED[x]}'s HRUs aren't soil_water__percolation_rate's")
        # the stream side's state: per channel, the mass due in SWAT+ tomorrow and what the aquifer sent today
        self.due: dict[str, float] = {}       # kg, the gain's mass set in SWAT+ tomorrow
        self.sent: dict[str, float] = {}      # kg, to the channel today (gains and drains)
        self.weights: dict[str, list] = {}    # per channel, per package: (entries, lost m3, discharged kg)
        self.in_swat: dict[str, float] = {}   # kg, what SWAT+ reported it received from the exchange today
        self.c_out: dict[str, float] = {}     # kg m-3, SWAT+'s outflow concentration on the day just run
        self.soil_due = None                  # kg per soil layer, SWAT+'s next day
        self.soil_split = 0.0

    # ---------- recharge ----------
    def recharge_mass(self, link, mf6) -> tuple[np.ndarray, float, float, float]:
        """Per recharge entry, the concentration of today's recharge (kg m-3); and the day's mass made in SWAT+, sent
        with water, and sent to cells with no water (kg). After the recharge link's transfer."""
        if self.source == "hru":
            kg_m2 = self.sw.get(LEACHED[self.x]) / 1e4  # kg ha-1 -> kg m-2, the day
            if np.isnan(kg_m2).any():
                raise ValueError(f"{LEACHED[self.x]} is NaN: SWAT+ made no record for the day (print.prt's nyskip "
                                 "years; the driver sets nyskip to 0 in the run's copy)")
        else:
            kg_m2 = self.sw.get(AQUIFER_RECHARGE.format(x=self.x)) * DAY / self.source_area
        made = float(np.sum(kg_m2 * self.source_area))
        area = mf6.cell_area
        cell_kg = link.map.to_cells(kg_m2, area) * area
        cell_m3 = link.map.to_cells(link.rate, area) * area
        wet = cell_m3 > 0
        c = np.divide(cell_kg, cell_m3, out=np.zeros_like(cell_kg), where=wet)
        nodes = self.t.entry_nodes(link.package)
        return (np.where(nodes >= 0, c[np.maximum(nodes, 0)], 0.0), made, float(cell_kg[wet].sum()),
                float(cell_kg[~wet].sum()))

    # ---------- SWAT+'s day ----------
    def before_swat(self, water_due: dict[str, float]):
        """Set the gains' mass for SWAT+'s coming day, with the gains' water already set (the channels whose due
        water is a gain), and the soil link's mass with its water. Returns the mass set (kg) per channel, and in
        the soil."""
        set_ = {}
        if self.stream:
            var = GAIN.format(x=self.x)
            names = self.sw.names(var)
            idx, vals = [], []
            for i, c in enumerate(names):
                if water_due.get(c, 0.0) > 0 and c in self.due:
                    idx.append(i)
                    vals.append(self.due[c] / DAY)
                    set_[c] = self.due[c]
            if idx:
                self.sw.bmi.set_value_at_indices(var, np.array(idx), np.asarray(vals, np.float64))
        soil = 0.0
        if self.soil and self.soil_due is not None:
            self.sw.set(SOIL.format(x=self.x), self.soil_due / DAY)
            soil = float(self.soil_due.sum())
        return set_, soil

    def after_swat(self, channels: list[str]) -> float:
        """Read what SWAT+'s day did: each exchanged channel's mass received (kg, negative where a loss took mass)
        and its outflow concentration; returns the soil's mass received (kg)."""
        if self.stream:
            got = dict(zip(self.sw.names(GAIN.format(x=self.x)), self.sw.get(GAIN.format(x=self.x)) * DAY))
            self.in_swat = {c: float(got[c]) for c in channels}
            self.c_out = dict(zip(self.sw.names(CONCENTRATION.format(x=self.x)),
                                  np.nan_to_num(self.sw.get(CONCENTRATION.format(x=self.x)))))
        if self.soil:
            return float(self.sw.get(SOIL.format(x=self.x)).sum() * DAY)
        return 0.0

    # ---------- MODFLOW 6's day ----------
    def put_back(self, streams) -> dict[str, float]:
        """Set today's mass sources on the channels' cells: what the aquifer sent each channel yesterday less what
        SWAT+ reports it received, spread as ``spread`` says over yesterday's losing entries (or gaining ones). After
        ``prepare_day``, before the solve. Returns kg per channel."""
        back = {}
        rates = {pkg.upper(): np.zeros(self.t.src_nbound(pkg)) for pkg, _, _ in streams if self.t.has_stream_mass(pkg)}
        for c in set(self.sent) | set(self.in_swat):
            a = self.sent.get(c, 0.0) - self.in_swat.get(c, 0.0)
            back[c] = a
            if a == 0.0:
                continue
            parts = self.weights.get(c)
            if not parts:
                raise ValueError(f"channel {c}: mass to put back ({a} kg) and no entries to put it on")
            lost = np.concatenate([p[1] for p in parts])
            out = np.concatenate([p[2] for p in parts])
            share = spread(a, lost, out, np.ones(len(lost)))
            i = 0
            for pkg, lost_p, _, entries in parts:
                if pkg not in rates:
                    raise ValueError(f"{self.t.model.lower()}: no mass source src_{pkg.lower()} on {pkg}'s cells, "
                                     f"for the mass channel {c}'s losses carry (add_transport's stream_mass)")
                rates[pkg][entries] += share[i:i + len(lost_p)]
                i += len(lost_p)
        for pkg, r in rates.items():
            self.t.set_stream_mass(pkg, r)
        return back

    def after_solve(self, streams, water_due: dict[str, float], budget) -> dict[str, float]:
        """From the day just solved: what the aquifer sent each channel (kg), each entry's weight for tomorrow's
        put-back, and the gain's mass due in SWAT+ tomorrow where the channel's net water is a gain: the gains' and
        drains' mass less the losing entries' water at SWAT+'s outflow concentration for the day just run, not below
        zero. Returns the due mass per channel."""
        self.sent = dict(budget.to_stream)
        self.weights = {}
        lost_m3: dict[str, float] = {}
        for pkg, entry_channel, channels in streams:
            q, m = self.t.last[pkg.upper()]
            e = np.asarray(entry_channel[: len(q)])
            order = np.argsort(e, kind="stable")  # each channel's entries together, found once (many channels)
            bounds = np.searchsorted(e[order], np.arange(len(channels) + 1))
            for i, c in enumerate(channels):
                idx = order[bounds[i]:bounds[i + 1]]
                if not len(idx):
                    continue
                lost = np.clip(q[idx], 0.0, None)
                self.weights.setdefault(c, []).append((pkg.upper(), lost, np.clip(-m[idx], 0.0, None), idx))
                lost_m3[c] = lost_m3.get(c, 0.0) + float(lost.sum())
        self.due = {}
        for c, w in water_due.items():
            if w > 0:
                self.due[c] = max(self.sent.get(c, 0.0) - self.c_out.get(c, 0.0) * lost_m3.get(c, 0.0), 0.0)
        return self.due

    def soil_after_solve(self, gw_soil, x: np.ndarray) -> float:
        """The soil link's mass for SWAT+'s next day, split as its water was; returns its total (kg)."""
        self.soil_due = gw_soil.mass_split(x)
        self.soil_split = float(self.soil_due.sum())
        return self.soil_split

    # ---------- restarts ----------
    def state(self) -> dict:
        return {"due": self.due, "sent": self.sent, "in_swat": self.in_swat, "c_out": self.c_out,
                "soil_due": None if self.soil_due is None else self.soil_due.tolist(),
                "weights": {c: [(p, lost.tolist(), out.tolist(), idx.tolist()) for p, lost, out, idx in parts]
                            for c, parts in self.weights.items()}}

    def load(self, s: dict):
        self.due, self.sent, self.in_swat, self.c_out = (dict(s[k]) for k in ("due", "sent", "in_swat", "c_out"))
        self.soil_due = None if s["soil_due"] is None else np.asarray(s["soil_due"], float)
        self.weights = {c: [(p, np.asarray(lost), np.asarray(out), np.asarray(idx, int)) for p, lost, out, idx in parts]
                        for c, parts in s["weights"].items()}
