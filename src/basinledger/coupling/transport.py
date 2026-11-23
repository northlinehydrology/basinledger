"""Groundwater transport: MODFLOW 6 GWT models on the flow model's grid, one per species, stepped with the flow.

**Building.** ``add_transport`` is called by an example's builder, after the flow model is made and before the
simulation is written. For each species it adds a GWT model named ``gwt_<species>`` with its own IMS (MODFLOW 6
solves flow first, then each transport model, every time step), the GWF6-GWT6 exchange, and the packages: mobile
storage with porosity and optional sorption (MST), advection (ADV), dispersion (DSP), starting concentrations (IC),
sources and sinks (SSM) and output control (OC). Each flow package that brings water in with a concentration is named
in ``sources`` with its auxiliary variable, which is added to the flow package (zero) if it lacks it. Water leaving
the aquifer through any package leaves at its cell's concentration. Water entering through a flow package not named
in ``sources`` enters at zero concentration (MODFLOW 6's rule).

**Each day** (``Transport``, from the driver), after ``Modflow6.prepare_day``:

    set_inflow("rch", c)                 the recharge's concentration, kg m-3 ([A] in stage 1)
    set_inflow("riv", c, entries=...)    that of the water the coupled river entries lose to the aquifer ([A])
    (the flow and transport solves, Modflow6.solve_day)
    concentration(), budget(...)         the day's concentrations and its mass budget, kg

A concentration is set on the transport model's own copy of the flow package's auxiliary values
(``GWT_<NAME>/FMI-FT<n>/AUXVAR``): MODFLOW 6 copies the flow package's values there when it prepares the step, so a
value set on the flow package itself (``GWF/RCH/AUXVAR``) after that has no effect, and one set on the flow package's
input before it takes effect a step late.

**Mass budget.** Every flow package's entries carry mass: inflow at the concentration SSM gives it (the auxiliary
value, or zero), outflow at the cell's concentration at the end of the day (MODFLOW 6's implicit step). The budget
splits it as the water ledger does (recharge, wells, coupled streams, the soil link, boundaries) and adds the release
from mobile and sorbed storage (MST's rates; the sorbed part on its own), the loss to first-order decay (MST's decay
rates, dissolved and sorbed) and GWT's own boundary packages (a constant concentration, a mass source). A mass source
named ``src_<package>`` on a coupled stream package's cells is the stream side's (coupling/surface_mass.py), booked to
the entries' channels. The sum closes to the transport solve's convergence; ``closure_share`` is the residual over the
day's gross mass flow.

**Steady state** (``track`` each day, then ``steady_state``): the concentrations each cell settles at under the mean
flows and inflows of the days tracked, with decay, by advection and mixing (dispersion left out [A]).

Concentrations in kg m-3 (= g L-1), masses in kg. A model in other units (``Modflow6``, Units) holds its
concentrations as kg per its own volume unit, and they are converted with its flows: a transport model in feet holds
kg ft-3. Its length units must be the flow model's.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

def _dispersivities(dispersivity) -> dict:
    """A longitudinal dispersivity (m), with transverse horizontal at a tenth and transverse vertical at a hundredth
    [A] (the usual ratios, e.g. Gelhar and others 1992), or the three given (longitudinal, transverse horizontal,
    transverse vertical)."""
    if np.ndim(dispersivity) == 0:
        al = float(dispersivity)
        return dict(alh=al, ath1=al / 10.0, atv=al / 100.0)
    al, at, av = (float(x) for x in dispersivity)
    return dict(alh=al, ath1=at, atv=av)


def _add_auxiliary(package, aux: str):
    """Give a flow package an auxiliary variable (zero), if it hasn't one by that name: an array package (READASARRAYS,
    as RCHA) gets a zero array in each period, a list package a zero column."""
    have = [a for row in (package.auxiliary.get_data() or []) for a in row if a.lower() != "auxiliary"]
    if aux.lower() in (a.lower() for a in have):
        return
    package.auxiliary = have + [aux]
    arrays = getattr(package, "readasarrays", None)
    if arrays is not None and arrays.get_data():
        old = package.aux.get_data() or {0: []}
        package.aux = {k: [*(v if v is not None else []), 0.0] for k, v in old.items()}
        return
    package.stress_period_data = {k: [(*r, 0.0) for r in rec.tolist()]
                                  for k, rec in package.stress_period_data.get_data().items() if rec is not None}


def _by_zone(value, zones, what: str):
    """A per-cell value given as a scalar, an array, or ``{zone: value}`` over ``zones`` (an integer array per cell,
    the grid's shape): every zone in the array must have a value."""
    if not isinstance(value, dict):
        return value
    if zones is None:
        raise ValueError(f"{what} given by zone needs zones (an integer array per cell)")
    z = np.asarray(zones)
    vals = {int(k): float(v) for k, v in value.items()}
    missing = sorted(set(np.unique(z).tolist()) - set(vals))
    if missing:
        raise ValueError(f"{what}: zones {missing} have no value")
    out = np.zeros(z.shape)
    for k, v in vals.items():
        out[z == k] = v
    return out


def add_transport(sim, gwf, species: str, porosity, dispersivity, sorption: dict | None = None, strt=0.0,
                  sources: list[tuple[str, str]] = (), scheme: str = "tvd", diffusion: float = 0.0,
                  ims: dict | None = None, save: str = "LAST", decay: float | None = None,
                  decay_sorbed: float | None = None, stream_mass: list[str] = ()):
    """A GWT model ``gwt_<species>`` on ``gwf``'s grid (DIS or DISV), linked to it and solved after it.

    porosity: effective porosity [A] (scalar or per cell). dispersivity: m [A], see ``_dispersivities``.
    sorption: None (a conservative species, such as chloride or nitrate), or ``{"distcoef": Kd (m3 kg-1),
    "bulk_density": kg m-3}`` for a linear isotherm [A] (soluble phosphorus); Kd may be a scalar, an array per cell,
    or ``{zone: Kd}`` with ``"zones"``, an integer array per cell. strt: starting concentration, kg m-3. sources:
    ``(flow package pname, auxiliary name)`` for each package whose inflow carries a concentration. scheme: ADV's
    (``upstream``, ``central`` or ``tvd``). diffusion: molecular diffusion coefficient, m2 d-1 [A]. ims: settings for
    the transport solution, over the defaults below. save: when OC saves concentrations and the budget (``LAST`` or
    ``ALL``). decay: first-order decay of the dissolved species, d-1 [A] (denitrification, for nitrate); with sorption,
    ``decay_sorbed`` (d-1 [A]) is the sorbed mass's, the dissolved rate unless given. stream_mass: flow packages
    (the coupled rivers and drains) whose cells get a mass source of their own, ``src_<package>``, one entry per the
    package's entry (its first stress period's cells), through which Basin Ledger puts back the mass a stream's
    losses carry (coupling/transport.py, the stream side)."""
    import flopy

    name = f"gwt_{species}"
    dis = next((p for p in (gwf.get_package("dis"), gwf.get_package("disv")) if p is not None), None)
    if dis is None or dis.package_type not in ("dis", "disv"):
        raise NotImplementedError("add_transport supports a structured (DIS) or vertex (DISV) grid")
    gwt = flopy.mf6.ModflowGwt(sim, modelname=name, save_flows=True)
    # the residual's criterion, 1e-8 kg d-1 (1e-10 before Nov 22, 2026): at 1e-10 the residual's round-off held it
    # above the criterion for good on the synthetic grid with sorption (Kd 2e-3 m3 kg-1: the sorbed store's terms are
    # large) and on Plaster Creek's 13,215 cells without, on day 47. The day's mass closes to the budget either way
    # (closure_share about 1e-15 on Plaster Creek at 1e-8), which the ledger checks
    settings = dict(complexity="simple", linear_acceleration="bicgstab", outer_maximum=50, outer_dvclose=1e-10,
                    inner_maximum=300, inner_dvclose=1e-12, rcloserecord=[1e-8, "strict"])
    settings.update(ims or {})
    solver = flopy.mf6.ModflowIms(sim, filename=f"{name}.ims", **settings)
    sim.register_ims_package(solver, [name])
    units = dis.length_units.get_data()  # the flow model's own, so the two grids are the same
    if dis.package_type == "dis":
        flopy.mf6.ModflowGwtdis(gwt, length_units=units, nlay=dis.nlay.get_data(), nrow=dis.nrow.get_data(),
                                ncol=dis.ncol.get_data(), delr=dis.delr.get_data(), delc=dis.delc.get_data(),
                                top=dis.top.get_data(), botm=dis.botm.get_data(), idomain=dis.idomain.get_data())
    else:  # the flow model's vertex grid, cell for cell
        flopy.mf6.ModflowGwtdisv(gwt, length_units=units, nlay=dis.nlay.get_data(), ncpl=dis.ncpl.get_data(),
                                 nvert=dis.nvert.get_data(), vertices=dis.vertices.get_data(),
                                 cell2d=dis.cell2d.get_data(), top=dis.top.get_data(), botm=dis.botm.get_data(),
                                 idomain=dis.idomain.get_data())
    mst = dict(porosity=porosity)
    if sorption:
        if set(sorption) - {"distcoef", "bulk_density", "zones"}:
            raise ValueError(f"sorption: unknown settings "
                             f"{sorted(set(sorption) - {'distcoef', 'bulk_density', 'zones'})}")
        mst.update(sorption="linear", distcoef=_by_zone(sorption["distcoef"], sorption.get("zones"), "distcoef"),
                   bulk_density=sorption["bulk_density"])
    if decay is not None:
        if decay < 0 or (decay_sorbed is not None and decay_sorbed < 0):
            raise ValueError("decay is a loss: a rate of 0 or more, d-1")
        mst.update(first_order_decay=True, decay=float(decay))
        if sorption:
            mst["decay_sorbed"] = float(decay if decay_sorbed is None else decay_sorbed)
    elif decay_sorbed is not None:
        raise ValueError("decay_sorbed needs decay")
    flopy.mf6.ModflowGwtmst(gwt, **mst)
    flopy.mf6.ModflowGwtadv(gwt, scheme=scheme)
    flopy.mf6.ModflowGwtdsp(gwt, diffc=diffusion, **_dispersivities(dispersivity))
    flopy.mf6.ModflowGwtic(gwt, strt=strt)
    src = []
    for pname, aux in sources:
        pkg = gwf.get_package(pname)
        if pkg is None:
            raise KeyError(f"add_transport: the flow model has no package {pname!r}")
        _add_auxiliary(pkg, aux)
        src.append((pname, "AUX", aux))
    flopy.mf6.ModflowGwtssm(gwt, sources=src or None)
    for pname in stream_mass:
        pkg = gwf.get_package(pname)
        if pkg is None:
            raise KeyError(f"add_transport: the flow model has no package {pname!r}")
        rec = pkg.stress_period_data.get_data(0)
        cells = [tuple(r[0]) if not np.isscalar(r[0]) else (r[0],) for r in rec]
        flopy.mf6.ModflowGwtsrc(gwt, pname=f"src_{pname}", filename=f"{name}.src_{pname}", maxbound=len(cells),
                                stress_period_data={0: [(c, 0.0) for c in cells]})
    flopy.mf6.ModflowGwtoc(gwt, budget_filerecord=f"{name}.cbc", concentration_filerecord=f"{name}.ucn",
                           saverecord=[("CONCENTRATION", save), ("BUDGET", save)], printrecord=[("BUDGET", save)])
    flopy.mf6.ModflowGwfgwt(sim, exgtype="GWF6-GWT6", exgmnamea=gwf.name, exgmnameb=name,
                            filename=f"{gwf.name}_{name}.gwfgwt")
    return gwt


@dataclass
class MassBudget:
    """One day of a transport model's mass, kg, positive as named. ``closure`` is the residual (into the aquifer
    less out of it less what storage gained less what decayed); ``closure_share`` divides it by the day's gross mass
    flow."""

    recharge: float = 0.0                                      # in with recharge
    wells: dict[str, float] = field(default_factory=dict)      # out with each coupled well
    to_stream: dict[str, float] = field(default_factory=dict)  # out to each coupled channel (rivers' gains, drains)
    from_stream: dict[str, float] = field(default_factory=dict)  # in from each coupled channel (rivers' losses)
    leakage_m3: dict[str, float] = field(default_factory=dict)  # the water those losses carried, m3
    returned: dict[str, float] = field(default_factory=dict)   # in through each channel's mass sources (stream side)
    gw_soil: float = 0.0                                       # out with the water up into the soil
    boundary: dict[str, float] = field(default_factory=dict)   # net, into the aquifer, by package (not coupled)
    storage_release: float = 0.0                               # released from mobile and sorbed storage
    sorbed_release: float = 0.0                                # of it, released from sorbed storage
    decay: float = 0.0                                         # lost to first-order decay (dissolved and sorbed)
    closure: float = 0.0
    gross: float = 0.0

    @property
    def closure_share(self) -> float:
        # against the day's gross mass flow, with a floor of 1 mg so a day with almost no mass isn't noise over noise
        return self.closure / max(self.gross, 1e-6)


class Transport:
    """One species' GWT model in a running simulation (``Modflow6``), as the driver uses it each day.

    ``recharge``: the recharge package the driver sets. ``wells``: name -> (package, entry). ``streams``: (package,
    entry_channel, channels) for each package coupled to the surface model's channels (the exchange's rivers and
    drains), with entry_channel -1 for entries not coupled. ``gw_soil``: the package of the water up into the soil.
    Everything else is a boundary. A mass source named ``src_<package>`` on a coupled stream package's cells
    (``add_transport``'s stream_mass) is the stream side's: its mass is booked to the entry's channel."""

    def __init__(self, mf6, species: str, model: str | None = None):
        self.mf6, self.species = mf6, species
        self.model = (model or f"gwt_{species}").upper()
        api = mf6.api
        try:
            api.get_var_address("X", self.model)
        except Exception:
            raise KeyError(f"the MODFLOW 6 simulation has no transport model {self.model.lower()!r} "
                           "(coupling.transport.add_transport)") from None
        names = set(api.get_input_var_names()) | set(api.get_output_var_names())
        self.copies: dict[str, int] = {}  # flow package -> its FMI copy's number
        n = 1
        while f"{self.model}/FMI-FT{n}/NAME" in names:
            self.copies[str(self._get(f"FMI-FT{n}", "NAME")[0]).strip().upper()] = n
            n += 1
        # SSM: per FMI copy, the auxiliary column a package's inflow takes its concentration from (1-based, 0 for none),
        # and how (0 none, 1 AUX; AUXMIXED and SPC files aren't built)
        self.iaux = self._get("SSM", "IAUXPAK").astype(int) if self.copies else np.zeros(0, int)
        kind = self._get("SSM", "ISRCTYPE").astype(int) if self.copies else np.zeros(0, int)
        if np.any(kind > 1):
            raise NotImplementedError(f"{self.model.lower()}: SSM sources other than AUX (AUXMIXED, SPC files) "
                                      "aren't built")
        idcy = int(self._get("MST", "IDCY")[0])
        if idcy not in (0, 1):
            raise NotImplementedError(f"{self.model.lower()}: zero-order decay isn't built; first-order decay is "
                                      "(add_transport's decay)")
        self.decays = idcy == 1
        self.sorbs = int(self._get("MST", "ISRB")[0]) != 0
        if self.sorbs and int(self._get("MST", "ISRB")[0]) != 1:
            raise NotImplementedError(f"{self.model.lower()}: only linear sorption is built")
        # GWT's own boundary packages (a constant concentration, CNC; a mass source, SRC): every package of the
        # transport model with simulated values, mass from outside the coupled models
        prefix = self.model + "/"
        self.boundaries = sorted({v.split("/")[1] for v in names if v.startswith(prefix) and v.endswith("/SIMVALS")})
        # the stream side's mass sources, src_<flow package>, by the flow package
        self.stream_sources: dict[str, str] = {}
        for p in self.boundaries:
            if p.startswith("SRC_") and p[4:] in self.copies and \
                    str(self._get(p, "PACKAGE_TYPE")[0]).strip().lower() == "src":
                self.stream_sources[p[4:]] = p
        self._src_checked: set[str] = set()
        self.last: dict[str, tuple[np.ndarray, np.ndarray]] = {}  # each package's entries' flow and mass, last day
        self.mean = None  # running sums for the steady state (``track``)
        from .modflow6 import LENUNI
        own = LENUNI.get(int(self._get("DIS", "LENUNI")[0]))
        if own is not None and own != mf6.length_units:
            raise ValueError(f"{self.model.lower()}: in {own}, the flow model in {mf6.length_units}")

    def _get(self, package: str, var: str) -> np.ndarray:
        """In meters, days and kilograms (``Modflow6.get``)."""
        return self.mf6.get(var, self.model, package)

    def _copy(self, package: str) -> int:
        try:
            return self.copies[package.upper()]
        except KeyError:
            raise KeyError(f"{self.model.lower()}: no flow package {package!r} (it has {sorted(self.copies)})") from None

    # ---------- concentrations ----------
    def concentration(self) -> np.ndarray:
        """Concentration per reduced node, kg m-3 (GWT's X): at the end of the last day solved, or the start."""
        return self.mf6.get("X", self.model, "")

    def set_concentration(self, c: np.ndarray):
        """Concentrations per reduced node, for a restart or a steady start: after initialize and before the first
        prepare_day."""
        if self.concentration().shape != np.shape(c):
            raise ValueError("saved concentrations don't match this model's nodes")
        self.mf6.set("X", "", np.asarray(c, dtype=np.float64), model=self.model)

    def takes_inflow(self, package: str) -> bool:
        """Whether SSM gives the package's inflow a concentration (an AUX source)."""
        return int(self.iaux[self._copy(package) - 1]) > 0

    def set_inflow(self, package: str, c, entries: np.ndarray | None = None):
        """The concentration (kg m-3) of water entering through a flow package today: a value for all its entries or
        those selected (a boolean mask over the entries), or an array with a value per entry; after
        ``prepare_day``, before the solve."""
        n = self._copy(package)
        col = int(self.iaux[n - 1]) - 1
        if col < 0:
            raise ValueError(f"{self.model.lower()}: {package!r} isn't an AUX source in SSM, so its inflow carries no "
                             "concentration (add_transport's sources)")
        nb = int(self._get(f"FMI-FT{n}", "NBOUND")[0])
        aux = self._get(f"FMI-FT{n}", "AUXVAR")
        a = aux.reshape(-1, aux.shape[-1]) if aux.ndim > 1 else aux.reshape(-1, 1)
        if np.ndim(c) == 1:
            if entries is not None or len(c) != nb:
                raise ValueError(f"{package}: {len(c)} concentrations for {nb} entries")
            a[:nb, col] = c
        else:
            sel = slice(0, nb) if entries is None else np.flatnonzero(np.asarray(entries[:nb], bool))
            a[sel, col] = c
        self.mf6.set("AUXVAR", f"FMI-FT{n}", a.reshape(aux.shape), model=self.model)

    def has_stream_mass(self, package: str) -> bool:
        return package.upper() in self.stream_sources

    def src_nbound(self, package: str) -> int:
        """The number of entries of a coupled stream package's mass source."""
        return int(self._get(self.stream_sources[package.upper()], "NBOUND")[0])

    def set_stream_mass(self, package: str, rate: np.ndarray):
        """The mass (kg d-1) the stream side puts back on each entry's cell of a coupled stream package today,
        through its mass source ``src_<package>``; after ``prepare_day``, before the solve. The source's cells are
        checked against the package's on the first day."""
        src = self.stream_sources.get(package.upper())
        if src is None:
            raise ValueError(f"{self.model.lower()}: no mass source src_{package.lower()} on {package}'s cells "
                             "(add_transport's stream_mass)")
        nb = int(self._get(src, "NBOUND")[0])
        if src not in self._src_checked:
            n = self._copy(package)
            fnodes = self._get(f"FMI-FT{n}", "NODELIST")[: int(self._get(f"FMI-FT{n}", "NBOUND")[0])].astype(int)
            snodes = self._get(src, "NODELIST")[:nb].astype(int)
            if not np.array_equal(fnodes, snodes):
                raise ValueError(f"{self.model.lower()}: {src.lower()}'s cells aren't {package}'s, entry for entry")
            self._src_checked.add(src)
        r = np.asarray(rate, float)
        if len(r) != nb:
            raise ValueError(f"{src.lower()}: {len(r)} rates for {nb} entries")
        self.mf6.set("SMASSRATE", src, r, model=self.model)

    # ---------- the day's mass ----------
    def _entries(self, package: str, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Each entry's flow (m3 d-1, positive into the aquifer) and mass flow (kg d-1) on the day just solved, with
        the concentrations ``x``."""
        n = self._copy(package)
        F = f"FMI-FT{n}"
        nb = int(self._get(F, "NBOUND")[0])
        q = self._get(F, "FLOW")[:nb].astype(float)
        node = self._get(F, "NODELIST")[:nb].astype(int) - 1  # -1 for an inactive column (RCHA): no flow
        c = np.where(node >= 0, x[np.maximum(node, 0)], 0.0)
        col = int(self.iaux[n - 1]) - 1
        if col >= 0:
            aux = self._get(F, "AUXVAR")
            a = aux.reshape(-1, aux.shape[-1]) if aux.ndim > 1 else aux.reshape(-1, 1)
            c = np.where(q > 0, a[:nb, col], c)
        else:
            c = np.where(q > 0, 0.0, c)  # inflow from a package SSM doesn't list carries none
        return q, q * c

    def entry_nodes(self, package: str) -> np.ndarray:
        """A flow package's entries' reduced nodes, 0-based (-1 for an inactive column)."""
        n = self._copy(package)
        F = f"FMI-FT{n}"
        return self._get(F, "NODELIST")[: int(self._get(F, "NBOUND")[0])].astype(int) - 1

    def budget(self, recharge: str | None = None, wells: dict[str, tuple[str, int]] | None = None,
               streams: list[tuple[str, np.ndarray, list[str]]] = (), gw_soil: str | None = None) -> MassBudget:
        """The day's mass, kg, from the transport model's flows and concentrations after the solve. Each package's
        entries' flow and mass are kept in ``last``."""
        b = MassBudget()
        net = 0.0
        x = self.concentration()
        flows = {pkg: self._entries(pkg, x) for pkg in self.copies}
        self.last = dict(flows)
        coupled: dict[str, np.ndarray] = {}
        for pkg, (q, m) in flows.items():
            net += float(m.sum())
            b.gross += float(np.abs(m).sum())
            coupled[pkg] = np.zeros(len(m), bool)
            if recharge and pkg == recharge.upper():
                b.recharge = float(m.sum())
                coupled[pkg][:] = True
            if gw_soil and pkg == gw_soil.upper():
                b.gw_soil = -float(m.sum())
                coupled[pkg][:] = True
        for name, (pkg, i) in (wells or {}).items():
            q, m = flows[pkg.upper()]
            b.wells[name] = -float(m[i])
            coupled[pkg.upper()][i] = True
        own_coupled: dict[str, np.ndarray] = {}  # the stream side's mass sources: their entries' channels
        for pkg, entry_channel, channels in streams:
            q, m = flows[pkg.upper()]
            e = entry_channel[: len(m)]
            src = self.stream_sources.get(pkg.upper())
            sm = None
            if src is not None:
                sm = self._get(src, "SIMVALS")[: int(self._get(src, "NBOUND")[0])].astype(float)
                self.last[src] = (np.zeros(len(sm)), sm)
                own_coupled[src] = e[: len(sm)] >= 0
            for i, ch in enumerate(channels):
                sel = e == i
                b.to_stream[ch] = b.to_stream.get(ch, 0.0) - float(m[sel & (m < 0)].sum())
                b.from_stream[ch] = b.from_stream.get(ch, 0.0) + float(m[sel & (m > 0)].sum())
                b.leakage_m3[ch] = b.leakage_m3.get(ch, 0.0) + float(q[sel & (q > 0)].sum())
                if sm is not None:
                    b.returned[ch] = b.returned.get(ch, 0.0) + float(sm[sel[: len(sm)]].sum())
            coupled[pkg.upper()] |= e >= 0
        for pkg, (q, m) in flows.items():
            rest = ~coupled[pkg]
            if rest.any():
                b.boundary[pkg.lower()] = float(m[rest].sum())
        for pkg in self.boundaries:  # GWT's own, kg d-1 into the aquifer
            m = self._get(pkg, "SIMVALS")[: int(self._get(pkg, "NBOUND")[0])].astype(float)
            net += float(m.sum())
            b.gross += float(np.abs(m).sum())
            rest = ~own_coupled[pkg] if pkg in own_coupled else np.ones(len(m), bool)
            if rest.any():
                b.boundary[pkg.lower()] = b.boundary.get(pkg.lower(), 0.0) + float(m[rest].sum())
        # MST's rates, positive where storage releases mass into the flow (MODFLOW 6's budget sign), and decay's,
        # negative where it takes mass out of it
        sto = float(self._get("MST", "RATESTO").sum())
        if self.sorbs:
            b.sorbed_release = float(self._get("MST", "RATESRB").sum())
            sto += b.sorbed_release
        b.storage_release = sto
        if self.decays:
            b.decay = -float(self._get("MST", "RATEDCY").sum())
            if self.sorbs:
                b.decay -= float(self._get("MST", "RATEDCYS").sum())
        b.closure = net + sto - b.decay
        b.gross += abs(sto) + abs(b.decay)
        return b

    # ---------- the steady state ----------
    def track(self):
        """Add the day just solved to the running sums ``steady_state`` takes its mean flows from: the flow between
        cells, each cell's saturation, and each cell's inflow of water and mass from the packages (and of mass from
        GWT's own mass sources). Call after ``budget``."""
        if self.mean is None:
            nja = len(self._get("FMI", "GWFFLOWJA"))
            n = len(self.concentration())
            self.mean = {"days": 0, "flowja": np.zeros(nja), "sat": np.zeros(n), "water_in": np.zeros(n),
                         "mass_in": np.zeros(n)}
        s = self.mean
        s["days"] += 1
        s["flowja"] += self._get("FMI", "GWFFLOWJA")
        s["sat"] += self._get("FMI", "GWFSAT")
        n = len(s["sat"])
        for pkg, (q, m) in self.last.items():
            if pkg in self.copies:
                node = self.entry_nodes(pkg)
            else:  # a mass source of GWT's own
                node = self._get(pkg, "NODELIST")[: len(m)].astype(int) - 1
            ok = node >= 0
            inflow = (q > 0) & ok
            s["water_in"] += np.bincount(node[inflow], weights=q[inflow], minlength=n)
            add = ok & (q >= 0)  # mass with water coming in, or with none (a mass source)
            s["mass_in"] += np.bincount(node[add], weights=m[add], minlength=n)
        for pkg in self.boundaries:  # GWT's own mass sources not already in ``last``
            if pkg not in self.last and str(self._get(pkg, "PACKAGE_TYPE")[0]).strip().lower() == "src":
                m = self._get(pkg, "SIMVALS")[: int(self._get(pkg, "NBOUND")[0])].astype(float)
                node = self._get(pkg, "NODELIST")[: len(m)].astype(int) - 1
                s["mass_in"] += np.bincount(node[node >= 0], weights=m[node >= 0], minlength=n)

    def steady_state(self, sweeps: int = 50, tol: float = 1e-12) -> np.ndarray:
        """The concentrations at steady state under the mean of the days ``track`` added: each cell mixes the water
        that comes into it (from upstream cells at their concentrations and from the packages at theirs) and loses
        mass to decay, ``c = mass in / (water in + decay coefficient)``, where the decay coefficient is the first-order
        rate times the cell's dissolved (porosity x saturated volume) and sorbed (bulk density x Kd x volume) capacity.
        Dispersion is left out [A]: the steady state of advection, mixing and decay, as cells upstream to downstream.
        A cell with no water coming in keeps its present concentration. Solved by sweeps from the highest-head cell
        down, the order water flows in, until no cell changes by more than ``tol`` of the largest."""
        if self.mean is None or not self.mean["days"]:
            raise ValueError("steady_state needs the mean flows of at least one day (track)")
        if any(p.startswith("CNC") for p in self.boundaries):
            raise NotImplementedError(f"{self.model.lower()}: a steady state with constant-concentration cells "
                                      "isn't built")
        s, d = self.mean, float(self.mean["days"])
        flowja, sat = s["flowja"] / d, s["sat"] / d
        water_in, mass_in = s["water_in"] / d, s["mass_in"] / d
        ia = self.mf6.api.get_value(self.mf6.api.get_var_address("IA", self.model, "CON")).astype(int) - 1
        ja = self.mf6.api.get_value(self.mf6.api.get_var_address("JA", self.model, "CON")).astype(int) - 1
        n = len(sat)
        vol = self._get("DIS", "AREA") * (self._get("DIS", "TOP") - self._get("DIS", "BOT")) * sat
        coef = np.zeros(n)
        if self.decays:
            coef = self._get("MST", "DECAY") * self._get("MST", "POROSITY") * vol
            if self.sorbs:
                coef = coef + (self._get("MST", "DECAY_SORBED") * self._get("MST", "BULK_DENSITY")
                               * self._get("MST", "DISTCOEF") * vol)
        # each cell's water in from its neighbors (MODFLOW 6's FLOWJA: positive into the row's cell)
        rows = np.repeat(np.arange(n), np.diff(ia))
        off = ja != rows
        inflow = off & (flowja > 0)
        water_in = water_in + np.bincount(rows[inflow], weights=flowja[inflow], minlength=n)
        denom = water_in + coef
        c = self.concentration().copy()
        order = np.argsort(-self._get("FMI", "GWFHEAD"), kind="stable")  # upstream to downstream, as a rule
        starts, ends = ia[:-1], ia[1:]
        live = denom > 0
        for _ in range(sweeps):
            change = 0.0
            for i in order:
                if not live[i]:
                    continue
                k = slice(starts[i], ends[i])
                f, j = flowja[k], ja[k]
                up = (f > 0) & (j != i)
                v = (mass_in[i] + float(np.dot(f[up], c[j[up]]))) / denom[i]
                change = max(change, abs(v - c[i]))
                c[i] = v
            if change <= tol * max(float(np.abs(c).max()), 1e-30):
                break
        return c
