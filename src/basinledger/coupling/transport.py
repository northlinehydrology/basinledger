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
splits it as the water ledger does (recharge, wells, coupled streams, boundaries) and adds the release from mobile
and sorbed storage (MST's rates) and GWT's own boundary packages (a constant concentration, a mass source). The sum
closes to the transport solve's convergence; ``closure_share`` is the residual over the day's gross mass flow.

Concentrations in kg m-3 (= g L-1), masses in kg, with MODFLOW 6 in meters and days.
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


def add_transport(sim, gwf, species: str, porosity, dispersivity, sorption: dict | None = None, strt=0.0,
                  sources: list[tuple[str, str]] = (), scheme: str = "tvd", diffusion: float = 0.0,
                  ims: dict | None = None, save: str = "LAST"):
    """A GWT model ``gwt_<species>`` on ``gwf``'s grid (DIS only), linked to it and solved after it.

    porosity: effective porosity [A] (scalar or per cell). dispersivity: m [A], see ``_dispersivities``.
    sorption: None (a conservative species, such as chloride), or ``{"distcoef": Kd (m3 kg-1), "bulk_density": kg
    m-3}`` for a linear isotherm [A]. strt: starting concentration, kg m-3. sources: ``(flow package pname, auxiliary
    name)`` for each package whose inflow carries a concentration. scheme: ADV's (``upstream``, ``central`` or
    ``tvd``). diffusion: molecular diffusion coefficient, m2 d-1 [A]. ims: settings for the transport solution, over
    the defaults below. save: when OC saves concentrations and the budget (``LAST`` or ``ALL``)."""
    import flopy

    name = f"gwt_{species}"
    dis = gwf.get_package("dis")
    if dis is None or dis.package_type != "dis":
        raise NotImplementedError("add_transport supports a structured grid (DIS) only")
    gwt = flopy.mf6.ModflowGwt(sim, modelname=name, save_flows=True)
    settings = dict(complexity="simple", linear_acceleration="bicgstab", outer_maximum=50, outer_dvclose=1e-10,
                    inner_maximum=300, inner_dvclose=1e-12, rcloserecord=[1e-10, "strict"])
    settings.update(ims or {})
    solver = flopy.mf6.ModflowIms(sim, filename=f"{name}.ims", **settings)
    sim.register_ims_package(solver, [name])
    flopy.mf6.ModflowGwtdis(gwt, length_units="meters", nlay=dis.nlay.get_data(), nrow=dis.nrow.get_data(),
                            ncol=dis.ncol.get_data(), delr=dis.delr.get_data(), delc=dis.delc.get_data(),
                            top=dis.top.get_data(), botm=dis.botm.get_data(), idomain=dis.idomain.get_data())
    mst = dict(porosity=porosity)
    if sorption:
        if set(sorption) - {"distcoef", "bulk_density"}:
            raise ValueError(f"sorption: unknown settings {sorted(set(sorption) - {'distcoef', 'bulk_density'})}")
        mst.update(sorption="linear", distcoef=sorption["distcoef"], bulk_density=sorption["bulk_density"])
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
    flopy.mf6.ModflowGwtoc(gwt, budget_filerecord=f"{name}.cbc", concentration_filerecord=f"{name}.ucn",
                           saverecord=[("CONCENTRATION", save), ("BUDGET", save)], printrecord=[("BUDGET", save)])
    flopy.mf6.ModflowGwfgwt(sim, exgtype="GWF6-GWT6", exgmnamea=gwf.name, exgmnameb=name,
                            filename=f"{gwf.name}_{name}.gwfgwt")
    return gwt


@dataclass
class MassBudget:
    """One day of a transport model's mass, kg, positive as named. ``closure`` is the residual (into the aquifer
    less out of it less what storage gained); ``closure_share`` divides it by the day's gross mass flow."""

    recharge: float = 0.0                                      # in with recharge
    wells: dict[str, float] = field(default_factory=dict)      # out with each coupled well
    to_stream: dict[str, float] = field(default_factory=dict)  # out to each coupled channel (rivers' gains, drains)
    from_stream: dict[str, float] = field(default_factory=dict)  # in from each coupled channel (rivers' losses)
    leakage_m3: dict[str, float] = field(default_factory=dict)  # the water those losses carried, m3
    boundary: dict[str, float] = field(default_factory=dict)   # net, into the aquifer, by package (not coupled)
    storage_release: float = 0.0                               # released from mobile and sorbed storage
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
    drains), with entry_channel -1 for entries not coupled. Everything else is a boundary."""

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
        if int(self._get("MST", "IDCY")[0]) != 0:
            raise NotImplementedError(f"{self.model.lower()}: decay isn't in the mass ledger yet")
        self.sorbs = int(self._get("MST", "ISRB")[0]) != 0
        # GWT's own boundary packages (a constant concentration, CNC; a mass source, SRC): every package of the
        # transport model with simulated values, mass from outside the coupled models
        prefix = self.model + "/"
        self.boundaries = sorted({v.split("/")[1] for v in names if v.startswith(prefix) and v.endswith("/SIMVALS")})

    def _get(self, package: str, var: str) -> np.ndarray:
        return self.mf6.api.get_value(self.mf6.api.get_var_address(var, self.model, package)).copy()

    def _copy(self, package: str) -> int:
        try:
            return self.copies[package.upper()]
        except KeyError:
            raise KeyError(f"{self.model.lower()}: no flow package {package!r} (it has {sorted(self.copies)})") from None

    # ---------- concentrations ----------
    def concentration(self) -> np.ndarray:
        """Concentration per reduced node, kg m-3 (GWT's X): at the end of the last day solved, or the start."""
        return self.mf6.api.get_value(self.mf6.api.get_var_address("X", self.model)).copy()

    def set_concentration(self, c: np.ndarray):
        """Concentrations per reduced node, for a restart: after initialize and before the first prepare_day."""
        addr = self.mf6.api.get_var_address("X", self.model)
        if self.mf6.api.get_value(addr).shape != np.shape(c):
            raise ValueError("saved concentrations don't match this model's nodes")
        self.mf6.api.set_value(addr, np.asarray(c, dtype=np.float64))

    def takes_inflow(self, package: str) -> bool:
        """Whether SSM gives the package's inflow a concentration (an AUX source)."""
        return int(self.iaux[self._copy(package) - 1]) > 0

    def set_inflow(self, package: str, c, entries: np.ndarray | None = None):
        """The concentration (kg m-3) of water entering through a flow package today, on all its entries or those
        selected (a boolean mask over the entries): after ``prepare_day``, before the solve."""
        n = self._copy(package)
        col = int(self.iaux[n - 1]) - 1
        if col < 0:
            raise ValueError(f"{self.model.lower()}: {package!r} isn't an AUX source in SSM, so its inflow carries no "
                             "concentration (add_transport's sources)")
        nb = int(self._get(f"FMI-FT{n}", "NBOUND")[0])
        aux = self._get(f"FMI-FT{n}", "AUXVAR")
        a = aux.reshape(-1, aux.shape[-1]) if aux.ndim > 1 else aux.reshape(-1, 1)
        sel = slice(0, nb) if entries is None else np.flatnonzero(np.asarray(entries[:nb], bool))
        a[sel, col] = c
        self.mf6.api.set_value(self.mf6.api.get_var_address("AUXVAR", self.model, f"FMI-FT{n}"), a.reshape(aux.shape))

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

    def budget(self, recharge: str | None = None, wells: dict[str, tuple[str, int]] | None = None,
               streams: list[tuple[str, np.ndarray, list[str]]] = ()) -> MassBudget:
        """The day's mass, kg, from the transport model's flows and concentrations after the solve."""
        b = MassBudget()
        net = 0.0
        x = self.concentration()
        flows = {pkg: self._entries(pkg, x) for pkg in self.copies}
        coupled: dict[str, np.ndarray] = {}
        for pkg, (q, m) in flows.items():
            net += float(m.sum())
            b.gross += float(np.abs(m).sum())
            coupled[pkg] = np.zeros(len(m), bool)
            if recharge and pkg == recharge.upper():
                b.recharge = float(m.sum())
                coupled[pkg][:] = True
        for name, (pkg, i) in (wells or {}).items():
            q, m = flows[pkg.upper()]
            b.wells[name] = -float(m[i])
            coupled[pkg.upper()][i] = True
        for pkg, entry_channel, channels in streams:
            q, m = flows[pkg.upper()]
            e = entry_channel[: len(m)]
            for i, ch in enumerate(channels):
                sel = e == i
                b.to_stream[ch] = b.to_stream.get(ch, 0.0) - float(m[sel & (m < 0)].sum())
                b.from_stream[ch] = b.from_stream.get(ch, 0.0) + float(m[sel & (m > 0)].sum())
                b.leakage_m3[ch] = b.leakage_m3.get(ch, 0.0) + float(q[sel & (q > 0)].sum())
            coupled[pkg.upper()] |= e >= 0
        for pkg, (q, m) in flows.items():
            rest = ~coupled[pkg]
            if rest.any():
                b.boundary[pkg.lower()] = float(m[rest].sum())
        for pkg in self.boundaries:  # GWT's own, kg d-1 into the aquifer
            m = self._get(pkg, "SIMVALS")[: int(self._get(pkg, "NBOUND")[0])].astype(float)
            b.boundary[pkg.lower()] = b.boundary.get(pkg.lower(), 0.0) + float(m.sum())
            net += float(m.sum())
            b.gross += float(np.abs(m).sum())
        # MST's rates, positive where storage releases mass into the flow (MODFLOW 6's budget sign)
        sto = float(self._get("MST", "RATESTO").sum())
        if self.sorbs:
            sto += float(self._get("MST", "RATESRB").sum())
        b.storage_release = sto
        b.closure = net + sto
        b.gross += abs(sto)
        return b
