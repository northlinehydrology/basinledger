"""What Basin Ledger takes over in a MODFLOW 6 model, and what it leaves as the model's (docs/design.md, "Extending an
existing MODFLOW 6 model").

A model someone else built comes with its own boundaries: recharge from its own estimate, rivers, drains, general
and fixed heads, its wells, evapotranspiration. Basin Ledger sets only what it couples:

- recharge, in the package ``[recharge] package`` names (every entry, from SWAT+);
- the stream exchange, on the entries of ``[exchange] package`` (and its drains) that the map gives a SWAT+ channel;
- groundwater up into the soil, in ``[groundwater_to_soil] package``;
- users' wells, the entries ``[[wells]]`` names.

Everything else stays the model's, and the ledger books its flow each day as water from or to outside the coupled
models (``groundwater:boundary:<package>``), as it always has. The start-up check lists each package with what
it was taken over for and how many of its entries, so a reviewer sees what the coupling changed. It also refuses
what would count water twice:

- another recharge package left as the model's (RCH): the model's own recharge would be added to SWAT+'s;
- a UZF package: it takes infiltration of its own, which is recharge again (coupling UZF is not built);
- an EVT package with groundwater up into the soil: SWAT+'s plants then take groundwater from the soil, and EVT
  would take it a second time.

And it notes, without refusing, an EVT package without that link (groundwater ET stays the model's), and an SFR or
LAK package left as the model's while SWAT+'s channels are coupled: their water is the model's, so if they stand for
the same streams or lakes as SWAT+'s channels, the river is carried twice.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


class TakeoverNote(UserWarning):
    """Something about the model's own packages a reviewer should know: a start-up warning, listed by a check."""


@dataclass
class PackageRole:
    package: str         # the package's name in the model, lower case
    type: str            # its type, lower case: rch, riv, drn, ghb, chd, wel, evt, sfr, lak, maw, uzf, ...
    entries: int         # its entries this period (reaches for SFR, connections for LAK and MAW)
    taken: int           # how many Basin Ledger sets
    what: str            # what Basin Ledger takes it over for, or "" (left as the model's)

    @property
    def role(self) -> str:
        if self.taken == 0:
            return "the model's"
        if self.taken == self.entries:
            return "Basin Ledger's"
        return "partly Basin Ledger's"

    def describe(self) -> str:
        left = f"booked as groundwater:boundary:{self.package}"
        if self.taken == 0:
            return f"the model's: water from or to outside the coupled models, {left}"
        if self.taken == self.entries:
            return f"Basin Ledger's: {self.what}"
        return (f"{self.taken} of {self.entries} entries Basin Ledger's ({self.what}); the other "
                f"{self.entries - self.taken} the model's, {left}")


def flow_packages(mf6, model: str | None = None) -> list[str]:
    """A groundwater model's packages with flows (SIMVALS), upper case, sorted: ``model``, or the one Basin Ledger
    couples."""
    prefix = (model or mf6.model).upper() + "/"
    names = set(mf6.api.get_input_var_names()) | set(mf6.api.get_output_var_names())
    return sorted({n.split("/")[1] for n in names if n.startswith(prefix) and n.endswith("/SIMVALS")})


def inventory(mf6, model: str | None = None) -> list[PackageRole]:
    """Each package with flows in a model as given, before anything is coupled: what a first look at someone's model
    lists (call after ``prepare_day``, when the packages have their entries). ``model``: any groundwater model of the
    simulation, the one ``mf6`` couples if not given."""
    from .exchange import package_type
    m = (model or mf6.model).upper()
    out = []
    for p in flow_packages(mf6, m):
        n = int(mf6.get("NBOUND", m, p)[0])
        t = str(mf6.get("PACKAGE_TYPE", m, p)[0]).strip().lower()
        out.append(PackageRole(p.lower(), t, n, 0, ""))
    return out


def roles(mf6, coupled: dict[str, tuple[np.ndarray, list[str]]]) -> list[PackageRole]:
    """Each package with flows (SIMVALS), in the model's order of names: ``coupled`` maps a package's upper-case name to
    (a mask of its entries Basin Ledger sets, the links that set them)."""
    from .exchange import package_type
    out = []
    for p, (mask, links) in coupled.items():
        out.append(PackageRole(p.lower(), package_type(mf6, p), int(len(mask)), int(mask.sum()),
                               ", ".join(dict.fromkeys(links))))
    return out


def problems(rows: list[PackageRole], groundwater_to_soil: bool) -> tuple[list[str], list[str]]:
    """What the model's own packages would count twice (problems), and what a reviewer should know (notes)."""
    bad, notes = [], []
    for r in rows:
        if r.taken:
            continue
        if r.type == "rch":
            bad.append(f"{r.package}: a recharge package left as the model's; its recharge would be added to SWAT+'s. "
                       "Name it in [recharge] package if it is the one SWAT+'s recharge replaces, or leave it out of "
                       "the run's copy of the model")
        elif r.type == "uzf":
            bad.append(f"{r.package}: UZF takes infiltration of its own, which is recharge again; coupling UZF to "
                       "SWAT+'s percolation isn't built")
        elif r.type == "evt" and groundwater_to_soil:
            bad.append(f"{r.package}: EVT takes groundwater ET, and with [groundwater_to_soil] SWAT+'s plants take it "
                       "from the soil too; leave EVT out of the run's copy of the model")
        elif r.type == "evt":
            notes.append(f"{r.package}: groundwater ET stays the model's (EVT), apart from SWAT+'s ET from the soil")
        elif r.type in ("sfr", "lak"):
            notes.append(f"{r.package}: the model's {r.type.upper()} carries water SWAT+ doesn't see; if it stands for "
                         "the same streams or lakes as SWAT+'s channels, that water is carried twice")
    return bad, notes


def table(rows: list[PackageRole]) -> str:
    """The start-up check's list, one package per line."""
    w = max([len(r.package) for r in rows] + [7])
    lines = [f"{'package':<{w}}  type  {'entries':>7}  taken over, or left as the model's"]
    for r in rows:
        lines.append(f"{r.package:<{w}}  {r.type:<4}  {r.entries:>7}  {r.describe()}")
    return "\n".join(lines)
