"""Operations of systems and reservoirs (docs/design.md, "Operating rules for systems and reservoirs"): transfers,
drought restrictions on a system's classes, and the reservoirs' normal pools the rules read.

**A transfer** (``[[transfers]]``) moves water from a SWAT+ intake (on a channel or a reservoir) to a point source
(into a reservoir or a channel), an outfall or offsite, at the rate the rules set for the day (``rules.py``:
``pumped_storage``, ``release``; 0 when none does), up to its capacity. The intake's delivery is known only after
SWAT+'s day and a point source is set before it, so what was delivered goes into ``to`` the next day, held overnight
(carried by a restart). Under ``transfer:<name>:`` the ledger books ``rate`` (m3 s-1), ``requested`` = ``delivered``
+ ``shortfall`` and ``delivered`` = ``sent`` + ``held_gain``.

**A restriction** (rule ``restriction``) cuts a system's customer classes by share. Each class's demand asked is its
demand times (1 - cut); under ``utility:<name>:`` the ledger books ``wanted:<class>`` = ``demand:<class>`` +
``restricted:<class>``, and ``restricted`` for the day, apart from what the sources didn't deliver (``unmet``).
"""

from __future__ import annotations

from pathlib import Path

from .units import DAY

FLOAT32 = 1e-6


def pools(project: str | Path) -> dict[str, float]:
    """Each reservoir's normal pool, m3, by its name in reservoir.con (SWAT+'s, as BMI reports it): its principal
    spillway volume (hydrology.res ``vol_ps``, ha-m), through reservoir.con's ``res`` row of reservoir.res and that
    row's hydrology. Empty for a project without reservoirs."""
    project = Path(project)
    if not all((project / f).exists() for f in ("reservoir.con", "reservoir.res", "hydrology.res")):
        return {}

    def table(name):
        lines = [ln.split() for ln in (project / name).read_text().splitlines()[1:] if ln.split()]
        head = lines[0]
        return [dict(zip(head, row)) for row in lines[1:]]

    hyd = {r["name"]: float(r["vol_ps"]) * 1e4 for r in table("hydrology.res")}
    res = table("reservoir.res")
    out = {}
    for c in table("reservoir.con"):
        i = int(c["res"]) - 1
        if 0 <= i < len(res) and res[i].get("hyd") in hyd:
            out[c["name"]] = hyd[res[i]["hyd"]]
    return out


class Transfer:
    def __init__(self, name: str, frm: str, to: str, capacity: float, held: float = 0.0):
        self.name, self.frm, self.to, self.capacity = name, frm, to, float(capacity)
        self.held = float(held)  # delivered yesterday, into ``to`` today
        self._req = self._sent = 0.0

    def declare(self, ledger):
        k = f"transfer:{self.name}:"
        ledger.pair(k + "requested", [k + "delivered", k + "shortfall"], rtol=FLOAT32)
        ledger.pair(k + "delivered", [k + "sent", k + "held_gain"])

    def before(self, ledger, router, rate: float):
        """Before the models' day: the day's withdrawal at ``rate`` (m3 s-1, up to the capacity), and yesterday's
        delivery into ``to``."""
        k = f"transfer:{self.name}:"
        rate = min(max(rate, 0.0), self.capacity)
        self._req, self._sent = rate * DAY, self.held
        ledger.enter(k + "rate", rate)
        ledger.enter(k + "requested", self._req)
        ledger.enter(k + "sent", self._sent)
        router.request(self.frm, self._req)
        router.request(self.to, self._sent)

    def after_surface(self, ledger):
        """After SWAT+'s day: its part of what the intake delivered, held for tomorrow."""
        k = f"transfer:{self.name}:"
        asked = ledger.today(f"route:{self.frm}:requested")
        got = self._req * min(ledger.today(f"route:{self.frm}:delivered") / asked, 1.0) if asked > 0 else 0.0
        ledger.enter(k + "delivered", got)
        ledger.enter(k + "shortfall", self._req - got)
        ledger.enter(k + "held_gain", got - self._sent)
        self.held = got


class Restricted:
    """A system's class of customers whose demand the day's restriction cuts by ``cut`` (a share)."""

    def __init__(self, base):
        self._base = base
        for a in ("name", "sewered", "septic", "consumed", "septic_to"):
            setattr(self, a, getattr(base, a))
        self.cut = 0.0
        self.today = (0.0, 0.0)  # (wanted, restricted), m3, the day planned

    def __getattr__(self, a):  # anything else the class carries (its mass, say) is the class's own
        return getattr(self.__dict__["_base"], a)

    def on(self, day) -> float:
        wanted = self._base.on(day)
        asked = wanted * (1.0 - self.cut)
        self.today = (wanted, wanted - asked)
        return asked


def restrict(chain) -> dict[str, Restricted]:
    """The system's classes, each behind a restriction (0 until a rule cuts it); returns them by name."""
    out = {}
    for i, c in enumerate(chain.customers):
        if not isinstance(c, Restricted):
            chain.customers[i] = c = Restricted(c)
        out[c.name] = c
    return out


def declare_restricted(ledger, system: str, classes):
    k = f"utility:{system}:"
    for c in classes:
        ledger.pair(k + f"wanted:{c}", [k + f"demand:{c}", k + f"restricted:{c}"])


def book_restricted(ledger, system: str, classes: dict[str, Restricted]):
    """After the system's plan: what each class wanted and what the restriction held back."""
    k = f"utility:{system}:"
    total = 0.0
    for name, c in classes.items():
        wanted, held_back = c.today
        ledger.enter(k + f"wanted:{name}", wanted)
        ledger.enter(k + f"restricted:{name}", held_back)
        total += held_back
    ledger.enter(k + "restricted", total)
