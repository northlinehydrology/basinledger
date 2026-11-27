"""The soil cover SWAT+'s soil evaporation sees: a process substitution (docs/design.md, Soil cover for soil
evaporation). Approved by Andrew on Nov 27, 2026, after the Rappahannock calibration (the Virginia study's SWAT+ issue 5).

SWAT+ 62.0.1's ``et_act.f90`` takes an HRU's soil cover as all its above-ground biomass plus its residue, and its
potential soil evaporation as PET x exp(-5e-5 x (cover + 0.1)). A forest's wood counts as cover, so its soil hardly
evaporates all year, while a dormant grass (LAI 0) evaporates near PET in winter. In the Rappahannock calibration the
observed forest-minus-open ET contrast was +93 mm a year and the model gave +16.

``rule = "leaves"``: for a tree, cover = its leaf mass + the residue (no wood); for any other plant, SWAT+'s own.
It applies to the HRUs and land uses a config names, and only to soil evaporation: SWAT+'s albedo and erosion compute
their own cover. Off unless a config has ``[soil_cover]``, so a run without it is unchanged.

The pieces here: ``Settings`` (the config's section, checked on read) and ``selected`` (the HRUs it names, by name or
land use, refusing names the project doesn't have).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

RULES = ("leaves",)


@dataclass(frozen=True)
class Settings:
    rule: str
    hrus: tuple[str, ...] = ()
    land_uses: tuple[str, ...] = ()


def settings(where: str, x: dict) -> Settings:
    """``[soil_cover]``: ``rule`` (required) and ``hrus`` and/or ``land_uses`` (at least one)."""
    if not isinstance(x, dict):
        raise ValueError(f"{where}: [soil_cover] is a table: rule, and hrus and/or land_uses")
    unknown = sorted(set(x) - {"rule", "hrus", "land_uses"})
    if unknown:
        raise ValueError(f"{where}: [soil_cover] unknown settings {unknown}")
    if x.get("rule") not in RULES:
        raise ValueError(f"{where}: [soil_cover] rule is one of {RULES} (required, no default)")
    out = {}
    for key in ("hrus", "land_uses"):
        v = x.get(key, [])
        if not isinstance(v, list) or not all(isinstance(s, str) and s for s in v):
            raise ValueError(f"{where}: [soil_cover] {key} is a list of names; got {v!r}")
        if len(set(v)) != len(v):
            raise ValueError(f"{where}: [soil_cover] {key} names one twice")
        out[key] = tuple(v)
    if not out["hrus"] and not out["land_uses"]:
        raise ValueError(f"{where}: [soil_cover] names no HRUs or land uses: give hrus, land_uses or both")
    return Settings(rule=x["rule"], **out)


def selected(s: Settings, hrus: list[str], landuse: list[str]) -> np.ndarray:
    """Indices into ``hrus`` (swatplus-bmi's HRU order) of the HRUs ``s`` names, directly or by their land use
    (``landuse``, each HRU's ``lu_mgt`` in ``hru-data.hru``). A name the project doesn't have is refused."""
    missing = sorted(set(s.hrus) - set(hrus))
    if missing:
        raise ValueError(f"[soil_cover] hrus {missing[:5]} aren't in the SWAT+ project")
    missing = sorted(set(s.land_uses) - set(landuse))
    if missing:
        raise ValueError(f"[soil_cover] land_uses {missing[:5]} aren't used by any HRU in the SWAT+ project "
                         f"(its land uses: {sorted(set(landuse))[:10]})")
    named, uses = set(s.hrus), set(s.land_uses)
    return np.array([i for i, (h, lu) in enumerate(zip(hrus, landuse)) if h in named or lu in uses], dtype=int)
