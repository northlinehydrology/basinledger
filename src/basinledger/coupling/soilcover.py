"""The soil cover SWAT+'s soil evaporation sees: a process substitution (docs/design.md, Soil cover for soil
evaporation). Approved by Andrew on Nov 27, 2026, after the Rappahannock calibration (the Virginia study's SWAT+ issue 5).

SWAT+ 62.0.1's ``et_act.f90`` takes an HRU's soil cover as all its above-ground biomass plus its residue, and its
potential soil evaporation as PET x exp(-5e-5 x (cover + 0.1)). A forest's wood counts as cover, so its soil hardly
evaporates all year, while a dormant grass (LAI 0) evaporates near PET in winter. In the Rappahannock calibration the
observed forest-minus-open ET contrast was +93 mm a year and the model gave +16.

SWAT+ has no marker for a tree, so a config names the HRUs whose plants are trees, directly or by land use (the
Rappahannock's forest land uses each have a one-plant community: frsd, frse, frst, wetw). On those HRUs, before each
day, Basin Ledger sets the cover (swatplus-bmi's ``soil_surface_cover__mass-per-area_density``, which replaces
``et_act``'s ``cover`` for that day only):

- ``rule = "leaves"``: the leaf mass + the residue, no wood;
- ``rule = "swat"``: the above-ground mass + the residue, SWAT+'s own cover set back: a check, the same run.

Every other HRU is left to SWAT+. The masses are read at the pause before the day, so they are the previous day's end,
and summed as SWAT+ sums them (kg/ha, 4-byte reals). On a run's first day SWAT+'s own cover is the residue alone: it
sums the plants' masses after ``et_act`` (``pl_community``), and the masses read before that day agree.
``et_act`` uses the same masses unless the HRU's decision-table actions (management.sch's auto operations, with a
fertilizer they schedule) change them earlier that day; there the cover is a day behind. An HRU with no management
schedule has none, as forests usually don't.

Only soil evaporation sees the set cover: SWAT+'s albedo and erosion compute their own (albedo.f90, ero_cfactor.f90).
Under snow (0.5 mm or more) SWAT+ takes half of PET whatever the cover. Off unless a config has ``[soil_cover]``, so a
run without it is unchanged.

The pieces here: ``Settings`` (the config's section, checked on read), ``selected`` (the HRUs it names, by name or
land use, refusing names the project doesn't have) and ``SoilCover`` (the surface model with the cover set each day).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

RULES = ("leaves", "swat")
COVER = "soil_surface_cover__mass-per-area_density"           # settable: et_act's cover for the next day, kg m-2
LEAF = "land_vegetation_leaf__mass-per-area_density"          # SWAT+'s leaf_com
ABOVE = "land_vegetation__above-ground_mass-per-area_density"  # SWAT+'s ab_gr_com: leaves, stems and wood
RESIDUE = "land_surface_residue__mass-per-area_density"       # SWAT+'s rsd_tot


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


def _has_cover(bmi):
    names = set(bmi.get_input_var_names()) | set(bmi.get_output_var_names())
    missing = [v for v in (COVER, LEAF, ABOVE, RESIDUE) if v not in names]
    if missing:
        import swatplus_bmi
        raise ValueError(f"[soil_cover] needs swatplus-bmi's soil cover ({', '.join(missing)}), which swatplus-bmi "
                         f"{swatplus_bmi.__version__} doesn't have")


def wrap(bmi, s: Settings, swat_dir) -> SoilCover:
    """SWAT+ (started, from ``swat_dir``) with the cover set on the HRUs ``s`` names."""
    from ..salt import _hru_data
    _has_cover(bmi)
    hrus = list(bmi.get_grid_object_names(bmi.get_var_grid(COVER)))
    return SoilCover(bmi, s.rule, selected(s, hrus, _hru_data(swat_dir, hrus)["landuse"]))


class SoilCover:
    """The surface model with the soil cover set on ``hrus`` (indices) before each day, by ``rule``; everything else
    goes to the wrapped model. ``days`` counts the days set."""

    def __init__(self, bmi, rule: str, hrus: np.ndarray):
        _has_cover(bmi)
        if rule not in RULES:
            raise ValueError(f"[soil_cover] rule is one of {RULES}")
        self._m, self.rule, self.hrus = bmi, rule, np.asarray(hrus, dtype=int)
        n = bmi.get_grid_size(bmi.get_var_grid(COVER))
        self._bufs = {v: np.empty(n, dtype=bmi.get_var_type(v)) for v in (LEAF, ABOVE, RESIDUE)}
        self.days = 0

    def cover(self) -> np.ndarray:
        """The cover the rule gives the selected HRUs from the masses now (the end of the previous day), kg m-2.

        Summed as SWAT+ sums its own: in kg/ha, in 4-byte reals. Each mass read back to kg/ha is SWAT+'s own number
        exactly, and their float32 sum is rounded as SWAT+'s is (to even on a tie), which a sum in float64 isn't, so
        ``"swat"`` gives SWAT+'s cover to the bit. The sum is representable in 4 bytes, so it survives the kg m-2
        it is set in."""
        def kg_ha(v):
            self._m.get_value(v, self._bufs[v])
            return (self._bufs[v][self.hrus].astype(np.float64) * 1e4).astype(np.float32)
        total = (kg_ha(LEAF) if self.rule == "leaves" else kg_ha(ABOVE)) + kg_ha(RESIDUE)
        return total.astype(np.float64) * 1e-4

    def update(self):
        if len(self.hrus):
            self._m.set_value_at_indices(COVER, self.hrus, self.cover())
            self.days += 1
        self._m.update()

    def __getattr__(self, name):
        return getattr(self._m, name)
