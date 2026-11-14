# Prior art: coupling SWAT and SWAT+ with MODFLOW

What earlier couplings exchange, how often, how they map SWAT units to grid cells, and whether their code is public.
Compiled Nov 4, 2026; SWATGenX added and Basin Ledger's row brought up to 0.7.0 on Nov 14, 2026. Where a detail comes from a tutorial, code or a later paper rather than the paper itself, it says
so. Details not confirmed are marked.

| Work | Models | Exchanged | How often | Units to cells | Code |
|---|---|---|---|---|---|
| Hughes et al. 2022 | MODFLOW 6 API (BMI + XMI) | anything in MODFLOW 6's memory | per time step or per outer iteration | n/a | MODFLOW 6, xmipy, modflowapi (public) |
| Bailey et al. 2016, SWAT-MODFLOW | SWAT 2012 + MODFLOW-NWT, one executable | recharge; stream-aquifer exchange (RIV) to subbasin channels; water table to HRUs; later pumping and EVT | daily | HRUs split into located parts (DHRUs), intersected with the grid, by area | public (swat.tamu.edu) |
| Bailey et al. 2020, SWAT+ gwflow | SWAT+ with its own gridded groundwater module | recharge, lateral flow, groundwater ET, pumping, streambed exchange, saturation excess | daily | HRU-cell overlap areas; channel cells by intersection | inside SWAT+ |
| Bailey et al. 2025, SWAT+MODFLOW (GMD 18:5681) | SWAT+ rev 61 + MODFLOW-NWT, one executable | recharge and soil transfer, groundwater ET, channel and canal exchange, tile drains, reservoirs, irrigation pumping | daily | HRU recharge × fraction of HRU in cell; channel length per cell | public (Zenodo) |
| Nikiel and Żurek 2026 (EGU26-22367) | SWAT+ + MODFLOW 6, through Python scripts | recharge and river flow, both ways | not stated | "upscaling" HRUs to the grid; method not described | none found |
| SWATGenX, SWAT+ with MODFLOW 6 (Rafiei, 2026) | a SWAT+ fork calling MODFLOW 6 (libmf6, XMI), one program | recharge; groundwater discharge to channels (SFR); PFAS in soil, aquifer and channels | daily | HRU to cell by area; reaches to channels by map file | engine public (LGPL-2.1); model generator not public |
| **Basin Ledger** | SWAT+ 62.0.1 (swatplus-bmi) + MODFLOW 6 (libmf6) + water-user BMI components, separate programs | recharge (per aquifer or per HRU); well pumping; stream-aquifer exchange both ways (RIV with a supply limit, or SFR); groundwater up into the soil | daily, discharge returned the next day, with a ledger check each day | area maps (aquifer or HRU to cell); flowpath lengths per cell | this repository |

## Hughes et al. 2022: the MODFLOW API

Hughes, J.D., Russcher, M.J., Langevin, C.D., Morway, E.D., McDonald, R.R., 2022. The MODFLOW Application Programming
Interface for simulation control and software interoperability. *Environmental Modelling & Software* 148, 105257.
https://doi.org/10.1016/j.envsoft.2021.105257

- **What it is.** MODFLOW 6 built as a shared library (`libmf6`) that implements BMI and an extension, XMI. XMI allows
  control inside a time step, which tight coupling needs.
- **Time step.** `prepare_time_step(dt)`, `do_time_step()` and `finalize_time_step()`.
- **Solve.** `get_subcomponent_count()`, then `prepare_solve(i)`, `solve(i)` (returns whether it converged) and
  `finalize_solve(i)`, so a caller can change values between outer iterations.
- **Variable access.** `get_var_address` and `get_value_ptr` give live pointers to internal arrays: a package's
  rates, `HCOF` and `RHS`, and heads (`X`). `set_value` writes to them.
- **[modflowapi](https://github.com/MODFLOW-ORG/modflowapi)** (USGS, CC0) wraps this for Python:
  - `ModflowApi`, a subclass of Deltares' [xmipy](https://github.com/Deltares/xmipy);
  - `run_simulation(dll, sim_path, callback)`, which calls back at the start and end of the simulation, each stress
    period, each time step and each iteration;
  - `ApiSimulation`, `ApiModel` and `ApiPackage` objects for reading and writing package data by name.
- **What Basin Ledger uses.** `ModflowApi`, driven step by step rather than through the callback, because
  Basin Ledger's driver controls time for every model. One detail found while building it: in MODFLOW 6.5 and later the
  values to set are the named input arrays (`RCH/RECHARGE`, `WEL/Q`), not `BOUND` (see `docs/design.md`).

## Bailey et al. 2016: SWAT-MODFLOW

Bailey, R.T., Wible, T.C., Arabi, M., Records, R.M., Ditty, J., 2016. Assessing regional-scale spatio-temporal patterns
of groundwater–surface water interactions using a coupled SWAT-MODFLOW model. *Hydrological Processes* 30(23),
4420–4433. https://doi.org/10.1002/hyp.10933 (Sprague River, Oregon, 1970–2003).

- **Exchanged:**
  - HRU deep percolation becomes MODFLOW recharge;
  - MODFLOW's stream-aquifer exchange (RIV) goes back to subbasin channels;
  - the water table goes back to HRUs.
  - Later versions added pumping for irrigation HRUs and groundwater ET. That is from the tutorial and later papers;
    whether the 2016 version had them is not confirmed.
- **How often:** daily. This comes from the account of the lineage in Bailey et al. 2025, not confirmed in the 2016
  text.
- **Mapping:** HRUs, which in SWAT have no location, are split into located polygons (DHRUs), intersected with the
  grid and weighted by area (`dhru_grid` and `grid_dhru` files). River cells come from intersecting the grid with
  subbasin channels. Source: the
  [SWAT-MODFLOW tutorial](https://swat.tamu.edu/media/115048/swat-modflow-tutorial.pdf).
- **Code:** https://swat.tamu.edu/software/swat-modflow (source and executables, versions 1–3). MODFLOW-NWT is compiled
  into SWAT as one program.
- **Tools:** QSWATMOD (Park et al. 2019, *Environ. Model. Softw.* 111, 493–497,
  https://doi.org/10.1016/j.envsoft.2018.10.017) and SWATMOD-Prep (https://doi.org/10.1111/1752-1688.12502).

## Bailey et al. 2020: SWAT+ gwflow

Bailey, R.T., Bieger, K., Arnold, J.G., Bosch, D.D., 2020. A new physically-based spatially-distributed groundwater
flow module for SWAT+. *Hydrology* 7(4), 75. https://doi.org/10.3390/hydrology7040075 (Little River, Georgia).

- **What it replaces:** SWAT+'s lumped aquifer, with a water balance per grid cell. It is its own solver, not MODFLOW.
- **Exchanged:**
  - HRU recharge;
  - lateral flow between cells;
  - groundwater ET and pumping;
  - streambed exchange and saturation-excess flow.
  - Later versions add tile drains, groundwater-to-soil transfer, reservoirs and floodplains.
- **How often:** daily, inside SWAT+'s daily loop.
- **Mapping:** HRU-cell overlap areas (`gwflow.hrucell`); channel cells from intersecting the grid with SWAT+'s
  channels (`gwflow.rivcells`).
- **Code:** part of SWAT+.
- **Later work:** Yimer, E.A., Bailey, R.T., et al., 2023. Regional evaluation of groundwater–surface water
  interactions using a coupled geohydrological model (SWAT+gwflow). *J. Hydrol.: Reg. Stud.* 50, 101532.
  https://doi.org/10.1016/j.ejrh.2023.101532

**For Basin Ledger:** gwflow is the simpler choice when a study needs only SWAT+ with a gridded aquifer. MODFLOW 6 is
the choice when a study needs a layered aquifer, a regional model that already exists (most utilities' well fields
have one), or MODFLOW's packages: wells screened in several layers, SFR, transport. swatplus-bmi's saved states don't
yet cover gwflow (swatplus-bmi README).

## Bailey et al. 2025: SWAT+MODFLOW (GMD 18:5681)

Bailey, R.T., Abbas, S., Arnold, J.G., White, M.J., 2025. SWAT+MODFLOW: a new hydrologic model for simulating
surface–subsurface flow in managed watersheds. *Geoscientific Model Development* 18, 5681–5697.
https://doi.org/10.5194/gmd-18-5681-2025 (discussion: egusphere-2025-300).

- **Exchanged:**
  - recharge and soil transfer (RCH);
  - groundwater ET (EVT);
  - channel and canal exchange (RIV);
  - tile drains (DRN);
  - reservoirs;
  - irrigation pumping (WEL).
- **How often:** daily.
- **Mapping:** links are built beforehand in GIS and read from files.
  - Recharge per cell is the sum of HRU recharge × the fraction of that HRU's area in the cell.
  - Channel exchange uses each channel's length in each cell, with one stage per channel.
  - Each drain cell is assigned to its nearest channel.
  - Pumping comes from the cell under the field, or from cells linked to the HRU with the demand.
  - The study grids use 500 m cells and one unconfined layer, and SWAT+'s aquifer is switched off.
- **Code:** https://doi.org/10.5281/zenodo.14674981. MODFLOW-NWT is compiled into SWAT+ revision 61, with MODFLOW
  called as a subroutine.

**For Basin Ledger:** this paper sets the mapping to follow: HRU recharge to cells by area, channels to cells by
length. Basin Ledger differs in keeping the models separate, using MODFLOW 6, and checking every exchange with a
ledger.

## Nikiel and Żurek 2026 (EGU26-22367)

Nikiel, M., Żurek, A.J., 2026. Coupling SWAT+ and MODFLOW 6 using Python: a flexible approach to representing
groundwater–surface water interactions in an agricultural catchment. EGU General Assembly 2026, EGU26-22367.
https://doi.org/10.5194/egusphere-egu26-22367 (AGH Kraków).

- **Exchanged:** groundwater recharge and river flow, both ways, passed by Python scripts.
- **How often:** not stated.
- **Mapping:** the abstract names "upscaling" irregular HRUs to the grid as the main difficulty and does not describe a
  method.
- **Model:** MODFLOW 6 on a DIS grid, migrated from a MODFLOW-96 model of a catchment of about 250 km² in south-central
  Poland. An earlier abstract by the group covers data preparation with flopy (EGU25-12133,
  https://doi.org/10.5194/egusphere-egu25-12133).
- **Code:** none found. Only an abstract so far.

**For Basin Ledger:** close prior work, and the only SWAT+ to MODFLOW 6 coupling through Python scripts found. It is
worth contacting the authors.

## SWATGenX: SWAT+ with MODFLOW 6 (Rafiei, 2026)

SWATGenX LLC (Michigan); the engine is by Vahid Rafiei. Sources, read Nov 14, 2026:
- the coupling page, https://www.swatgenx.com/swat-modflow-6;
- the engine, https://github.com/rafiei-vahid/swatplus (a SWAT+ fork, LGPL-2.1);
- the PFAS study repository, https://github.com/SWATGenX/swatplus-mf6-pfas (created June 2026);
- the SWAT+ model generator's announcement on the SWAT+ user group, March 2026.

Related paper, not read: Rafiei, V., Nejadhashemi, A.P., 2023. *Water Research* 240, 120073 (watershed-scale PFAS
fate and transport).

- **What it is.** SWAT+ changed to load `libmf6` and call MODFLOW 6 through XMI, as one program. An `mf6.con` file
  switches the coupling on: the MODFLOW 6 workspace, the step in whole days, a recharge multiplier and the library
  path. Without it the engine runs as plain SWAT+.
- **Exchanged:**
  - SWAT+ deep percolation becomes MODFLOW 6 recharge;
  - MODFLOW 6 groundwater discharge (SFR) goes back to SWAT+ channels;
  - PFAS in the soil, the aquifer (GWT) and the channels. Passing soil PFAS down into MODFLOW 6 transport is built but,
    per the page, has not yet carried PFAS in a run.
- **How often:** daily. The order within the day is not stated.
- **Mapping:** HRU to cell by area (`mf6_recharge.map`); SFR reaches to SWAT+ channels (`mf6_baseflow.map`).
- **Automatic builds:** SWAT+ models for USGS gages (over 600, self-reported). The groundwater model is built
  steady-state from annual-mean percolation; the page says automatic builds do not yet run the two-way coupling.
- **Checks published:** MODFLOW 6's percent discrepancy, and simulated against observed heads at wells. Recharge and
  hydraulic conductivity are calibrated to heads.
- **Code:** the engine is public; the model generator is not.

**For Basin Ledger:** the nearest work in its parts: the same MODFLOW 6 interface, HRU recharge by area, and SFR back to
channels. The design differs. SWATGenX changes SWAT+ so that it calls MODFLOW 6 inside one program. Basin Ledger keeps
stock SWAT+ (through swatplus-bmi) and MODFLOW 6 as separate components under its own driver, links water-user
components the same way, and stops a run when an exchange does not close. SWATGenX was not used in Basin Ledger's
design; it is listed here for comparison.

## A groundwater model for a Michigan example

Feinstein, D.T., Hunt, R.J., Reeves, H.W., 2010. Regional groundwater-flow model of the Lake Michigan Basin in support
of Great Lakes Basin water availability and use studies. USGS Scientific Investigations Report 2010-5109.
https://pubs.usgs.gov/sir/2010/5109/

- **Code:** SEAWAT-2000 (variable density, built on MODFLOW-2000), with a uniform-density MODFLOW-2000 version for
  comparison. Calibrated with PEST.
- **Grid:** 391 rows × 261 columns × 20 layers, about 2 million cells. Inner cells are 5,000 ft (about 1.5 km); outer
  cells grow by a factor of 1.3.
- **Period:** transient, 1864–2005.
- **Files:**
  - A data release for the 2010 model itself was not found. Whether one exists is not confirmed.
  - MODFLOW-NWT inset models derived from it are archived (Haserodt, Kauffman, Feinstein): https://doi.org/10.5066/F76D5R5V,
    with SIR 2018-5038 (https://pubs.usgs.gov/publication/sir20185038).
  - A local grid refinement study is in SIR 2010-5117 (https://pubs.usgs.gov/publication/sir20105117).
- **MODFLOW 6:** no MODFLOW 6 version of the regional model was found. USGS's nearest MODFLOW 6 product is a 2024 set
  of 75 steady-state models of the surficial aquifer for the contiguous U.S.
  ([data release](https://www.usgs.gov/data/modflow-6-models-used-simulate-long-term-average-surficial-groundwater-system-contiguous)).
- **Michigan's Water Withdrawal Assessment Tool** estimates streamflow depletion analytically, not with MODFLOW. It
  was checked against a MODFLOW model of Kalamazoo County (Reeves, H.W., Hamilton, D.A., Seelbach, P.W., Asher, A.J.,
  2009, USGS SIR 2009-5003, https://pubs.usgs.gov/sir/2009/5003).

**For a Michigan example:** start from the MODFLOW-NWT insets where they cover the basin, or ask USGS's Upper Midwest
Water Science Center for the regional model's files. Convert to MODFLOW 6 with flopy, refine around the study basin,
and take boundary heads from the regional model.
