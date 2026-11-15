# Changelog

## Unreleased

The second way in begins: an existing MODFLOW 6 model, taken as given (docs/design.md, "Extending an existing
MODFLOW 6 model").

- The design pass on its open questions: LAK's rain, evaporation, runoff, inflow and withdrawal, MAW's rate and SFR's
  inflows, status, stage and upstream fractions can all be set each step through libmf6 (kept as tests,
  `tests/test_existing_api.py`); how SFR routes and what driving it means; three options for who owns the water of a
  model whose SFR was calibrated as the river (not decided); candidate published models for the test.
- What Basin Ledger takes over in MODFLOW 6 (`coupling/existing.py`): at start-up each package with flows is listed
  with the entries a link sets and the rest, the model's, booked as water from or to outside the coupled models.
  `basinledger check` prints the list, and `Basin.check(..., takeover=[])` returns it. A second recharge package
  left as the model's, a UZF package, and an EVT package with groundwater up into the soil are refused: each would
  count water twice. An EVT, SFR or LAK package left as the model's is noted.
- Maps on a grid Basin Ledger didn't write (`gridmaps.py`, the `gis` extra): HRUs or aquifers to cells by area on
  each column's top active layer, with ground outside the grid, over inactive columns or in slivers reported; SWAT+
  channels to a RIV, DRN or SFR package's entries by length, with entries no channel reaches reported (or reached
  within a distance the study sets), shared cells reported, and each channel's length checked. DIS and DISV grids,
  rotated and georeferenced.
- The SFR line-up rule, the same for a package Basin Ledger built and one a study brought: a connection between
  channels is allowed from a channel's last reach into another's top reach (or from the model's own reaches into a
  top reach) and is cut each day (the top reach's upstream fraction set to 0), so SWAT+ alone carries water between
  channels. A confluence inside a channel is refused by name; `gridmaps.sfr_confluences` lists where SWAT+ needs an
  outlet. Packages built from channel lines are unchanged.
- Who owns a calibrated SFR's water (Andrew): A, SFR routing within each channel and SWAT+ between, in both paths;
  and B, a fixed stage from SWAT+'s channel depth, for channels that can't be lined up. `[exchange] sfr = "route"`
  (the default: a channel that doesn't line up is refused), `"route_or_stage"` or `"stage"`, chosen per study and
  applied per channel. A channel at a fixed stage takes status simple, its stage its streambed top plus SWAT+'s
  depth (`exchange:<channel>:sfr_depth_m`), and its water in at its entries by reach length. The start-up check
  lists which channels run as A and which as B, and why.
- The first published test model: Pleasant Lake, Wisconsin Central Sands (USGS, doi:10.5066/P9BVFSGJ;
  `examples/central-sands/fetch.py`, 16 MB, into the cache). Read as given, its packages listed
  (`existing.inventory`), its SFR networks lined up with SFRmaker's own lines as SWAT+'s channels (they line up;
  `gridmaps.sfr_aux`), its confluences listed and its mover cut at run time; what stops a coupled run is refused by
  name. Its tests skip until it is fetched.
- The second published test model: the Kalamazoo River basin, Michigan (USGS, doi:10.5066/F76D5R5V), its MODFLOW-NWT
  inset converted to MODFLOW 6 by mf5to6 (`examples/kalamazoo/fetch.py`, from a model.zip downloaded by hand). Feet
  and days, steady with no storage, 7.2 million cells, SFR, RIV, DRN, GHB, IBOUND heads and MNW2 wells as MAW.
  Its tests read the converted files and skip until it is fetched.
- Lakes follow the stream precedent (Andrew): SWAT+'s reservoir owns the water, MODFLOW 6's LAK is the face, active
  within the day from SWAT+'s stage (A) or at a fixed stage (B). Both checked through the API; not built until
  swatplus-bmi exposes reservoirs (docs/requests.md, 6).
- Large users and irrigation, MODFLOW first (Andrew): groundwater withdrawals come out of MODFLOW 6's wells. In
  two-way mode a SWAT+ project that irrigates or allocates water from its own aquifer is refused
  (`swatplus.aquifer_withdrawals`): that aquifer only holds what it passed to MODFLOW 6. Irrigation pumped from the
  model's wells needs swatplus-bmi to report HRU irrigation (docs/requests.md, 7).
- A MODFLOW 6 model in its own units, converted at each link (Andrew's choice over rescaling the model): feet,
  meters or centimeters, and seconds, minutes, hours or days, from DIS and TDIS, or stated in `[modflow6]
  length_units` and `time_units` where the model leaves them undefined (as mf5to6 does); a run that doesn't state
  them, or states them otherwise than the model, is refused. `Modflow6.get` and `set` convert every value a link
  reads or sets by its dimension (recharge, stream exchange stages, conductances and flows on RIV, DRN and SFR,
  wells, the soil link, heads, every package's flow and storage, transport concentrations and mass), so the ledger
  stays in cubic meters; a real value of unknown dimension is refused, and a model in meters and days converts
  nothing. Maps on a grid in feet come out in square meters and meters. SFR's Manning's constant is checked for the
  model's units. `basinledger check` prints the units and their factors (`Basin.check(..., units=[])`). Twins in
  feet, seconds and centimeters give their meters twin's values and, coupled to SWAT+, its ledger.
- A model's own stress periods in daily steps: a run rewrites the TDIS of its copy so every step is one day and
  every period is kept (`coupling/timesteps.py`, `daily_steps`): monthly or seasonal periods of one step, or of
  steps growing by a factor, run day by day with each period's input as given. Periods that aren't whole days end
  on the nearest day (at most 12 hours off, reported); a period shorter than that is refused; adaptive time stepping
  is left out. The start-up says what was rewritten, and notes a TDIS start date that isn't the run's. A one-period
  model of longer steps now restarts. The synthetic foreign model in five growing steps a month gives its daily
  twin's run.
- A synthetic model in someone else's style (`examples/synthetic/existing.py`, `existing.toml`): a rotated two-layer
  grid with inactive cells, monthly stress periods of daily steps, its own recharge, creek, tile drains, regional
  boundaries, city wells and riparian ET. Coupled both ways, the ledger closes every day across the periods, on DIS
  and DISV alike. It can be written in other units (`length_units`, `time_units`, and `state_units=False`
  to leave them undefined), as twins for the unit conversion.

## 0.7.0 (2026-11-11)

SFR's reaches built from SWAT+'s channel lines, and `basinledger check`. Needs swatplus-bmi 0.9.0, as 0.6.0
did.

- SFR reaches for the exchange, built from SWAT+'s channel lines (`coupling/sfr.py`): `channel_reaches` cuts each
  channel's line into reaches on a structured or vertex grid, in order downstream; `add_sfr` writes the package the
  exchange takes (reaches chained within each channel only, Manning's constant for meters and days); and
  `write_reach_map` writes the `object, reach` map. The streambed's numbers come from the study; a streambed rising
  downstream is refused, a flat one needs `min_slope`. Needs shapely: the new `gis` extra. The synthetic SFR example
  is built with it.
- `basinledger check basin.toml` (`Basin.check`): starts the models, has MODFLOW 6 read its first day, builds every
  link, and lists every problem the start-up finds, without running a day; exits with 1 if there are any. A rule's
  cell is now checked at start-up rather than on the first day it's read.
- The design note's principles: Basin Ledger builds and checks the links between models, never the models.

## 0.6.0 (2026-11-11)

The two-way exchange through SFR, and calibration parameters with PEST++ IES files. Needs swatplus-bmi 0.9.0, as
0.5.1 did.

- The two-way exchange through MODFLOW 6's streamflow routing package (`[exchange] package` an SFR package): each
  coupled channel's water goes in at its top reach, SFR computes the stages (Manning's equation) and the leakage
  inside the solve, so a reach can't lose more than it carries, and the channel's net exchange goes to SWAT+ the next
  day, as with RIV (`exchange:<channel>:sfr_inflow`). SWAT+ stays the only model routing between channels: reaches
  connected across channels, a channel with more than one top reach, diversions, the mover, SFR's storage option and
  a Manning's constant for other units than meters and days are refused at start-up, as are `limit = false`,
  `[rivers]` on the SFR package and transport with an SFR exchange. Maps can give reaches by number (`object,
  reach`). The synthetic example `sfr.toml`; restarts match the unbroken run bit for bit.
- Calibration parameters written into a run's copy of the models' files before they start
  (`calibration.apply_parameters`): SWAT+ through its own `calibration.cal` (swatplus-bmi 0.9.0, which refuses
  changes SWAT+ 62.0.1 wouldn't make as asked), with targets like `cn2:pctchg` or `awc:absval:layers=1-2`; MODFLOW 6
  as multipliers on `npf.k`, `npf.k22`, `npf.k33`, `sto.sy`, `sto.ss`, `riv.cond`, `drn.cond`, `ghb.cond`, `sfr.rhk`
  and `sfr.man`, everywhere, by layer or by zone (`[calibration] zones`). The example's own set (`[calibration]
  parameters`) is applied on every run, with a calibration's values in place of any it names; `Basin(...,
  parameters=...)`. Water users' parameters, and parameters for a model the example doesn't run, are refused.
- PEST++ IES files (`calibration.write_pest_files`): the control file, a template of the parameter table and the
  instructions for the simulated values, with `python -m basinledger forward` (new) as the model command, each
  worker running in its own folder (`--run-folder`).

## 0.5.1 (2026-11-11)

Needs swatplus-bmi 0.9.0. Its exchange protocol changed in 0.8.0, so a SWAT+ built by an earlier swatplus-bmi has to
be built again (`uv run swatplus-bmi build --jobs 2`). No changes to what Basin Ledger computes.

- swatplus-bmi 0.9.0 required (0.7.0 before); the full suite passes on it unchanged.
- A SWAT+ built by another swatplus-bmi: the protocol error now says to build it again.

## 0.5.0 (2026-11-10)

MODFLOW 6 vertex grids. Needs swatplus-bmi 0.7.0, as 0.4.0 did.

- MODFLOW 6 vertex grids (DISV): maps (`object, layer, icell2d, ...`), wells, head outputs and drawdown rules
  (`[layer, icell2d]`) name a vertex grid's cells, next to a structured grid's `layer, row, col`; existing configs are
  unchanged. `Modflow6.node` reads the grid from MODFLOW 6's discretization and refuses a cell in the other grid's
  form, outside the grid (a column past the edge used to wrap to the next row) or inactive; a config's cells must be
  two or three integers. `add_transport` copies a DISV grid. The synthetic example runs on a refined DISV grid
  (`examples/synthetic/disv.toml`, `model.build_small_disv`), and the Plaster Creek geometry helpers write icell2d
  maps from a vertex grid.

## 0.4.0 (2026-11-09)

Weather from a run's first day, and groundwater up into the soil. Needs swatplus-bmi 0.7.0; a SWAT+ built with an
earlier swatplus-bmi has to be built again (`uv run swatplus-bmi build --jobs 2`).

- Weather into SWAT+'s gage records (`[weather]`, `weather.py`; swatplus-bmi's `set_weather`, and 0.7.0's
  `project.write_weather`): a series per gage from a CSV per variable (precipitation, temperatures, solar radiation,
  humidity, wind), and a scenario's change (`[weather.change]`: factors and temperature shifts, one value or 12
  monthly) on the series or the project's own records. A fresh run's first day (two for precipitation), which SWAT+
  reads before it first pauses, is written into the run's copy of the gage files before it starts, and the rest set
  through the BMI, so the run is the series from its first day; after a restart, the series is set from the first day
  SWAT+ hasn't read, so a case from a warm-up takes its own weather a day later (two for precipitation). A series that
  ends before the run stops it, rather than falling back to the files (decided Nov 9, 2026, as was writing the first
  days). Water users on a SWAT+ station's weather read the series. `from_grid` maps gridded daily values to gages, by
  nearest cell or by the cells under each gage's HRUs.
- Groundwater up into the soil (`[groundwater_to_soil]`, `coupling/gwsoil.py`; swatplus-bmi 0.7.0's
  `soil_layer_water__groundwater_inflow_rate`): where MODFLOW 6's water table rises into an HRU's soil profile, water
  leaves the cell by a drain at the profile's bottom (`rule = "drain"`) or by SWAT+ gwflow's rule (`rule = "gwflow"`),
  optionally no more than the layers' room (`cap = "room"`), and goes into the HRUs' layers below the water table the
  next day. The rule and the cap have no default; the examples use `drain` and `none`, chosen Nov 9, 2026 (the drain
  departs from gwflow's published rule; its conductance is calibrated [A]). Booked as `gw_soil:sent`, `due` and
  `received`; the whole-basin balance and the capture of pumping count it.

## 0.3.0 (2026-11-08)

Recharge HRU by HRU. Needs swatplus-bmi 0.6.1, whose `finalize` waits for a branched run's SWAT+ to exit, so its
summary files, which the whole-basin balance reads, are written before it returns. A SWAT+ built with a
swatplus-bmi before 0.6 has to be built again (`uv run swatplus-bmi build --jobs 2`): 0.6's exchange protocol is 4.

- Recharge HRU by HRU (`[recharge] source = "hru"`): each HRU's percolation out of its soil profile (swatplus-bmi
  0.6.0's `soil_water__percolation_rate`) on the cells under its ground, in place of one rate per aquifer. Each HRU's
  map rows are scaled to its area in SWAT+ (`AreaMap.conserving`), so the volume sent is percolation times HRU area,
  exactly, and the driver holds it to SWAT+'s own aquifer recharge every day (`recharge:swat_aquifers`). The aquifer
  path stays the default.
- `mapping.pixel_rows`: raster HRUs (land cover, soil and slope classes) on a structured grid, each pixel split
  exactly among the cells it overlaps.
- A SWAT+ project with no channels runs (Ames_sub1), and the empty-soil warning reads `soil_plant.ini` headers in
  capitals.

## 0.2.0 (2026-11-07)

Conservation across every link, the checks that hold a coupled run to it, stream temperature and groundwater
transport. Needs swatplus-bmi 0.2 and, for water users, waterledger-models 3.11.

Conservation across links (docs/design.md):
- Recharge keeps the surface model's volume: the area check is on by default (1e-6), and what a map adds or can't
  place is booked as `recharge:by_the_map`.
- Water from outside the models is its own account: each MODFLOW 6 package's flow on entries no link carries
  (`groundwater:boundary:<package>`) and the release from storage (`groundwater:storage_release`).
- The two-way exchange goes through swatplus-bmi 0.2.0's channel groundwater exchange, with river stage from SWAT+'s
  channel depth (`stage = "surface_depth"`), and is limited by supply inside MODFLOW 6's solve: a channel never loses
  more than it carries that day. On by default.
- Drains coupled to channels (`[exchange] drains`): seepage reaches the streams the next day, booked as
  `exchange:<channel>:seepage`.
- The coupling's settings are checked at start-up (SWAT+ on a daily step, two-way aquifers passing recharge through,
  no channel transmission loss, MODFLOW 6's recharge package one entry per cell), and a wrong one stops the run.
- `Basin.spin_up`: a start consistent with the coupling, cycling the first year until the year's release from
  storage is within a tolerance of its recharge.
- Warnings, not stops, for a start that can drift: soils that start empty (`StartWarning`), and a warm-up whose soils
  are still changing by more than 12% of its last year's precipitation (`warm_up.json`).

Checks and analysis:
- `balance.basin_balance`: SWAT+'s land, channel and lake budgets, the links and MODFLOW 6's close together over a
  run; `balance.soil_water_trend`.
- `analysis.capture`: a case's extra pumping as the extra release from storage, the water the streams lose and the
  extra boundary inflow, every day; `analysis.concentration_change` at low flow.
- `analytic`: Glover and Balmer (1954) stream depletion, by superposition, and Ogata and Banks (1961) for a column.
- Scenarios from a shared warm-up (`scenarios.warm_up`, `scenarios.run_case`): a case with no change matches the
  unbroken run bit for bit.
- A study can wrap the surface model (`wrap_surface`) to hold or replace a value through BMI each day.

Stream temperature and water quality:
- Stream temperature as a reach heat balance (`temperature.py`): inflows mixed by volume, groundwater discharge at a
  stated temperature, drains' seepage at the runoff's temperature, and SWAT+'s equilibrium exchange with the air
  integrated exactly. A user's return carries its heat (`route:<point source>:heat_m3C`).
- Groundwater transport (`[transport.<species>]`, `coupling/transport.py`): a MODFLOW 6 GWT model per species, its
  mass in the ledger every day and its closure held to 1e-6; restarts carry the concentrations.
- A user's intake concentration (`water_user_intake_water_<species>__mass_concentration`) is the aquifer's at its
  well, the day before.

Water users:
- Users run on SWAT+'s weather: hourly air temperature from a station's daily minimum and maximum, and its humidity.
- A user's other outputs are booked each day: a power integrated to MWh, an amount of time as the day's increase in a
  running total, anything else as the day's mean.
- `generic.py`: a user from a reported daily series (withdrawal, consumption, return, provenance).
- A user's input names are read after the component is initialized, as BMI intends.

SWAT+ and MODFLOW 6:
- Point sources can be added to a SWAT+ project that has none.
- Print settings are edited with `use_obj_labels` on, since SWAT+ 62.0.1 otherwise reads object rows by position.
- `Modflow6.solve_day` with a callback now solves every solution, not only the first, so a transport model is solved.

NextGen:
- MODFLOW 6 under a NextGen replay: recharge from CFE's percolation per catchment, river stage from t-route's depth
  (`[rivers]`), base flow per flowpath.
- Area maps run in reverse (cells to catchments) and list cells outside the basin.

Calibration: USGS daily values from the Water Data API, recession ratios and the Eckhardt base flow filter.

Examples: Plaster Creek's groundwater model (MODFLOW 6 from public data) and demo, the campus on its own well against
utility water; the Michigan cold-stream screen; the synthetic transport and intake examples.

## 0.1.0 (2026-11-04)

First release.

- A daily driver for coupled models through BMI, with a ledger that enters every transfer as sent and as received
  and stops a run on any day they differ.
- Surface models: SWAT+ live through swatplus-bmi, or a finished NextGen (NGIAB) run replayed day by day.
- MODFLOW 6 through libmf6 and modflowapi: SWAT+ recharge in, wells, and the two-way stream-aquifer exchange through
  the river package, with SWAT+'s aquifer passing recharge through.
- The water-user interface: any BMI component reporting withdrawal, consumption and return (and optionally gain,
  storage, loads, return temperature, electricity), integrated from hourly to daily, its balance checked every step.
  CoolingWaterBmi plugs in by role.
- Routes (wells, SWAT+ intakes and point sources, outfalls, offsite), utility chains, and feedback rules
  (curtailment, dry mode, drawdown limit).
- Restarts from a saved day, matching unbroken runs bit for bit.
- Calibration records per model, statistics, and a forward model; applying parameters and PEST++ are planned.
- Examples: synthetic (Osu_1hru with small MODFLOW 6 models) and Plaster Creek, Michigan (calibrated NextGen
  streamflow).
