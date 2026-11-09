# Changelog

## Unreleased

- Weather into SWAT+'s gage records (`[weather]`, `weather.py`; swatplus-bmi 0.6's `set_weather`): a series per gage
  from a CSV per variable (precipitation, temperatures, solar radiation, humidity, wind), and a scenario's change
  (`[weather.change]`: factors and temperature shifts, one value or 12 monthly) on the series or the project's own
  records. Written from the first day SWAT+ hasn't read, after a start or a restart, so a case from a warm-up takes
  its own weather; a series that ends before the run stops it. Water users on a SWAT+ station's weather read the
  series. `from_grid` maps gridded daily values to gages, by nearest cell or by the cells under each gage's HRUs.

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
