# Basin Ledger: design

Basin Ledger is the study engine behind Northline's water studies. It couples watershed, groundwater and water-user
models through BMI and steps them together one day at a time. It is a sibling of
[Water Ledger](https://waterledger.org). Water Ledger accounts for one campus's water: what it withdraws, consumes and
returns, from public data, as a calculator. Basin Ledger accounts for a whole basin's water. It balances quantity and
quality day by day across the coupled models, so a study can say where a large user's water comes from, where it goes,
and what changes in the streams and the aquifer.

Status (Nov 4, 2026). Built and tested on synthetic examples:

- the water-user interface;
- routes through SWAT+ and MODFLOW 6, and the utility chain;
- feedback rules;
- the one-way and two-way SWAT+–MODFLOW 6 couplings;
- restarts;
- the ledger, the command line, and the calibration module's structure.

The first real example is built too: Plaster Creek, Michigan, on a calibrated NextGen model replayed day by day.

Not built yet: the SWAT+ hydrofabric builder, AORC forcing, applying calibration parameters, and PEST++. Each section
says which parts are built.

## Principles

1. **A general engine, with places as examples.** Nothing in `basinledger` is tied to a place. A location is an
   example: one config file that wires real models for one basin and its water users, a script that fetches its public
   data, and a short write-up. An example is also a study's starting point. A public example shows the engine working.
   A client study is a private example built the same way, kept outside this repository. The engine's own tests run on
   a small synthetic example.
2. **Water users are pluggable.** A water user is a BMI component behind one common interface. Data center cooling
   (`CoolingWaterBmi`, from `waterledger-models`) is the first. Mines, battery and EV plants, chip fabs, food processing
   and hydrogen will follow in `waterledger-models`. Each one plugs into the same driver, utility chain and couplings
   with no engine changes. No user type is special-cased.
3. **Every transfer balances, and tests hold it to that.** Each volume or mass that crosses between models is entered
   in the ledger twice, as sent and as received. A day where the two differ stops the run. Each model's own budget
   must close too.
4. **Models stay whole.** Each model runs as its developers ship it. SWAT+ runs through
   [swatplus-bmi](https://github.com/northlinehydrology/swatplus-bmi)
   ([doi:10.5281/zenodo.23963201](https://doi.org/10.5281/zenodo.23963201)), whose output is byte-for-byte identical
   to stock SWAT+. MODFLOW 6 runs through USGS's own shared library. Couplings set fluxes between days; they don't edit
   model code.
5. **Packages first.** Water-user physics lives in `waterledger-models`, SWAT+ access in `swatplus-bmi`. Basin Ledger
   holds only the coupling, the accounting and the rules. Changes those packages need are written as requests
   (`docs/requests.md`).
6. **Calibrated models, or a screening run.** Using Basin Ledger depends on calibrated models. Every example carries
   its calibration and validation record. An example without calibrated models is labelled a screening run, never a
   study result (see Calibration).
7. **NextGen geography throughout.** SWAT+ projects are built on the NextGen hydrofabric, and everything the engine
   reports is keyed to hydrofabric ids, so results compare with the National Water Model (see NextGen geography).

## Layout

```
src/basinledger/
  config.py            an example's TOML config: models, maps, wells, users, routes
  driver.py            Basin: steps every model one day at a time and keeps the ledger
  ledger.py            daily accounts, paired as sent and received, checked every day
  users.py             the water-user interface: roles, units, hourly-to-daily integration, the user's balance
  mapping.py           area maps between SWAT+ objects (catchments, aquifers, HRUs) and MODFLOW 6 cells
  units.py             unit factors to the ledger's SI units
  coupling/
    modflow6.py        MODFLOW 6 through libmf6 (modflowapi): prepare, set, solve, read, restart
    swatplus.py        SWAT+ project edits (places of use, point sources, pass-through aquifers), reads and sets
    recharge.py        SWAT+ aquifer recharge -> MODFLOW 6 recharge
    exchange.py        MODFLOW 6 river exchange <-> SWAT+ channels (two-way)
    nextgen.py         a finished NextGen (NGIAB) run as the surface model, replayed day by day
  forcing.py           hourly weather for users from the surface model's forcing
  routes.py            where a day's water goes: wells, SWAT+ places of use and point sources, outfalls
  chain.py             the utility chain: production, losses, sources, sewer, outfall
  rules.py             feedback rules: curtailment, dry mode, drawdown limits, and study-defined rules
  calibration.py       records, parameter and observation tables, statistics, the forward model (PEST++ stubbed)
  __main__.py          the command line: run, status
  hydrofabric.py       SWAT+ projects and maps from the NextGen hydrofabric, results keyed to its ids (stub)
  testing.py           a synthetic water user for the tests
examples/
  synthetic/           SWAT+'s Osu_1hru + small MODFLOW 6 models + synthetic users (the tests use it)
    basin.toml         one-way: recharge and a pumping user
    surface.toml       a river user, a utility with a well and a plant outfall, three rules
    two_way.toml       the stream-aquifer exchange, over an aquifer of Osu_1hru's own area
    model.py  user.toml  mill.toml  small_user.toml  fetch.py  notebook.py  README.md
  plaster-creek/       Plaster Creek (Kent County, Michigan): calibrated NextGen streamflow, an assumed campus
  <place>/             any example:
    basin.toml         the config
    fetch.py           public data and observations into a cache (outside git)
    README.md          the write-up
    maps/              catchment, aquifer and HRU to cell areas; flowpath to cell lengths (CSV)
    calibration/       record.toml (status, periods, statistics), parameters.csv, observations.csv
tests/                 unit tests, and the one-way coupling run on the synthetic example
docs/
  design.md            this note
  prior-art.md         earlier SWAT/SWAT+ and MODFLOW couplings
  requests.md          changes asked of waterledger-models and swatplus-bmi
```

Each example folder holds:

- `basin.toml`, the config;
- `fetch.py`, which downloads its public data into a cache folder (outside git) and builds or places the model files;
- `README.md`, a short write-up: the basin, the models and their versions, calibration status, what the example shows;
- `maps/`, the maps between the hydrofabric and the MODFLOW 6 grid (see NextGen geography);
- `calibration/`, its calibration and validation record (see Calibration). Without one, the example is a screening
  run;
- optionally a builder module for its MODFLOW 6 model and a notebook.

## Surface models: SWAT+ or NextGen

An example names one surface model: `[swatplus]` or `[nextgen]`.

- **SWAT+, live.** swatplus-bmi pauses SWAT+ between days, so the driver sets withdrawals and returns for each day
  and reads the result the same day. Everything with feedback needs this: surface routes, the two-way exchange, and
  rules whose actions change the flow they read.
- **NextGen, replayed.** ngen runs every catchment over the whole period and t-route routes afterward, so a NextGen
  run can't be paused between days. Basin Ledger reads a finished NGIAB run (`coupling/nextgen.py`): daily means of
  t-route's hourly flow and depth per flowpath, each flowpath's nexus from the hydrofabric subset, and the run's
  forcing for users. The run appears to the driver through the same BMI calls as SWAT+, so outputs, rules, utility
  chains and the ledger work unchanged. Objects are hydrofabric flowpaths, so every result is keyed to NextGen ids.
- **What a replay can't do.**
  - There is no feedback: users in Basin Ledger don't change NextGen's flows. That's exact for a user whose water
    comes from and goes outside the basin. Users that take from or return to a stream must be inside the NextGen run,
    as nextgen-michigan's WaterUsers component is.
  - There are no surface routes and no two-way exchange.
  - MODFLOW 6 needs the run to output CFE's percolation (`SOIL_TO_GW_FLUX`). The Plaster Creek example reruns
    NextGen with it. MODFLOW 6 then runs one-way under the replay: recharge in, and river stages from t-route's
    depth. Its base flow is reported beside NextGen's flows and doesn't change them.

  The driver refuses each of these with a message.
- **Why both.** NextGen is the calibrated streamflow model where one exists, on the frame the National Water Model
  uses. SWAT+ is the model for water quality and for feedback. Where a study needs both, the plan is SWAT+ built on
  the same hydrofabric (NextGen geography), calibrated to the same gages, with NextGen's calibrated flow as a check.

## The models and how each is driven

| Model | Access | Step | What Basin Ledger reads and sets |
|---|---|---|---|
| SWAT+ 62.0.1 | `swatplus-bmi` (BMI 2.0, SWAT+ in its own process, commands over standard input and output) | 1 day | reads channel flow and storage, aquifer storage, depth to water, aquifer outflow and recharge; sets point source flow and water allocation demand for the next day |
| MODFLOW 6 (6.8.2 tested) | `libmf6` through BMI/XMI, from Python with USGS `modflowapi` (in process) | 1 day (one TDIS step) | sets recharge (RCH `RECHARGE`) and well rates (WEL `Q`) between prepare and solve; reads simulated flows (`SIMVALS`), heads (`X`) |
| NextGen (NGIAB v2.1.0 runs) | a finished run's t-route output, replayed (`coupling/nextgen.py`) | 1 day (means of hourly values) | reads flow and depth per flowpath; nothing is set |
| Water users | BMI components from `waterledger-models` (in process) | the component's own (hourly for `CoolingWaterBmi`) | sets weather inputs and limits; reads withdrawal, consumption, return, gain, constituents, temperature, power |

**The day.** For day *d* the driver:

1. runs the feedback rules on yesterday's state (SWAT+ channel flow, MODFLOW 6 heads) and sets each user's controls
   for today;
2. steps each water user through the day (24 hourly steps) and integrates its volumes;
3. routes each user's withdrawal and return, directly or through its utility chain, along with yesterday's
   stream-aquifer exchange;
4. sets today's SWAT+ demands and point sources, runs SWAT+ for day *d* (`update()`), and reads back what it withdrew
   and received. Values read at the pause are day *d*'s;
5. prepares MODFLOW 6's step for day *d*, sets day *d*'s recharge from SWAT+ and the day's well pumping, solves, and
   reads back the wells and the river exchange;
6. enters both sides of every transfer in the ledger, checks them, and records the day's flows and heads.

**Outputs.** `python -m basinledger run basin.toml` (or `Basin.write_outputs()`) writes three files:

- `ledger.csv`: every account, every day;
- `series.csv`: channel flows and heads at named cells;
- `summary.json`: the calibration status, totals, and water still in transit at the end.

**A SWAT+ detail.** SWAT+ 62.0.1 keeps one point of diversion per channel (`sd_ch%wallo_pod`), so several places of
use on one channel share one row in `point_of_diver.wro`, which lists them all. SWAT+ holds its values as 4-byte reals,
so transfers into and out of SWAT+ agree to a relative 1e-6 in the ledger. Transfers between the 8-byte models agree
to 1e-9.

**A MODFLOW 6 detail.** MODFLOW 6.5 and later keep a package's input in named arrays (`RCH/RECHARGE`, `WEL/Q`) and
copy them into `BOUND` when the step is prepared. A value written to `BOUND`, as in the examples in Hughes et al.
(2022), is overwritten before the solve. Basin Ledger sets the named arrays, after `prepare_time_step` and before
`do_time_step`. List packages are filled only when a period is read, so entries are found after the first prepare.
The step length is known only after prepare, so the one-day check happens then. libmf6 holds one simulation per
process. Runs can follow one another in a process but can't overlap, so scenarios run in separate processes.

## Water users

A water user is any BMI component that reports its water through these roles. A component that uses other names is
mapped role by role in the example's config (`variables = {withdrawal = "...", gain = "..."}`). Nothing in the engine
knows what kind of user it is.

| Role | Common variable name | Units | Required |
|---|---|---|---|
| withdrawal | `water_user_water__withdrawal_volume_flux` | m3 s-1 | yes |
| gain | `water_user_water__gain_volume_flux` | m3 s-1 | no: water entering by other ways (condensate from air, ore moisture, water made in a process) |
| consumption | `water_user_water__consumption_volume_flux` | m3 s-1 | yes, if any: evaporated, or leaving in product or waste |
| return | `water_user_water__return_volume_flux` | m3 s-1 | yes, if any: discharged to a sewer or an outfall |
| storage | `water_user_water__volume` | m3 | no: water held on site (ponds, tailings) |
| temperature | `water_user_return_water__temperature` | degC | no |
| constituents | `water_user_return_water_<name>__mass_flux` | kg s-1 | no: one per constituent, in the return |
| power | `water_user__electric_power` | MW | no: facility electricity, for the power side of a trade-off |
| limit (input) | `water_user_water__withdrawal_volume_flux_limit` | m3 s-1 | no: a cap a feedback rule sets; infinity means none |

Units are converted on reading (`units.py`): m3 d-1, L s-1, kW and kg d-1 are accepted, and anything else is
refused. Names follow the [CSDMS Standard Names](https://csdms.colorado.edu/wiki/CSDMS_Standard_Names) pattern.

**The user's balance.** Every step, withdrawal + gain = consumption + return + change in storage. The interface checks
it on each step, to a relative 1e-9, and refuses a component that doesn't close. `CoolingWaterBmi` closes with its
condensate mapped as a gain.

**Hourly to daily.** The hydrologic models step daily, and users usually step hourly. `WaterUser.advance_day(d)` steps
the component from wherever it is to the end of day *d*, setting its inputs before each step (hourly weather from the
shared forcing, limits from the rules). After each step it reads the rates, which describe the step just ended, and
multiplies by the step length:

- volumes (m3) are summed: the daily models get the day's total, or the day's mean rate as total / 86,400 s;
- constituent masses (kg) are summed. The day's concentration in the return is mass / return volume, never a mean of
  hourly concentrations;
- return temperature is flow-weighted: the sum of return volume × temperature, over the return volume (a heat
  balance);
- electricity is summed to MWh.

A component with a daily or other step works the same way. The driver only asks it to reach the end of the day.

**Routing.** Each user's config names where its withdrawal comes from and where its return goes:

- `well:<name>`: pumped from a MODFLOW 6 well;
- `swat:place_of_use:<name>`: a SWAT+ water allocation demand, taken from its channel the same day. What the channel
  couldn't supply is entered as a shortfall. The places of use are added to the run's copy of the project from the
  config (`[swatplus] places_of_use`);
- `swat:point_source:<name>`: a SWAT+ daily point source into its channel (`[swatplus] point_sources`);
- `outfall:<name>`: a return entered in the ledger only;
- `utility:<name>`: the utility chain.

Each route gets `requested` and `delivered` accounts in the ledger, paired every day.

**The utility chain** (`chain.py`, built in its first form) is accounting, not a model:

- **source:** the utility's withdrawal, split among its wells, intakes, and sources outside the basin;
- **distribution losses;**
- **users:** every water user, plus the system's existing demand, as a daily series;
- **sewer:** with inflow and infiltration;
- **wastewater plant;**
- **outfall:** a SWAT+ point source, with loads and temperature.

Each link is a ledger pair. A user on municipal supply is routed `utility:<system>`, and the chain places the water.
A user with its own wells or intake is routed directly. Still to add:

- losses returned to the aquifer as recharge;
- inflow and infiltration taken from the aquifer;
- a daily demand series;
- the plant's own consumption and loads.

**Feedback rules** (`rules.py`) run at the pause before each day, on yesterday's state:

- **curtailment:** cap a user's withdrawal after low flow in a channel;
- **dry mode:** run a user dry after low flow;
- **drawdown limit:** cap a user's withdrawal after a low head in a cell.

A study can add its own rule as a function in its example folder. When rules disagree, the strictest control wins.
The ledger records each rule's decision every day. A rule asking a user for a control the component doesn't take is
refused when the run starts.

## NextGen geography

To be usable with NextGen, Basin Ledger adopts NextGen's geography. The spatial frame for every real example is
NOAA's NextGen hydrofabric
([v2.2](https://lynker-spatial.s3-us-west-2.amazonaws.com/hydrofabric/index.html)). It has three kinds of object:

- catchments (`cat-`);
- flowpaths (`wb-`), one per catchment;
- nexuses (`nex-`), where flowpaths join.

Each example pins the hydrofabric version it was built on.

**SWAT+ on the hydrofabric.** A SWAT+ project for an example is built from the hydrofabric, not delineated anew:

- **Subbasins:** one per catchment. SWAT+'s routing units and aquifers follow, one of each per catchment.
- **Channels:** one per flowpath, connected as the hydrofabric's flowpath-nexus-flowpath network connects them. Channel
  geometry starts from the hydrofabric's flowpath attributes (length, slope, bankfull width and depth).
- **HRUs:** split within each catchment by land cover (NLCD), soils (gSSURGO) and slope, as SWAT+'s tools do. Each HRU
  keeps its catchment's id.
- **Nexus output:** each nexus's flow is the outflow of the channels draining to it, read every day through
  swatplus-bmi and reported under the nexus id.
- **The builder** (`hydrofabric.py`, planned) writes SWAT+'s connection and object files from the hydrofabric
  GeoPackage and a lookup table (`swat_ids.csv`) of SWAT+ object name against hydrofabric id. SWAT+'s usual tools
  (QSWAT+, SWAT+ Editor) delineate their own watersheds, so this builder is new work. Weather and parameters still go
  through SWAT+'s own input files.

**MODFLOW 6 to and from catchments.** MODFLOW 6 keeps its own grid, which usually reaches past the basin to its own
boundaries. Maps tie the grid to the hydrofabric. They are made in GIS (flopy's `GridIntersect` on the hydrofabric
layers), written to `examples/<place>/maps/` as CSV, and named in the example's config. Coupling code never sees
geometry, and nothing is hard-coded.

| Map | Rows | Used for |
|---|---|---|
| `catchment_cells.csv` | catchment (or SWAT+ aquifer, or HRU), layer, row, col, area_m2 | recharge in; drawdown and storage out |
| `flowpath_cells.csv` | flowpath, layer, row, col, length_m | stream-aquifer exchange (base flow) out, and in for losing reaches |
| wells in `basin.toml` | name, cell (or MAW connections) | pumping in; heads at wells out |

- **Recharge in.** Each catchment's recharge from SWAT+ (from its aquifer now; from its HRUs once swatplus-bmi exposes
  HRU percolation) spreads over cells by area of overlap. Volume is kept exactly (built: `mapping.py`).
- **Base flow out.** Each day, MODFLOW 6's stream-aquifer exchange (RIV or SFR `SIMVALS`) is summed per flowpath
  through the length map. That gives groundwater discharge to each flowpath, reported under its `wb-` id. In two-way
  mode this is what enters SWAT+'s channel (see The MODFLOW 6 coupling).
- **Well effects out.**
  - Streamflow depletion per flowpath is the difference in base flow between a run with the user and one without.
  - Drawdown per catchment is the area-weighted mean and the maximum over the catchment's cells. Drawdown at named
    cells (utility wells, sentinel wells) is reported by cell.
  - Each is keyed to the flowpath or catchment id.
- **Back to catchments.** The same area map runs in reverse: a catchment's value is the sum of cell values × overlap
  area, over the overlap area. Cells outside every catchment are reported together as "outside the basin", so the
  grid's water still balances.
- **Checks.** Every map is checked when the run starts:
  - no cell may be covered by more than its area;
  - each catchment's mapped area must match its hydrofabric area within the example's tolerance (`area_rtol`);
  - each flowpath's mapped length must match its hydrofabric length within the same tolerance.

**Keyed to hydrofabric ids.** Everything the engine reports is keyed to hydrofabric ids:

- ledger accounts (`recharge:cat-…`, `baseflow:wb-…`);
- daily output tables, with an id column;
- figures.

Streamflow at a nexus or flowpath then compares directly with the National Water Model's retrospective and
operational output for the same id. Water users are points and are reported by name, with the catchment they sit in.
The utility chain places their water by source.

**The synthetic example** has no hydrofabric. Osu_1hru's objects keep their SWAT+ names, and its map is
`all:aqu001`. The engine treats a hydrofabric id and a SWAT+ name the same way, as an object name in a map, so
nothing in the code depends on which kind it is.

**ngen.** NextGen's own framework runs one model per catchment. SWAT+ and MODFLOW 6 each solve their whole domain at
once, so Basin Ledger uses NextGen's geography, not its runtime. Because results are keyed to nexus ids, a Basin Ledger
result can be handed to ngen at a nexus, as an upstream inflow.

## Forcing: AORC

All models use one weather source. Built now: a user can take its hourly weather from a NextGen run's own forcing
(`forcing.py`, `forcing = {source = "nextgen", catchment = "cat-…"}`). That is temperature, and relative humidity from
specific humidity, temperature and pressure (Bolton 1980), for the user's catchment. So the campus and NextGen see the
same hours. For SWAT+ examples, NOAA's Analysis of Record for Calibration (AORC) is the natural pick: hourly, about
1 km, 1979 onward, the forcing NextGen and the National Water Model are calibrated with, and open on AWS. It is
chosen for this engine on its own merits, independently of Water Ledger's weather (ERA5).

- SWAT+ gets daily precipitation, minimum and maximum temperature, solar radiation, humidity and wind, aggregated from
  AORC over each catchment, in SWAT+'s weather files.
- Water users get hourly air temperature and relative humidity (from AORC's specific humidity, temperature and
  pressure) at their site, through the `inputs` hook.
- MODFLOW 6 needs weather only if an example uses its ET or unsaturated zone packages.

## The MODFLOW 6 coupling

Prior work is summarized in `docs/prior-art.md`. Every published SWAT or SWAT+ coupling with MODFLOW compiles
MODFLOW-NWT into one executable, or uses SWAT+'s own gwflow module. The only SWAT+ to MODFLOW 6 coupling in Python
we found is a 2026 conference abstract (Nikiel and Żurek), with no code. Basin Ledger's coupling differs in four
ways:

- the two models stay separate programs, each unmodified;
- they exchange water through standard interfaces;
- every exchange is checked by a ledger;
- the coupling is general over places and water users.

**One-way: recharge and wells (built).**

- After SWAT+'s day, the driver reads each aquifer's recharge (`soil_water_table__recharge_rate`, m s-1).
- The area map spreads it over cells. A cell's rate is the sum over the objects covering it of rate × overlap area,
  divided by the cell's area. Volume is kept exactly.
- The rate goes into MODFLOW 6's recharge package for the same day.
- Users' well withdrawals are summed per well and set in the well package.
- In the ledger, recharge sent (rate × mapped area) equals recharge received (RCH `SIMVALS`), and well pumping
  requested equals pumped, every day.

In one-way mode SWAT+ keeps its own aquifer and baseflow, and MODFLOW 6 is a reading of what the same recharge does to
heads and to stream leakage. Each model keeps its own calibration. This is enough to show drawdown at a utility's wells
and when depletion reaches a stream. It is not enough to change the streamflow SWAT+ reports.

**Recharge from HRUs.** SWAT+'s aquifer is lumped per subbasin, so mapping aquifers to cells spreads recharge evenly
within each subbasin. Bailey et al. (2025) map HRU recharge to cells by area. That needs HRU percolation through BMI,
which swatplus-bmi doesn't expose yet (`docs/requests.md`). The map format already takes any object type.

**Two-way: stream-aquifer exchange (built with RIV; `coupling/exchange.py`, `examples/synthetic/two_way.toml`).**
MODFLOW 6 computes the exchange between streams and the aquifer (RIV now, SFR for routed streams later), and that
exchange replaces SWAT+'s aquifer return flow.

- **No double counting.** SWAT+'s aquifer must stop returning water once MODFLOW 6 does.
  - This is done in the run's copy of the SWAT+ inputs, not in SWAT+'s code. SWAT+ 62.0.1 returns flow when the depth
    to water is at or above `flo_min` (m) and allows revap above `revap_min`, so both are set to a depth the water
    table can't reach (−1e6 m), and `rchg_dp` is zero. SWAT+'s aquifer then only receives recharge and holds it.
  - The ledger holds SWAT+'s aquifer outflow (`aquifer_water__volume_outflow_rate`) at zero every day.
  - A test checks that SWAT+'s aquifer storage grows by exactly the recharge sent to MODFLOW 6.
  - Recharge is counted once, in MODFLOW 6.
- **Per channel, net.** The exchange on day *d* is summed per SWAT+ channel through the channel map (`all:<channel>`,
  or a CSV of RIV cells by channel). The driver adds one point source (`gw_<channel>`) and one place of use
  (`seep_<channel>`) per coupled channel to the run's copy of the project.
- **Into SWAT+.** Since swatplus-bmi 0.2.0, the day's net exchange is set directly as the channel's groundwater
  exchange (`channel_water__groundwater_volume_exchange_rate`) on day *d*+1. A gain joins the channel's inflow before
  routing. A loss comes out of the day's outflow, no more than is there, and counts as SWAT+'s channel seepage, so
  SWAT+'s channel balance closes. The ledger enters what SWAT+ reports it actually exchanged.
  - Before 0.2.0, the exchange went through one point source and one water allocation place of use per channel.
  - A negative point source would not work: SWAT+ 62.0.1 records it in full even when the channel can't supply it.
- **Limited by supply (built Nov 6, 2026; `[exchange] limit`, on by default).** A channel's net exchange (its
  river cells' gains and losses and its drains' seepage) can't take more than the channel has without today's
  exchange, which won't be there tomorrow when this one is taken. The limit works inside MODFLOW 6's solve: after
  each outer iteration, a channel taking more than that at the heads as they are has its stages lowered by one
  depth, no entry below its bottom (`limit_to_supply`), and the day is done when the solve has converged and no
  stage moves more than 1 mm. The ledger records the depth (`exchange:<channel>:stage_lowered_m`). Two earlier
  forms failed on the Clam: a limit set before the solve, at the day's starting heads (wells and drains beside a
  river pull its heads down during the day), and a limit that counted the channel's own gain from the day before,
  so it lost on alternate days what it had gained.
- **Shortfalls.** What the channel still can't supply is entered in the ledger: water MODFLOW 6 took from a stream
  that, in SWAT+, wasn't there. Most of what remains is the one-day lag. Over 120 days at the Clam the shortfall is
  16.3 million m3 without the limit (35% of the exchange) and 0.2 million m3 with it (0.7%); the run takes about
  twice as long. In the synthetic two-way example it is 186 m3 over three years without the limit and 76 m3 with
  it, of about 61,000 m3. Each day starts from its own stages, so a restarted run matches.
- **Stage.** With `[rivers] stage = "surface_depth"`, MODFLOW 6's river stage on day *d* is the river bottom plus
  the bed thickness plus SWAT+'s channel depth for day *d* (`channel_water__depth`, on SWAT+'s rating curve,
  swatplus-bmi 0.2.0). Otherwise it is fixed in the input. Under a NextGen replay the same setting uses t-route's
  depth.
- **Timing.** The exchange is explicit with a one-day lag. Volume is conserved (ledger pairs): each day's exchange is
  due in SWAT+ the next day, and the last day's is reported as in transit. Timing shifts by a day, which is small
  against the time scales of groundwater response. Iterating within the day by branching SWAT+
  is possible but costs a SWAT+ day per iteration. It will be tried if the lag matters for an example's question.
- **Mapping.** Channels map to cells by length (`object, layer, row, col, length_m`). Conductance comes from the
  MODFLOW 6 model, not from the map.

**Wells.** A utility's wells are MODFLOW 6 wells named in the config, each by cell. A well screened over several layers
will use the MAW package. The chain splits a utility's groundwater withdrawal among its wells by share or by capacity.
A drawdown-limit rule reads the head at a well's cell or at a sentinel cell.

## Balances: quantity and quality

The ledger enters, every day, each transfer as sent and as received, plus each component's own balance terms.

- **Held now, by tests on the synthetic example:**
  - recharge sent = received, against SWAT+'s printed aquifer recharge and MODFLOW 6's budget file;
  - every route: requested = delivered (+ shortfall for SWAT+ withdrawals), and never more delivered than asked;
  - the channel's outflow changes, the same day, by what the users' routes put in and took out;
  - the utility chain: production = delivery + losses = the sum of its sources; sewer = outfall;
  - each user's withdrawal + gain = consumption + return, every hour and every day;
  - two-way: SWAT+'s aquifer returns nothing and stores exactly its recharge; each day's river exchange arrives in
    SWAT+ the next day, whole, with any shortfall reported;
  - MODFLOW 6's budget (storage, recharge, wells, river) closes every day to 1e-7 of the gross flow;
  - restarts: a run saved and carried on matches the unbroken run bit for bit.
- **Quality.** Constituents are masses (kg), never concentrations, across every transfer:
  - a user's return carries its loads (dissolved solids concentrated by evaporation, a cooling tower's blowdown) to the
    outfall;
  - the outfall becomes a SWAT+ point source with those loads;
  - temperature is carried as heat (volume × temperature) and becomes a temperature only at the receiving water.
  - Needs: point source loads settable through swatplus-bmi, and intake concentrations as an input to user components
    (both in `docs/requests.md`).
- **Whole basin.** For each example, a basin-wide check closes over the run:
  - precipitation = evapotranspiration + streamflow out + change in soil, aquifer and channel storage + net user
    consumption + transfers out of the basin;
  - the same for each constituent's mass, from SWAT+'s and MODFLOW 6's own budgets plus the ledger.

## Conservation across links (Andrew, Nov 6, 2026; written before the fixes)

The models must never fight. Water may enter from outside (a supply, recharge where no surface model runs, a fixed
head), but every crossing of a model boundary is an account in the ledger, and no link makes or loses water. Sent =
received is necessary, not sufficient: the Clam's first coupled run (clam-basin-study, stage 3) balanced every pair
and still had four leaks, each of which the pairs hid:

| Leak | Size at the Clam | Why the pairs passed |
|---|---|---|
| Recharge mapped onto more area than the surface model makes it on (lakes in the catchments; a catchment outside SWAT+ given the mean rate) | 28 million m3, 9% of recharge | "sent" was figured on the mapped area, not SWAT+'s |
| Losing reaches taking water the channel didn't have (headwater channels left nearly dry once SWAT+'s aquifers pass recharge through) | 84 million m3, 30% of recharge | booked as a shortfall, while MODFLOW 6 kept the water |
| A starting state not consistent with the coupling (a steady state solved with every reach able to lose) | 49 million m3 of storage drained | storage is inside MODFLOW 6's budget, not a link |
| River cells outside the coupled channels at a fixed stage | an unlimited source or sink | not mapped, so not in the ledger |

Rules, each to be built and tested:

1. **One owner per process.** Each flux is computed by one model. The other models' copy of it is switched off or
   masked: set through BMI each day to the owner's value or to zero (swatplus-bmi's settable inputs, as the soil
   temperature in clam-basin-study). SWAT+'s aquifer return flow, revap and deep recharge are MODFLOW 6's once it is
   coupled; channel transmission loss is the river exchange's. Masking through BMI is preferred to editing inputs,
   because a later input (`calibration.cal`) can undo an edit.
2. **Recharge conserves the surface model's volume (built Nov 6, 2026).** A recharge map's area per object equals the
   object's area in the source model: the area check is on by default (`area_rtol`, 1e-6), and an example that maps
   more or less on purpose turns it off and says why. What the surface model made is booked beside what was sent;
   the difference (`recharge:by_the_map`) is water the map added or couldn't place.
3. **Water from outside is its own account (built Nov 6, 2026).** Each day, every MODFLOW 6 package's flow on entries
   no link carries (fixed or general heads, drains and rivers not coupled, recharge packages Basin Ledger doesn't
   set, wells not in the config) is booked as `groundwater:boundary:<package>`, positive into the aquifer. The
   release from storage (STRGSS + STRGSY) is `groundwater:storage_release`, and the model's closure against the
   day's gross flow (`groundwater:closure_share`) must be within 1e-6 every day.
4. **An exchange is limited by what the giver has (built Nov 6, 2026).** A losing reach can't take more than its
   channel carries that day. The day's river stage is lowered, smoothly, until the predicted loss at the day's
   starting heads fits the channel's flow. What still can't be supplied after the solve stays a shortfall, and should
   be small (see The MODFLOW 6 coupling).
5. **Seepage reaches the streams (built Nov 6, 2026).** Drains (seepage at land surface) are coupled to channels
   like the rivers' gains (`[exchange] drains`): their discharge joins the channel's exchange the next day and is
   booked as `exchange:<channel>:seepage`. Without it a drain is a boundary flow, water leaving the basin.
6. **The start is consistent with the coupling (built Nov 6, 2026).** `Basin.spin_up` runs the first year again and
   again, SWAT+ from its own start and MODFLOW 6 from the heads the last cycle ended with, until the year's net
   release from storage is within a tolerance of its recharge; `Basin(..., heads=...)` starts the run there. The
   ledger's `groundwater:storage_release` shows any drift that remains.
7. **Each model is checked for the settings coupling needs, at start-up (built Nov 6, 2026).** Below; a run with a
   wrong setting stops before its first day.
8. **The whole basin closes** (see Balances), from each model's own budget plus the ledger.

Settings the coupling needs, checked when a run starts:

- **SWAT+, two-way:** aquifers pass recharge through (no return flow, revap or deep recharge), and nothing in
  `calibration.cal` turns them back on (`flo_min`, `revap_min`, `rchg_dp`, `chk`); channel transmission loss (`k`
  in `hyd-sed-lte.cha`) is zero. SWAT+ 62.0.1 (`ch_rtmusk.f90`) takes that loss out of the channel, records it as
  `seep` and sends it nowhere, so in two-way mode it would leave the basin unbooked; Basin Ledger sets it to zero
  in the run's copy, as it does the aquifers (`no_transmission_loss`). Osu_1hru's channels had it on, so the
  synthetic two-way example lost water this way until Nov 6, 2026. Daily time step (`time.sim` step 0).
- **MODFLOW 6:** a one-day time step; a recharge package Basin Ledger sets has one entry per cell and no auxiliary
  multiplier (an entry's rate is its cell's rate); river and drain entries a map couples are one per cell; packages
  not coupled are reported as boundary flows.
- **Water users (Water Ledger's models):** to be listed with the first coupled user in a real example.

## Water quality across links (Andrew, Nov 6, 2026; written before the build)

The coupling's real benefit is that chemistry follows the water through every link, as volume does now. A campus's
intake is the aquifer's own water at its well, and its return goes into the stream's own water quality, not into
an analysis after the run. Andrew: "we do need to do the work to get the real benefit of the coupling."

**Constituents, to start:**
- chloride: conservative, and road salt is Michigan's main source;
- dissolved solids: the sum of SWAT+'s eight salt ions, conservative as a total;
- phosphorus: sorbs in the aquifer, by a linear isotherm [A];
- nitrate: later, since SWAT+ carries it by default.

**Each link, with its mass in the ledger as sent and received:**

| Link | From | To | Needs |
|---|---|---|---|
| recharge's mass | SWAT+'s leached mass per aquifer (salt module `perc_*`; soluble P leached) | MODFLOW 6 GWT: the recharge package's concentration (SSM) | swatplus-bmi: recharge mass fluxes readable |
| the aquifer to the streams | GWT: discharge times the cell's concentration (river and drain cells) | SWAT+'s channel, with the gain | swatplus-bmi: a gain's mass settable |
| the streams to the aquifer | SWAT+'s channel concentration | GWT: the river package's concentration | swatplus-bmi: channel concentrations readable |
| a well | GWT: removed at the cell's concentration | the user's intake concentration, set each day | built in Basin Ledger |
| a return | the user's returned mass | SWAT+'s point source loads (salt recall loads; soluble P) | swatplus-bmi: point-source loads settable (listed already) |

**MODFLOW 6:**
- a groundwater transport model (GWT) on the flow model's grid, linked by a GWF-GWT exchange;
- packages: advection; dispersion [A]; mobile storage with porosity [A] and, for phosphorus, sorption [A]; sources
  and sinks (SSM) reading each boundary's concentration from an auxiliary variable Basin Ledger sets each day;
- Basin Ledger reads the concentration (GWT's X) at the wells and the river and drain cells.

**Starting concentrations:** a steady-state transport solution under the run's mean recharge concentrations. It is
checked against public well samples (Water Quality Portal) and the streams' late-summer chemistry, which is mostly
groundwater. Transport settles over decades, not the year a flow warm-up takes, so the steady state is the start.

**Stages, each with its tests before the next:**
1. GWT in Basin Ledger with stated recharge concentrations [A]: the mass ledger; the well's intake from GWT; tests
   that mass closes every day, that a steady state with uniform recharge concentration gives that concentration
   everywhere, and a one-dimensional transport case against Ogata and Banks (1961).
2. The streams' side through swatplus-bmi (recharge mass readable, a gain's mass and point-source loads settable,
   channel concentrations readable): requested of the swatplus-bmi session once stage 1 shows which variables are
   needed. Until then the return's concentration change is worked out after the run (`analysis.concentration_change`).

**Stage 1, built (Nov 6, 2026).** `coupling/transport.py`, `[transport.<species>]` in an example's config, and
`examples/synthetic/transport.toml`.
- **Building:** `add_transport(sim, gwf, species, porosity, dispersivity, sorption, strt, sources)` in an example's
  builder adds `gwt_<species>` with its own IMS, the GWF6-GWT6 exchange and MST, ADV (TVD by default), DSP
  (transverse dispersivities a tenth and a hundredth of the longitudinal [A]), IC, SSM and OC. Each source package
  is given its auxiliary variable if it lacks it.
- **Each day:** after MODFLOW 6 prepares the day, the recharge's concentration and that of the coupled rivers'
  losses are set, both stated [A]. They are set on the transport model's copy of the flow package's auxiliary values
  (`GWT_<NAME>/FMI-FT<n>/AUXVAR`): set on the flow package itself after prepare they do nothing, and before it they
  act a day late. The solve with the supply limit's callback now solves the transport solutions after the flow's
  (it had solved only the flow's).
- **Intake:** after the solve, the concentration at each configured well becomes the next day's
  `water_user_intake_water_<species>__mass_concentration` of each user withdrawing from it. The users run before
  MODFLOW 6 each day, so the intake lags the aquifer by a day, as the exchange does; the mass that lag moves is the
  difference between `transport:<sp>:well:<w>:withdrawn` and the users' `intake:<sp>:kg` (0.03% over two years on
  the synthetic run). A user on a well through a utility chain isn't given the intake yet.
- **Ledger (kg):** recharge as sent (stated concentration × the volume sent) and received (the transport model's
  intake of it); each well's mass withdrawn and its concentration; each coupled channel's mass to the stream (gains
  and drains) and from it (losses, as sent and received); boundaries not coupled, by package, and GWT's own (CNC,
  SRC); release from mobile and sorbed storage (MST's rates); and `transport:<sp>:closure_share`, held to 1e-6 of
  the day's gross mass flow every day (it closes to about 1e-12). Inflow carries the auxiliary value, outflow the
  cell's concentration at the end of the day, as MODFLOW 6's implicit step has it.
- **Refused at start-up:** a transport model with no `[transport]` section, a recharge or coupled river package
  that isn't an AUX source in SSM, a coupled exchange with no stated stream concentration, a constant intake the
  aquifer would replace, SSM's AUXMIXED and file sources, and decay (not yet in the ledger).
- **Restarts** save and restore each species' concentrations with the heads, bit for bit.
- **Not yet:** the steady-state start (the tests start from a stated concentration), concentrations at the river and
  drain cells reported per channel for stage 2, and the stage 2 links below.

**What stage 2 needs of swatplus-bmi** (SWAT+ 62.0.1's own variables; names to be settled there):
- recharge mass per aquifer, readable: each salt ion's `obcs(icmd)%hin(1)%salt(m)` (kg a day, kept as
  `asaltb_d(iaq)%salt(m)%rchrg`; chloride is the sixth of the eight ions), soluble P `ob(icmd)%hin%solp`, nitrate
  `ob(icmd)%hin%no3`;
- a gain's mass, settable with the gain's water: the channel's `obcs(icmd)%trans%salt(m)`, which
  `sd_channel_control3` adds to the inflow's constituents as it adds `ob(icmd)%trans` to its water, and
  `ob(icmd)%trans%solp`;
- a loss's mass taken with the loss: `exch_channel` takes the water but not the salts (`hcs2%salt`), and the mass
  taken should be readable;
- channel concentrations, readable: `ch_water(ich)%salt(m) / ch_stor(ich)%flo` (or `chsalt_d(ich)%salt(m)%conc`, mg/L
  at the end of the day) and `ch_stor(ich)%solp / ch_stor(ich)%flo`;
- point-source loads, settable: the recall's `rec_salt(irec)%hd_salt(day, yr)%salt(m)` (copied to
  `obcs(icmd)%hd(1)%salt(m)` by `recall_salt`) and `recall(irec)%hd(day, yr)%solp`;
- the aquifer's salt passed through as its water is: with return flow off, `asaltb_d(iaq)%salt(m)%saltgw` stays zero
  but `cs_aqu(iaq)%salt(m)` keeps the recharge's mass, a store the ledger would have to hold. SWAT+'s aquifer P
  outflow is the recharge volume times a fixed concentration (`aqu_d%minp = hin%flo * aqu_dat%minp`), not mass
  carried from recharge.
3. In each study, SWAT+'s salt module set up (rain chemistry, fertilizer, road salt), and starting concentrations
   from public samples.

**Stage 2 in swatplus-bmi, and what SWAT+ can't carry (Nov 6, 2026).** The SWAT+ BMI session built swatplus-bmi
0.5.0 on its branch `chem-050`, not yet merged or released. For each species (the eight salt ions, nitrate and
soluble P) it has:
- `channel_water_<sp>__mass_concentration`, the day's outflow, kg m-3;
- `channel_water_<sp>__groundwater_mass_exchange_rate`, kg s-1, set after the day's gain and read back as the mass
  in or out. A loss now takes the outflow's nitrate, P and salt at its concentration, counted as seepage;
- `aquifer_water_<sp>__recharge_mass_rate`, kg s-1;
- `point_source_water_<sp>__mass_flow_rate`, kg s-1, daily recalls.

But stock SWAT+ 62.0.1's salt module can't run in a project with aquifers or point sources: it stops at start-up or
on day 1 (soil arrays left unallocated without a pesticide, a double allocation for aquifers, `salt_recall.rec`
never read). So the salt ions, chloride and dissolved solids, can't come from SWAT+ in a coupled project. Basin
Ledger carries them in the streams itself, as it carries stream temperature:
- the surface's load, from a stated or measured concentration on runoff, including road salt [A];
- the aquifer's discharge at GWT's concentration;
- users' returns;
- mixed and routed reach by reach, the mass closed in the ledger.

Nitrate and soluble P can go through SWAT+ once 0.5.0 is released. Two more findings from that session:
- SWAT+'s aquifer P outflow is the day's recharge volume times `minp`, not the return flow, and recharge P isn't
  kept (its design.md, item 11).
- In 0.4.0, a groundwater gain in a project with any constituent on (a pesticide included) hangs SWAT+, because the
  transfer's constituent arrays are never allocated; fixed in 0.5.0.

## Stream temperature: a reach heat balance (built Nov 6, 2026)

SWAT+ 62.0.1's `ch_temp.f90` can't give a stream temperature that responds to groundwater: no aquifer routine fills
the groundwater term it reads (`hdsep%flo_gwsw`), so aquifer return flow and a gain set through swatplus-bmi both mix
in at the transfer hydrograph's temperature, and it mixes HRU water yields in mm with flows in m3/s. Andrew's choice
(Nov 6): a small Basin Ledger component, `temperature.py`, run after a coupled run.

- **Each reach, each day, upstream first:** the water entering is mixed by volume, each source at its temperature
  [A]: upstream reaches at theirs; lakes at `5 + 0.75 x` air temperature (SWAT+'s value for an inflow that isn't a
  channel); runoff, lateral flow and snowmelt at `ch_temp`'s lagged air temperatures; groundwater from MODFLOW 6
  (the ledger's applied exchange, where positive) at a stated temperature, by default the station's mean annual air
  temperature over the run.
- **Heat exchange with the air** over the reach's travel time, by `ch_temp`'s equilibrium-temperature method
  (Edinger and others 1974) and coefficients [A], with a Manning velocity at the Fuller and Peck peak rate [A].
- **Energy closes:** `reach_heat_balance` returns the heat the sources bring, the heat exchanged (surface flux x
  the parcel's surface x travel time) and the heat leaving, plus what the 0.001 C floor adds (ice isn't modelled);
  `tests/test_temperature.py` holds them equal.
- **Dew point:** `tave - (100 - RH%) / 5` (Lawrence 2005). SWAT+'s `tave - (1 - rhum) / 5` takes RH as a fraction;
  `ch_temp`'s coefficients were set with it, which is one reason the coefficients need checking against measured
  temperature in each example.
- An example asks for it with `[outputs] temperature = true`: SWAT+ then prints `lsunit_wb` daily in the run's copy.

## Calibration

Using Basin Ledger depends on calibrated models. Calibration is part of every example, not something a study adds
later. `calibration.py` has the structure:

- the record and its status;
- parameter and observation tables;
- simulated values at observations;
- the statistics;
- a forward model that runs the driver.

Applying parameters to the models' files, PEST++ files and USGS downloads are stubbed. This section is the whole
design.

**Parameters.** The engine exposes each model's adjustable parameters as one named list, with bounds:

- **SWAT+:** parameters changed in its input files before a run, through `calibration.cal` (absolute or relative
  changes by object). swatplus-bmi refuses SWAT+'s soft and hard calibration runs, but a run with changed parameters
  is an ordinary run. Whether every `calibration.cal` change applies in exchange mode is to be checked.
- **MODFLOW 6:** package arrays (hydraulic conductivity, specific yield, streambed conductance, recharge multipliers),
  written with flopy, by zone or by pilot points through pyEMU.
- **Water users:** parameters in each component's config file, e.g. `CoolingWaterBmi`'s cycles of concentration.
  These are usually set from the applicant's data, not calibrated.

**Observations.** An example's `fetch.py` loads them into its cache, with their source and retrieval date:

- **USGS streamflow:** daily values from the National Water Information System. Each gage is matched to its
  hydrofabric flowpath or nexus through the hydrofabric's gage crosswalk.
- **USGS groundwater levels:** each well is matched to a MODFLOW 6 cell, with the layer taken from its screened depth.
  State records (in Michigan, Wellogic static water levels) are added where USGS has few wells, weighted for their
  lower quality.
- **Where a question needs them:** water quality from the Water Quality Portal, and base flow separated from gage
  records.

**How models are calibrated.**

- **One-way coupling:** each model is calibrated on its own:
  - SWAT+ against streamflow, with low flows weighted, since low flows are where a large user shows;
  - MODFLOW 6 against heads and base flow, with SWAT+'s calibrated recharge as its input.
- **Feedback** (two-way exchange, curtailment, drawdown limits): the models are calibrated jointly.
  - The individual calibrations give the prior.
  - [PEST++ IES](https://github.com/usgs/pestpp) runs Basin Ledger's driver as the forward model,
    `calibration.run_forward_model(basin.toml, parameters.csv, observations.csv, simulated.csv)`. It applies a
    parameter set through the driver's `prepare` hook, on the run's copy of each model's files before it starts. It
    then runs every model and writes the simulated value at every observation for PEST++ to read. Applying
    parameters is the part still to build.
  - The adjustable parameters are those the feedback touches: streambed conductance, aquifer hydraulic conductivity
    and specific yield near streams, and SWAT+'s percolation and channel parameters.
  - IES gives an ensemble, so results can carry its spread.
- **Shared machines:** calibration runs are parallel, and are held to the cores that are free.

**The order, for an example with SWAT+ and MODFLOW 6** (Andrew, November 5, 2026; Plaster Creek first):

1. **SWAT+ alone,** against streamflow: land, soil and channel parameters. SWAT+'s own aquifer parameters (alpha_bf,
   delay) are a stand-in that the coupled run replaces.
2. **MODFLOW 6 alone,** one-way under SWAT+'s recharge, against heads and separated base flow: streambed
   conductance, specific yield and storage, and K multipliers.
3. **Linked both ways,** with MODFLOW 6 supplying base flow. A joint PEST++ IES pass on the parameters the feedback
   touches, against low flows, heads and base flow together.

Recharge comes from SWAT+, because SWAT+ carries water quality. Where a calibrated NextGen run exists, it is a
comparison, a check on the National Water Model's frame, not the source of recharge. Until an example's SWAT+
project exists, stage 2 may be prepared against NextGen's recharge, flagged as an interim forcing.

**The record.** Each example keeps `calibration/` beside its config, named in `basin.toml` as
`[calibration] record = "calibration/record.toml"`. A record gives a status per model (`[models.nextgen]`,
`[models.modflow6]`, ...) with what it was calibrated to (`target`). An example is calibrated only if every hydrologic
model it runs is. Water users' settings come from assumptions or applicant data and are reported as such, never as
calibrated. The record holds:

- `record.toml`:
  - status (`calibrated` or `screening`) and the date;
  - the calibration and validation periods;
  - model versions and the hydrofabric version;
  - what was calibrated (separately or jointly) and the objective;
  - statistics by site, for calibration and validation. For streamflow: NSE, KGE, percent bias and a low-flow
    statistic. For heads: mean error, RMSE and RMSE over the range of observed heads. Each is judged against
    published criteria (Moriasi et al., 2015).
- `parameters.csv`: the final values, or the ensemble.
- `observations.csv`: every site used (id, source, period, weight) and the hydrofabric id or cell it was matched to.

**Screening runs.** An example with no record, or with a record whose status isn't `calibrated`, is a screening run.
The driver reports its status (`Basin.status`), and every table and figure it writes carries it. A screening run can
show how the engine and the scenarios work. It is never a study result. The synthetic example is a screening run.

## Scenarios from a shared warm-up

A study runs a baseline and several scenarios (no user; the user as proposed; with curtailment; dry cooling; another
source). All of them start from one warm-up:

- **SWAT+:** `branch()` copies the paused SWAT+ process with fork (Mac, Linux), and `save_state()` writes it to a
  folder (any platform). Both are exact (swatplus-bmi's tests).
- **MODFLOW 6:** its state is its heads, plus the stages and storages of advanced packages (SFR, MAW, UZF) if used.
  At the end of the warm-up the driver saves heads (`X`). Each scenario starts a new MODFLOW 6 simulation covering
  only the days left: the example's builder is called for them, or a copied model's single stress period is
  shortened. The saved heads are then put straight into memory after initialize, not into the IC file, which would
  round them. Restarted runs match the unbroken run bit for bit (built and tested). Models with several stress periods
  need their period data re-keyed, which isn't built yet. Advanced packages will need their own state saved.
- **Water users:** a component's state, if any, is saved through `get_state`/`set_state` where the component offers
  them (CSDMS bmi PR #112). `CoolingWaterBmi` keeps no state between hours.
- **Driver:** the ledger, the series, the rules' memory and the exchange due are saved with the rest
  (`Basin.save`, `Basin(..., restart=...)`, `--save-at`/`--restart` on the command line). The users are stepped
  forward to the restart day, which reproduces a deterministic component exactly.

Scenarios then run in parallel, each in its own process. `branch()` (fork) isn't used by the driver yet, because
libmf6 lives in the driver's process. A branch would need MODFLOW 6 restarted from heads in the child process.

## Examples and private studies

- **Synthetic** (`examples/synthetic`, built): SWAT+'s bundled Osu_1hru project (10 ha, 2010–2023), with three
  configs. In `basin.toml` its recharge rate is spread over a 41 × 41 MODFLOW 6 grid of 50 m cells with a river along
  one edge, and a synthetic user pumps one well. `surface.toml` adds a river user, a utility and rules. `two_way.toml`
  couples Osu_1hru both ways with a 10 ha aquifer. It exists for the tests and the notebook. Its places and numbers
  mean nothing.
- **Plaster Creek** (`examples/plaster-creek`, built): Plaster Creek above USGS 04119055 in Kent County, Michigan.
  - **Built:**
    - nextgen-michigan's calibrated NextGen run (16 catchments, hydrofabric v2.2, two-stage calibration), replayed;
    - an assumed 300 MW campus on utility water from outside the basin, stepped by `CoolingWaterBmi` on the run's own
      forcing.
  - **Calibration status:** streamflow at the gage is calibrated, and the record carries both nextgen-michigan's
    published fit and Basin Ledger's recomputation of it. A test checks the two agree with the record. Low flow is
    approximate (7-day lows 8–18% below observed in 2020–2022), and the calibration is for the gage, not the ungaged
    headwaters.
  - **Groundwater (built, screening):** a MODFLOW 6 model under the calibrated surface, from public data
    (`examples/plaster-creek/groundwater`).
    - 250 m cells, 3 layers after SIR 2010-5109's scheme, K by mapped glacial unit and bedrock subcrop.
    - Recharge is CFE's daily percolation from a NextGen rerun.
    - River stages follow t-route's depth (`[rivers]`), and base flow is reported per flowpath.
    - The fit: heads have an RMSE of 7.9 m against Wellogic; mean base flow is within 4% of the gage's separation, but
      too steady.
  - **The demo** (`docs/plaster-creek-demo.md`): the campus on its own well against utility water, by year type;
    Glover against MODFLOW 6; one-at-a-time ranges; the utility's capacity.
  - **Next:**
    - fit the groundwater with PEST++ IES;
    - extend the forcing to include 2012;
    - replay nextgen-michigan's scenario runs;
    - SWAT+ on the same catchments for water quality.
  - **The groundwater model's starting point:** USGS's regional Lake Michigan Basin model (Feinstein et al., 2010) is
    SEAWAT-2000/MODFLOW-2000, 20 layers with about 1.5 km inner cells. The MODFLOW-NWT inset models derived from it
    are archived (doi:10.5066/F76D5R5V). Either can be converted to MODFLOW 6 with flopy and refined around the basin.
    The model's own files for 2010 weren't found in a data release; this needs checking with USGS.
  - **Calibration data for what's next:** Michigan well records for groundwater, and the Plaster Creek TMDL record
    (E. coli, sediment and biota) for water quality. Until those models are calibrated, they are screening.
- **Private studies** live in the private `northline` repository (`studies/<name>/`), never here. A study:
  - depends on `basinledger` at a pinned version;
  - follows the example layout (`basin.toml`, `fetch.py`, `maps/`, `calibration/`, a write-up), with client data in
    its own folders under the client's confidentiality terms;
  - runs the same engine, so a result can be reproduced from its config and its data.

  Anything a study needs from the engine is added here first, generally, with a test on the synthetic example.
  Nothing from a client study comes back into this repository or into Water Ledger.

## What's next

Done on the synthetic examples (Nov 4, 2026):

- the command line and outputs;
- the utility chain;
- the rules;
- surface routes;
- the two-way RIV exchange;
- MODFLOW 6 restarts;
- the calibration module's structure.

Next:

1. `CoolingWaterBmi` limit and dry-mode inputs (asked of `waterledger-models`), then the rules on a real cooling user.
2. Loads on surface routes: point source constituents and temperature (asked of swatplus-bmi), and intake
   concentrations on users.
3. Two-way refinements:
   - stage from SWAT+'s channel depth (asked of swatplus-bmi);
   - SFR for routed streams;
   - leaks and infiltration in the chain exchanged with MODFLOW 6.
4. Restarts for MODFLOW 6 models with several stress periods; scenarios in parallel processes from one saved
   warm-up.
5. Applying calibration parameters (SWAT+ `calibration.cal`, MODFLOW 6 arrays with pyEMU), PEST++ IES files, USGS
   observations.
6. AORC forcing for SWAT+ and the users; the SWAT+ project builder on the hydrofabric, with outputs at every nexus.
7. Plaster Creek, in the calibration order above:
   - the SWAT+ project on the hydrofabric, which blocks stage 1;
   - stage 2 prepared on NextGen's recharge as an interim forcing;
   - extend NextGen's forcing to include the 2012 drought;
   - add existing pumping once it can be placed;
   - SWAT+ on the same catchments for water quality.

## References

- Bailey, R.T., Abbas, S., Arnold, J.G., White, M.J., 2025. SWAT+MODFLOW: a new hydrologic model for simulating
  surface–subsurface flow in managed watersheds. Geosci. Model Dev. 18, 5681–5697.
  https://doi.org/10.5194/gmd-18-5681-2025
- Feinstein, D.T., Hunt, R.J., Reeves, H.W., 2010. Regional groundwater-flow model of the Lake Michigan Basin in
  support of Great Lakes Basin water availability and use studies. USGS SIR 2010-5109.
  https://pubs.usgs.gov/sir/2010/5109/
- Hughes, J.D., Russcher, M.J., Langevin, C.D., Morway, E.D., McDonald, R.R., 2022. The MODFLOW Application Programming
  Interface for simulation control and software interoperability. Environ. Model. Softw. 148, 105257.
  https://doi.org/10.1016/j.envsoft.2021.105257
- Bolton, D., 1980. The computation of equivalent potential temperature. Mon. Weather Rev. 108, 1046–1053.
  https://doi.org/10.1175/1520-0493(1980)108<1046:TCOEPT>2.0.CO;2
- Moriasi, D.N., Gitau, M.W., Pai, N., Daggupati, P., 2015. Hydrologic and water quality models: performance measures
  and evaluation criteria. Trans. ASABE 58(6), 1763–1785. https://doi.org/10.13031/trans.58.10715
- Sommerlot, A., 2026. swatplus-bmi: run SWAT+ one day at a time through BMI. Zenodo.
  https://doi.org/10.5281/zenodo.23963201
- Other couplings: `docs/prior-art.md`.
