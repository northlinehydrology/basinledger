# Basin Ledger: design

Basin Ledger is the study engine behind Northline's water studies. It couples watershed, groundwater and water-user
models through BMI and steps them together one day at a time. It is a sibling of
[Water Ledger](https://waterledger.org). Water Ledger accounts for one campus's water: what it withdraws, consumes and
returns, from public data, as a calculator. Basin Ledger accounts for a whole basin's water. It balances quantity and
quality day by day across the coupled models, so a study can say where a large user's water comes from, where it goes,
and what changes in the streams and the aquifer.

Status (Nov 26, 2026, release 0.18.0). Built and tested on synthetic examples:

- the water-user interface, with a user's other outputs (power, time, heat not removed) booked each day;
- routes through SWAT+ and MODFLOW 6, and the utility chain;
- feedback rules;
- the one-way and two-way SWAT+–MODFLOW 6 couplings, conserving water across every link (Conservation across links):
  recharge on each aquifer's own area, water from outside as its own account, the exchange limited by what a channel
  carries, drains' seepage to the streams, and the coupling's settings checked at start-up;
- recharge HRU by HRU, each HRU's percolation on the cells under its ground and held to SWAT+'s own aquifer recharge;
- groundwater up into the soil where the water table rises into it, by a drain at the profile's bottom or SWAT+
  gwflow's rule (swatplus-bmi 0.7.0; the drain chosen);
- a start consistent with the coupling (`Basin.spin_up`), and warnings for a start that can drift;
- the soil cover SWAT+'s soil evaporation sees, on the HRUs a study names: a tree's leaves and residue, not its wood
  (Soil cover for soil evaporation; swatplus-bmi's soil cover);
- the whole-basin balance, the capture of pumping, and depletion against Glover and Balmer (1954);
- intakes on any channel or reservoir, several to one served senior first, a user returning only the share its intake
  delivered, the next day, and SWAT+'s own irrigation from its channels and reservoirs booked by source
  (swatplus-bmi 0.14.0; Irrigation from surface water), and users' returns to the ground, into an aquifer cell;
- water systems with customers (towns as systems): customers by class, sources by share or priority with reuse,
  shortfall and next-day returns (Water systems with customers);
- scenarios from a shared warm-up, and restarts;
- stream temperature as a reach heat balance, with users' returns carrying their heat;
- water quality: groundwater transport (MODFLOW 6 GWT) with each species' mass in the ledger, and users' intake
  concentrations from the aquifer; nitrate and soluble phosphorus across every SWAT+-MODFLOW 6 link with decay and
  sorption (Water quality across links, stage 2; a user's returned loads still to come); dispersion with or without
  XT3D's cross terms;
- salt through land, aquifer and streams, and its replay from a recorded flow run (Replay), for calibration;
- weather written into SWAT+'s gage records from a run's first day, from a series or a scenario's change
  (`[weather]`, swatplus-bmi 0.7.0);
- MODFLOW 6 on a structured grid (DIS) or a vertex grid (DISV, such as a quadtree refined around a site), with maps,
  wells, head outputs and rules naming cells in the grid's own form;
- an existing MODFLOW 6 model taken as given (Extending an existing MODFLOW 6 model), also tested on two published
  USGS models (Pleasant Lake, Wisconsin; the Kalamazoo River basin, Michigan);
- lakes: at a fixed stage through their outlet channels, or a MODFLOW 6 LAK lake on its SWAT+ reservoir
  (swatplus-bmi 0.11.0); and irrigation pumped from MODFLOW 6's wells (swatplus-bmi 0.11.0);
- the ledger, the command line (with `basinledger check`, every problem the start-up finds), and calibration:
  records, statistics, the forward model, parameters written into SWAT+'s `calibration.cal` and MODFLOW 6's
  packages, and PEST++ IES files (Calibration).

The first real example is built too: Plaster Creek, Michigan, on a calibrated NextGen model replayed day by day. A
study's coupled SWAT+ and MODFLOW 6 model, kept in its own repository, runs every function (testing.md, tier 5).

Not built yet: water quality's second stage (through SWAT+'s streams), the SWAT+ hydrofabric builder and AORC
forcing. Each section says which parts are built.

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
8. **The links, not the models (Andrew, Nov 11, 2026).** Basin Ledger builds and checks the links between models
   and never builds the models themselves. A link's part is in: its map (which HRU's water reaches which cell, which
   reach belongs to which channel), the settings the coupling depends on, a helper that builds what the link needs
   (`mapping.py`'s HRU map, `coupling/sfr.py`'s reaches), and its checks at start-up (`basinledger check`). A model's
   own content is out: grids, hydrostratigraphy, DEM sampling, streambed elevations, bed properties, SWAT+ projects
   and calibration choices belong to the study or to tools made for them (QSWAT+, flopy, SFRmaker). Where a helper
   needs a model's numbers, they come in from the study, and a helper refuses what it can't use rather than fixing
   it. Each new link (a reservoir, a lake, the grid) gets its map, its helper and its checks, and no more.

## Layout

```
src/basinledger/
  config.py            an example's TOML config: models, maps, wells, users, routes
  driver.py            Basin: steps every model one day at a time and keeps the ledger
  ledger.py            daily accounts, paired as sent and received, checked every day
  users.py             the water-user interface: roles, units, hourly-to-daily integration, the user's balance
  mapping.py           area maps between SWAT+ objects (catchments, aquifers, HRUs) and MODFLOW 6 cells
  gridmaps.py          maps on a MODFLOW 6 grid Basin Ledger didn't write, from geometry (HRUs, channel lines)
  units.py             unit factors to the ledger's SI units
  coupling/
    modflow6.py        MODFLOW 6 through libmf6 (modflowapi): prepare, set, solve, read, restart
    swatplus.py        SWAT+ project edits (places of use, point sources, pass-through aquifers), reads and sets
    recharge.py        SWAT+ aquifer recharge -> MODFLOW 6 recharge
    exchange.py        MODFLOW 6 river exchange <-> SWAT+ channels (two-way)
    existing.py        what Basin Ledger takes over in a MODFLOW 6 model, and what it leaves as the model's
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
- `swat:place_of_use:<name>`: a SWAT+ water allocation demand, taken the same day from its channel, or from its
  reservoir (`{name, reservoir = "res30"}`, Nov 25, 2026). A channel's is taken from its outflow after routing. A
  reservoir's is taken from its storage after the day's inflow and before its release. What the source couldn't
  supply is entered as a shortfall. The places of use are added to the run's copy of the project from the config
  (`[swatplus] places_of_use`). Several on one channel or reservoir share its point of diversion, and SWAT+ serves
  them in the config's order, so the first listed is the senior: it takes what it asks, and the next takes what it
  asks of the rest. Places of use on more than one channel or reservoir need swatplus-bmi 0.14.0 or later. SWAT+
  62.0.1's wallo_withdraw.f90 gives a place of use on any point of diversion but the first another's demand or none,
  which 0.14.0 patches, and an older swatplus-bmi is refused them at start-up (Nov 25, 2026).

  A user drawing from one only learns what it got after SWAT+'s day, but its return enters SWAT+ before that day. So
  its day scales by the share delivered: consumption, gain, storage, and the return with its loads and heat, as a
  withdrawal limit scales a day. Its return is then held overnight and leaves the next day (Andrew, Nov 25, 2026).
  The ledger books `user:<name>:delivered`, `:unmet`, `:delivered_share`, `:return_due`, `:return_sent` and
  `:return_held`, and closes two pairs every day: what it asked = delivered + unmet, and delivered = consumption met
  + return due + storage less gain met. So no user returns water its intake never delivered. Before this, the
  Virginia study's Lake Pelham intake got 8,697 m3 of the 870,385 it asked, and its plant still returned 829,326 m3.
  A return held on the run's last day stays in `:return_held`. A restart carries it. A user on a utility, or drawing
  from a well, returns as before;
- `swat:point_source:<name>`: a SWAT+ daily point source into its channel (`[swatplus] point_sources`);
- `ground:<name>`: a return into an aquifer cell (`[[ground]]`; Nov 25, 2026), for septic drainfields, a system's
  leaks and, later, land application.
  - Basin Ledger adds a WEL entry at the cell, in a package of its own (`wel_ground`), injecting the return the same
    day. MODFLOW 6 delivers it exactly, and its budget, and so the groundwater closure, take it in. The ledger pairs
    requested and delivered. The whole-basin balance books it as `returns_to_ground`, and the capture analysis as
    `from_returns`.
  - A house pumping its own well and returning to a drainfield near it is two entries in different cells: the
    well's withdrawal and the field's injection, each counted once.
  - Chosen over a return into an HRU's soil: a drainfield is below the root zone, so its water belongs to the aquifer
    cell, and swatplus-bmi sets no soil water. The cell's layer is the study's to give (the top active layer, as a
    rule).
  - With transport and salt (0.17.0), a return carries a concentration on its WEL entry (Mass through the
    human-system routes).
  - Over several cells (`cells = [[layer, row, col, weight], ...]` or a CSV; `[[wells]]` too; Nov 25, 2026): the
    route's request splits by weight into a WEL entry per cell, and its delivered is their sum. This serves a
    town's septic water and leaks across its service area, and domestic wells grouped by HUC-12. With transport
    (0.17.0) a well over several cells takes each cell's concentration by what it gave;
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
A user with its own wells or intake is routed directly.

**Water systems with customers (designed and built Nov 25, 2026, for the Virginia study's towns; `chain.py`, `tests/test_chain.py`, `tests/test_systems.py`).** Andrew's principle:
Basin Ledger is the link between the natural and human systems, so every water use goes through it, because a
scenario changes users. A town is one system: its sources, its customers, and where their water goes. A new
customer, a data center, adds to the system's demand and is drawn from the same sources. The `[[utilities]]` table
grows, and an old config means what it did (its `demand` a class of existing customers, `sewered` its sewered
share, the rest consumed).

- **Customers** (`customers = [{name, demand, sewered, septic, septic_to, consumed}]`). A class's demand is a
  constant (m3 s-1), 12 monthly values (m3 s-1, the same each day of the month) or a daily CSV (date, m3). Its water
  splits three ways, summing to 1:
  - `sewered`, to the plant;
  - `septic`, to the ground through `septic_to` (a `ground:` route);
  - `consumed`, outdoor use.

  Users on the system (`utility:<name>`) add their withdrawal to the demand, and their return goes to the sewer.
- **Losses.** `losses` is a share of production. `leaks = {share, to}` sends that share of the losses to the ground,
  and the rest is unaccounted water, booked as consumed.
- **Sources.** Either shares (a constant, or 12 monthly values that sum to 1 each month), so a reported history
  reproduces, or `priority` and `capacity` (m3 s-1): the first source up to its capacity, then the next. Shares
  can also come from a dated table (`shares = "file.csv"`: a date, then a column per source's route, each row
  summing to 1 and held until the next date; Nov 25, 2026). A system whose split changes over the years uses it:
  Culpeper's wells start in 2015, Lake Mooney's plant in 2018, and Motts Run's river and reservoir split moves
  year to year. A day before the table's first date is refused. Shares and
  priorities don't mix within a system. Production asked = demand / (1 - losses). A source can be any withdrawal
  route (`swat:place_of_use:`, `well:`, `offsite:`, an import), or `effluent:<system>`, another system's treated
  plant outflow, taken before it reaches its outfall (reuse: a data center supplied from a town's plant).
- **Shortfall.** What the sources delivered is known only after the models' day. So the system's day scales by
  production delivered / asked:
  - each class and user gets that share, and the rest is unmet;
  - losses, sewer, septic and consumption scale with it;
  - a user on the system scales as a user on a short intake does (0.13.0).

  A priority system can't make up a short source the same day; the shortfall is booked.
- **Returns leave the next day,** as a short intake's do (Andrew's choice for 0.13.0): the plant's outflow, septic
  water and leaks are held overnight. So yesterday's sewer is today's plant outflow, less what other systems take
  from it as reuse, which is known at the day's start. The plant's inflow beyond the customers' sewered water
  (inflow and infiltration) stays out by default, since the models already carry it. `infiltration` (m3 s-1) still
  adds it as water from outside where a study wants it.
- **Ledger** (`utility:<name>:`), every day, closing:
  - production asked and delivered by source, and the share;
  - each class's demand, delivered, sewered, septic and consumed;
  - losses as leaks and unaccounted;
  - the sewer, the plant's outflow, its reuse and its outfall, and what is held overnight.

  Pairs: production = delivered + losses; delivered = the classes' and users' delivered; each class's delivered =
  sewered + septic + consumed; losses = leaks + unaccounted; and the plant's outflow = outfall + reuse. A restart
  carries what is held.
- **Constituents** (0.17.0): the sources' concentrations mixed by what each delivered, each class's returns at a
  stated concentration or the mixed water's (Mass through the human-system routes).
- **Not built:** a short source made up by the next in priority the same day; a plant's own consumption.

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
| `catchment_cells.csv` | catchment (or SWAT+ aquifer, or HRU), cell, area_m2 | recharge in; drawdown and storage out |
| `flowpath_cells.csv` | flowpath, cell, length_m | stream-aquifer exchange (base flow) out, and in for losing reaches |
| wells in `basin.toml` | name, cell (or MAW connections) | pumping in; heads at wells out |

**Cells (built Nov 10, 2026).** A cell is named in the grid's own form, 0-based as flopy numbers it: `layer, row, col`
on a structured grid (DIS), `layer, icell2d` on a vertex grid (DISV). A map's columns say which (`object, layer, row,
col, ...` or `object, layer, icell2d, ...`), or give MODFLOW 6's reduced `node`; wells, head outputs and drawdown rules
take `[layer, row, col]` or `[layer, icell2d]`. The coupling reads the grid from MODFLOW 6's discretization (its shape
has three numbers for DIS, two for DISV) and turns a cell into its reduced node (`Modflow6.node`), refusing a cell in
the other grid's form, outside the grid or inactive, so a map made for one grid can't be read silently on another.
Before this, a column past the grid's edge wrapped to the next row's cells. Everything after the node (recharge, the
exchange, drains, groundwater into the soil, transport) was already by node, so it runs unchanged on either grid.
A vertex grid lets a study refine its cells around a site and its streams (a quadtree from 200 m to 25 m, say) without
a structured grid's long thin cells. An unstructured grid (DISU) isn't supported: its maps could give nodes, but its
wells and heads have no cell form yet.

The test (`tests/test_disv.py`) runs the two-way synthetic example on a vertex grid of 28 cells, 79 m refined to 40 m
around the well, with the recharge and river maps, the well, the head output and a drawdown rule all by icell2d and
chloride carried: MODFLOW 6's budget closes every day to 1e-7 of the gross flow and matches the ledger to 1e-12, and
the transport's mass closes. The same aquifer as a uniform DISV grid, cell for cell the structured one, gives the
structured run's ledger to 1e-7 relative: MODFLOW 6 reads a vertex grid's coordinates to about ten digits, so its cell
areas differ from the structured grid's by about 4e-10.

- **Recharge in.** Each catchment's recharge from SWAT+ spreads over cells by area of overlap, from its aquifer or,
  with `[recharge] source = "hru"`, from each HRU on the cells under its ground. Volume is kept exactly (built:
  `mapping.py`, `coupling/recharge.py`).
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
  AORC over each catchment, through `[weather]` (below).
- Water users get hourly air temperature and relative humidity (from AORC's specific humidity, temperature and
  pressure) at their site, through the `inputs` hook.
- MODFLOW 6 needs weather only if an example uses its ET or unsaturated zone packages.

**Weather into SWAT+ (built Nov 8, 2026; `[weather]`, `weather.py`).** swatplus-bmi writes a series into a
gage's daily records ahead of the day SWAT+ reads them (`set_weather`, 0.6), and into a run's files before it starts
(`project.write_weather`, 0.7.0). The run is then the one a weather file with
those values would give: SWAT+ makes the mean temperature, PET, heat units and stream temperature's air lag from it,
as it always does. Basin Ledger's side:

- **Series.** A CSV per variable (precipitation, maximum and minimum temperature, solar radiation, humidity, wind): a
  date column, then one column per gage, named as the gage's file in the project, in the file's units (mm/day, degC,
  MJ m-2 day-1, a fraction, m s-1). SWAT+ reads a day's weather before the pause that comes before it, and
  precipitation a day earlier still, so a run's first day (two for precipitation) can't be set through the BMI.
  For a fresh run the driver writes those days into the run's copy of the gage files before SWAT+ starts (decided
  Nov 9, 2026), then sets the rest after it starts, to the run's end, so the run is the series from its first day.
  It stops the run on a day past the series' end rather than leave SWAT+ on its files (decided Nov 9, 2026: a
  silent fallback to the files would mix two sources in one run).
- **Changes** (`[weather.change]`): a factor for precipitation, solar radiation, humidity (capped at 1) and wind, a
  shift for both temperatures, each one value or 12 monthly values (the delta-change method). A change applies to
  the series, or to the project's own records where no series is given, so a wetter or warmer case needs no files.
- **Scenarios.** A case from a saved warm-up keeps the warm-up's records up to its restart and writes its own from
  the first day it can (`run_case(..., change=...)` editing `cfg.weather`). A case with the warm-up's weather matches
  the unbroken run bit for bit.
- **Water users** on a SWAT+ station's weather (`forcing = {source = "swatplus"}`) read the series where it replaced
  the files, so the campus and SWAT+ see the same days.
- **Gridded data to gages** (`from_grid`): each gage takes its nearest cell, or the area-weighted mean of the cells
  nearest the HRUs whose station reads it (`hru_points`, from `hru.con` and `weather-sta.cli`). Fetching ERA5, AORC
  or a climate model's output belongs in a study; Basin Ledger takes the daily grid.

Limits, all SWAT+'s or swatplus-bmi's:
- A case from a warm-up takes its own weather from the day after its first (two days after, for precipitation):
  SWAT+ has read those days, and the saved state holds the warm-up's. The summary gives each variable's first and
  last day of weather.
- Only inside each gage file's period (SWAT+ generates the days outside it) and only for daily gages; missing days
  (-99) in a changed record are refused, since SWAT+ generates them.
- Gages are per file, not per HRU. A series changes a gage's values, not which gage a station uses, so finer weather
  needs more gages in the project's files (one per grid cell, say), with stations and HRUs pointing at them. Editing
  a project's stations isn't built: it is the project builder's job (the hydrofabric builder, or QSWAT+).
- Not built: weather from a model stepped with the basin (a component setting each day's record a day or two
  ahead through `set_value`).

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

- After SWAT+'s day, the driver reads each aquifer's recharge (`soil_water_table__recharge_rate`, m s-1), or each
  HRU's percolation (Recharge from HRUs, below).
- The area map spreads it over cells. A cell's rate is the sum over the objects covering it of rate × overlap area,
  divided by the cell's area. Volume is kept exactly.
- The rate goes into MODFLOW 6's recharge package for the same day.
- Users' well withdrawals are summed per well and set in the well package.
- In the ledger, recharge sent (rate × mapped area) equals recharge received (RCH `SIMVALS`), and well pumping
  requested equals pumped, every day.

In one-way mode SWAT+ keeps its own aquifer and baseflow, and MODFLOW 6 is a reading of what the same recharge does to
heads and to stream leakage. Each model keeps its own calibration. This is enough to show drawdown at a utility's wells
and when depletion reaches a stream. It is not enough to change the streamflow SWAT+ reports.

**Recharge from HRUs (built Nov 8, 2026; `[recharge] source = "hru"`).** SWAT+'s aquifer is lumped per subbasin, so
mapping aquifers to cells spreads recharge evenly within each subbasin: a cell under row crops on outwash takes the
same rate as one under forest on till. Bailey et al. (2025) map HRU recharge to cells by area. Basin Ledger does the
same with swatplus-bmi 0.6.0's `soil_water__percolation_rate`: each HRU's percolation out of the bottom of its soil
profile, the day just done.

- **The map** gives each HRU's ground in each cell (object, cell, area_m2; the cell as layer, row, col or layer,
  icell2d), from HRU polygons (flopy's GridIntersect) or, where HRUs are raster classes on a structured grid, from their
  pixels (`mapping.pixel_rows`, which splits each pixel exactly among the cells it overlaps and takes a pixel shared by
  several HRUs in their shares).
- **Volume.** Each HRU's rows are scaled to add to its area in SWAT+ (`hru.con`; `AreaMap.conserving`), so a cell
  takes each HRU's percolation by its share of that HRU's ground, and the volume sent is percolation times HRU area,
  exactly. The map places water; it never adds or loses any, so `recharge:by_the_map` is zero. The area check
  (`area_rtol`) compares the map's ground with SWAT+'s areas before scaling.
- **Against SWAT+'s aquifers.** SWAT+ 62.0.1 puts the same day's percolation into the HRU's aquifer: the HRU's
  recharge hydrograph is percolation times its area (`hru_hyds.f90`), and the aquifer takes the sum of what reaches
  it (`aqu_1d_control.f90`; the recharge delay is commented out in 62.0.1). Where every HRU recharges an aquifer, the
  two paths carry the same water each day, so the driver books SWAT+'s own aquifer recharge
  (`recharge:swat_aquifers`) and holds `recharge:made` to it (1e-5, for float32 reads and areas rounded to 1e-5 ha).
  A project whose HRUs send recharge elsewhere would fail that pair, and should say so. A project with no aquifers
  has nothing to hold it to.
- **What it changes.** Where the water lands within a catchment, not how much.
- **Not built:** an unsaturated zone between the soil profile and the water table (MODFLOW 6's UZF); SWAT+ 62.0.1
  has none either, with its delay commented out. Water back up from the water table into the soil is below.

**Groundwater up into the soil (built Nov 8–9, 2026, on swatplus-bmi 0.7.0;
`[groundwater_to_soil]`, `coupling/gwsoil.py`).** Recharge HRU by HRU sends water down. Where MODFLOW 6's water table
rises into an HRU's soil profile (wetlands, low ground beside streams, outwash with a shallow water table), water
should also come back up, so the soil and its plants can use it and a saturated profile sheds it as runoff. On
glacial-drift basins this is likely the larger effect: a shallow water table supplies wetland ET through the summer,
and one-way recharge can't represent it.

swatplus-bmi's `gw-soil` branch adds `soil_layer_water__groundwater_inflow_rate`: per HRU soil layer, m s-1 over the
HRU's area, set for the next day, at least zero, and refused in a project that runs SWAT+'s own gwflow. SWAT+ adds it
to the layers where gwflow puts its own (`gwflow_soil`, before the day's infiltration and percolation) and counts it
in `gwsoilq` (`hru_wb`, `basin_wb`). Then SWAT+ moves it as any soil water: to the plants, out as lateral flow, down
as percolation, or to the surface where the profile can't hold it. Read back after the day, it is what went in. The
same branch reads each layer as SWAT+ keeps it (bottom depth, water, saturation, field capacity), which sizes and
places the transfer.

```toml
[groundwater_to_soil]
rule = "drain"        # chosen Nov 9, 2026; or "gwflow". Required, no default
cap = "none"          # chosen Nov 9, 2026; or "room". Required, no default
package = "soil"      # a DRN (drain) or WEL (gwflow) package with one entry on each mapped cell
```

It needs SWAT+, MODFLOW 6, the HRU map (`[recharge] source = "hru"`) and the two-way `[exchange]`, so SWAT+'s
aquifers pass recharge through and groundwater is MODFLOW 6's alone. The rule and the cap have no default: they are
choices with consequences (below), so a config states them; the examples use the ones chosen, `drain` and `none`.

- **The water table under each HRU.** The HRU map gives each HRU's ground in each cell. Depth to water is MODFLOW 6's
  top less its head, so the depth and the head come from one model. A layer is below the water table where its
  bottom is deeper than that.
- **How much leaves a cell** (`rule`):
  - `drain`: a drain at the bottom of the soil profile under the cell (MODFLOW 6's top less the profile depths of the
    HRUs over it, by their shares of the cell), at the conductance in the package's input [A]. MODFLOW 6 computes the
    outflow inside its solve, at the heads as they converge, and books it in its budget.
  - `gwflow`: SWAT+ gwflow's rule. Each day, specific yield × the water table's height above each HRU's profile
    bottom × that HRU's ground in the cell, at the day's starting heads, taken by a well.
- **The cap** (`cap = "room"`): no cell gives more than the room left (saturation less water, as SWAT+ had it at the
  end of its day) in its HRUs' layers below the water table. For `gwflow` the well's rate is cut. For `drain` the
  drain's elevation is raised before the solve so its outflow at the day's starting heads fits. A head that rises
  during the day takes a little more, which SWAT+ sends to the surface. Raising the elevation between MODFLOW 6's
  outer iterations, as the exchange's supply limit moves stages, was tried first: the elevations crept by fractions
  of a millimetre each iteration, MODFLOW 6 stopped converging, and on some runs it aborted.
- **Into the layers.** What left a cell is split among the HRUs over it by their ground, and each HRU's part among
  its layers below the cell's water table by their thickness below it (gwflow's split). Water from a cell whose
  water table is below one HRU's profile (the drain sits at the cells' share-weighted profile bottom) goes to that
  HRU's bottom layer, and the ledger says how much (`gw_soil:to_bottom_layer`). The next day it goes to SWAT+.
- **Ledger.** `gw_soil:sent` (out of MODFLOW 6 today, the soil package's budget), `gw_soil:due` (into SWAT+ the next
  day, the same water) and `gw_soil:received` (what SWAT+ reports went in), with `due` = `received` checked every
  day to SWAT+'s 4-byte reals. The last day's water is in transit (`in_transit_gw_soil` in the summary) and is saved
  with a state, so a restart carries it. The soil package is a coupled link, not a boundary. The whole-basin balance
  counts it on both sides: out of MODFLOW 6, and into SWAT+'s land budget as `gwsoilq`. Capture
  (`analysis.capture`) has a `from_soil` part: pumping can draw the water table out of the soil.

Tested on Osu_1hru over the small aquifer with its land surface 0.5 m above the starting water table (the profile is
0.75 m), for a year, with each rule and each cap:
- what leaves the cells, by MODFLOW 6's budget, is what SWAT+'s layers receive the next day;
- the whole basin closes, and SWAT+'s printed `gwsoilq` matches what was received;
- the cap holds water back under both rules;
- a restart carries what is due, and matches the unbroken run bit for bit;
- a water table below every profile moves nothing, and the run matches the one without the link bit for bit.

What each choice does, and what was chosen. Decided Nov 9, 2026 (Andrew): `rule = "drain"` and `cap = "none"`. The
examples use them; the config still requires both, and `gwflow` and `room` stay as options. The drain departs from
SWAT+ gwflow's published rule (Bailey et al. 2020, `docs/prior-art.md`), which takes the whole saturated thickness
inside the profile each day by specific yield; the drain's conductance is a new parameter [A], calibrated with the
rest. With no cap, water the profile can't hold runs off at the surface.
1. **The rule.**
   - `drain`: head-dependent and implicit in MODFLOW 6's solve, so it is stable at a daily step. Its one parameter,
     the conductance, decides how fast the water table is held at the profile's bottom: high, and the cell can't rise
     above it; low, and water stays in MODFLOW 6 and reaches the streams instead. It is calibrated with the rest
     [A]. On the test, a year moved 71,300 m3 from 10 ha at a conductance of 10 m2 d-1 per cell.
   - `gwflow`: no new parameter, and what SWAT+ users know. It is explicit, though: it takes the whole saturated
     thickness inside the profile every day, at the day's starting heads, so a cell's head can swing between days.
     On the test it moved 89,300 m3, and its first day took 3,750 m3, the storage above the profile's bottom
     that had built up before the link.
2. **The cap.** `none` lets SWAT+ decide: what the profile can't hold runs off at the surface (saturation-excess
   runoff, the process at seeps and wetlands). `room` leaves that water in MODFLOW 6, where it reaches the streams as
   base flow, later. Both close the balance; they differ in timing and path. On the test the cap cut a year's water
   by about a quarter (drain) and four tenths (gwflow).
3. **Lag.** The water reaches SWAT+ the next day, as the stream exchange's does. Same-day would need SWAT+ to run
   after MODFLOW 6, which the day's order (SWAT+ first, its recharge into MODFLOW 6) doesn't allow.
4. **Cycling.** Water that comes up and percolates back down returns to MODFLOW 6 as recharge the same day (SWAT+
   adds it before percolation). The ledger shows both directions, gross: net changes the heads, gross is what a
   wetland's ET draws on.

Not built: a check at start-up that the flow model has no EVT or UZF package (SWAT+'s soil now takes ET from
groundwater, and a second ET would count it twice); an example should leave them out.

**Recirculation: gross and net recharge (built Nov 27, 2026, for the Virginia study; bookkeeping only).** With
groundwater up into the soil, part of what MODFLOW 6 sends up percolates again, the same day or later, and is counted
again as recharge. In the Rappahannock about 60% does: about 282 mm a year of gross percolation against 165-180 mm
net. Recharge alone, then, overstates what the land gives the aquifer. Nothing here changes a flow.
- **The whole basin** (`balance.basin_balance`): `recharge_gross` (what MODFLOW 6 received from the land),
  `groundwater_to_soil` (what it sent up, its soil package's budget), `recharge_net` (the difference) and
  `recirculation_share` (sent up over gross). The run's `summary.json` gives the same from the ledger
  (`recirculation`) when the link is on. Every existing key and ledger account is unchanged.
- **By HRU group** (`[outputs] recirculation = true`, `balance.recirculation(run, groups)`). The driver reads each
  HRU's day from SWAT+ after it runs (`hru_wb`: percolation, `gwsoilq`, plant and soil ET, lateral and tile flow,
  saturation excess, with precipitation, irrigation, ET and runoff for reference), sums them over the run in m3
  (`balance.HruWater`, carried by a restart), and writes `outputs/recirculation_hru.csv`. SWAT+'s average-annual
  HRU file can't serve: it is written at the end of `time.sim`'s period, which a coupled run rarely reaches, and a
  daily HRU print is gigabytes for thousands of HRUs. The reads need SWAT+'s daily records from the first year
  (`records_every_year`, which changes what SWAT+ prints, not what it simulates). `recirculation` groups the HRUs by
  land use (each HRU's `lu_mgt`), one by one, or by a dict, and gives each group and the basin its percolation,
  groundwater to soil and net recharge, in m3 and mm a year.
- **Where the water that came up went** is an estimate [A]. SWAT+ mixes a layer's water, so it can't say which drop
  percolated. The split takes the soil's water as well mixed: over the run, the water that came up leaves the soil
  in the proportions all the soil's water left it, by percolation, ET (plant and soil), lateral and tile flow, and
  saturation excess. The four parts add up to what came up. Saturation excess is only what SWAT+ records as `satex`,
  which SWAT+ 62.0.1 does only with its own gwflow on (`swr_satexcess.f90`); otherwise what a full profile sheds is
  in its surface runoff with the rain's, and that part reads 0. An exact split would need a tracer carried through
  SWAT+'s soil.
- Tested (`tests/test_recirculation.py`, Osu_1hru over the small aquifer with its water table in the soil, as the
  link's own tests): net = gross - groundwater to soil, and the summary agrees with the ledger; the HRU's sums agree
  with the ledger's (`gw_soil:received`, `recharge:made`) to 4-byte reads; the split adds up to what came up; the
  ledger is the same to the bit with the table on and off; a restart carries the sums exactly; without the link the
  recharge is all net. Without the table, the ledger is the same bit for bit as 0.18.0's, with the link on and off.
  On that test (a paddy whose water table is held in its profile), more comes up than goes down: 711 mm against 696
  mm, so net recharge is slightly negative; about 68% of what came up went back down, 32% to ET.

**Two-way: stream-aquifer exchange (built with RIV, and with SFR since Nov 11, 2026; `coupling/exchange.py`,
`examples/synthetic/two_way.toml` and `sfr.toml`).** MODFLOW 6 computes the exchange between streams and the aquifer,
through a river (RIV) or a streamflow routing (SFR) package, and that exchange replaces SWAT+'s aquifer return flow.

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
  river cells' gains and losses, its drains' seepage and its lakes' exchange) can't take more than the channel has
  without today's exchange, which won't be there tomorrow when this one is taken. The limit works inside MODFLOW 6's
  solve: after each outer iteration, a channel taking more than that at the heads as they are has its stages
  lowered by one depth, no entry below its bottom (`limit_to_supply`), and the day is done when the solve has
  converged and no stage moves more than `limit_tolerance` (1 mm by default). The ledger records the depth
  (`exchange:<channel>:stage_lowered_m`). The stages move for at most half the solver's `outer_maximum` or 30 rounds
  after convergence; the rest of the day solves them as they stand (0.10.1). Each day's search starts from yesterday's
  lowered stages (`limit_warm_start`, on by default from 0.10.1), saved with a restart: the depths are lowered from
  the day's own stages either way, so only the number of outer iterations changes (a 1,157-channel study basin's
  summer days went from 35-54 to 7-13).
  Two earlier forms failed on a study basin: a limit set before the solve, at the day's starting heads (wells and
  drains beside a river pull its heads down during the day), and a limit that counted the channel's own gain from
  the day before, so it lost on alternate days what it had gained.
- **Shortfalls.** What the channel still can't supply is entered in the ledger: water MODFLOW 6 took from a stream
  that, in SWAT+, wasn't there. Most of what remains is the one-day lag. Over 120 days on a 600 km2 glacial-drift
  basin the shortfall is
  16.3 million m3 without the limit (35% of the exchange) and 0.2 million m3 with it (0.7%); the run takes about
  twice as long. In the synthetic two-way example it is 186 m3 over three years without the limit and 76 m3 with
  it, of about 61,000 m3. Each day starts from its own stages, so a restarted run matches.
- **Through SFR (built Nov 11, 2026; `[exchange] package` an SFR package).** SFR computes each reach's stage from its
  flow by Manning's equation and its leakage at that stage, inside MODFLOW 6's solve, and a reach can't lose more than
  it carries. So it replaces both RIV's fixed stage and the stage-lowering limit, and suits headwater drains that go
  dry.
  - Each day, a coupled channel's water without today's exchange (what the limit uses) goes in as the inflow of the
    channel's top reach (`exchange:<channel>:sfr_inflow`). The reaches below carry it, less what they lose and plus
    what the aquifer gives. The channel's net exchange goes to SWAT+ the next day, as with RIV.
  - SFR's outflow from a channel's last reach isn't entered anywhere: SWAT+ routes the channel's water on, and is
    the only model routing water between channels. SFR's reaches therefore connect only within a channel. A reach
    flowing into a reach of another channel, or one not coupled, is refused at start-up (cut that connection when
    building the model from the SWAT+ channels), as is a channel with more than one top reach.
  - Also refused: diversions and the mover where it doesn't line up (water SWAT+ doesn't see; below, "the mover in
    the line-up rule"), SFR's storage option (kinematic-wave routing:
    water held in the reaches from day to day, which a restart would need), `[exchange] limit = false`, `[rivers]` on
    the SFR package, and transport with an SFR exchange (MODFLOW 6's SFT package, not built).
  - Manning's constant must be for the model's units: `UNIT_CONVERSION 86400` (or `TIME_CONVERSION 86400`) for meters
    and days, 1.486 × 86,400 for feet and days (to a thousandth: the 1.486 written for feet is 1.4859). MODFLOW 6's
    default, 1, is for meters and seconds; on it, a reach of the synthetic example a few centimeters deep came out
    eight meters deep. It is refused.
  - Maps give each reach's cell, or its 0-based number (`object, reach`), since several reaches can share a cell.
  - **Built from the channels (Nov 11, 2026; `coupling/sfr.py`).** `channel_reaches` cuts each SWAT+ channel's line,
    drawn downstream, into reaches on a structured or vertex grid (flopy's GridIntersect; a line back into a cell
    makes a second reach there), in order downstream. `add_sfr` writes the package the exchange takes: reaches
    chained within each channel only, Manning's constant for meters and days, no inflow of its own. `write_reach_map`
    writes the map. The streambed top, width, thickness, conductivity, roughness and layer are the study's [A],
    given per reach, per channel or once. Slopes come from the tops; a streambed rising downstream is refused (a DEM
    shows culverts, ponds and bridges, and correcting them is the study's), a flat one needs `min_slope` [A]. USGS's
    SFRmaker builds SFR from hydrography on structured grids, with DEM sampling and smoothing; a package from it, or
    from anywhere, is checked the same way at start-up. The synthetic example's package is built with the helper.
  - Without storage, SFR keeps no state between days, so a restarted run matches the unbroken one bit for bit
    (tested).
  - In the synthetic example (`sfr.toml`: `two_way.toml` with the river as ten SFR reaches of the RIV entries' bed
    conductance), the reach gains on 1,074 days and loses on 22 in three years, never more than its inflow; the
    shortfall from the one-day lag is 6 m3 of about 62,000 m3 exchanged.
- **Lakes at a fixed stage (built Nov 13, 2026; `[exchange] lakes`, `examples/synthetic/lakes.toml` and
  `lakes_sfr.toml`).** On a glacial-drift study basin the lakes are SWAT+
  reservoirs, each between a flowpath's channel and the next channel down, but MODFLOW 6 had none: the aquifer was
  asked to supply the lakes' water, and lake cells drained only through seepage drains. swatplus-bmi had no reservoir
  exchange variable then, so a lake couples through a channel, as the drains do. (swatplus-bmi 0.11.0 has one: a
  MODFLOW 6 LAK lake coupled to its SWAT+ reservoir is `[[lakes]]`, below.)

  ```toml
  [exchange]
  lakes = {package = "lakes", map = "maps/lakes.csv"}   # a RIV package; object, layer, row, col (or layer, icell2d)
  ```

  - The package is a RIV package of the lakes' own. Its entries are lake cells at a fixed stage, the lake's surface,
    set by the study in the package's input [A], with the lakebed's bottom and conductance [A]. Basin Ledger changes a
    stage only through the supply limit below, and puts it back after the day, as for the rivers.
  - The map gives each lake cell its outlet channel (the channel the reservoir flows into), which must be one of the
    exchange's coupled channels (the rivers' or the drains'); the lakes may use only some of them.
  - The lakes' net exchange per channel (positive to the stream) joins that channel's exchange the next day, as the
    drains' seepage does, and is booked as `exchange:<channel>:lakes`, the lakes' part of `to_stream`. It is in
    `due`, `applied` and `shortfall` with the rest. The package's coupled entries are a link, not a boundary flow.
  - **Supply.** A lake can't take more from its channel than the channel has (the `available` of the limit and of
    SFR's inflow). With a RIV exchange the limit works on the river and lake entries together: a channel taking more
    than it has has its river and lake stages lowered by one depth (the drains' seepage still adds to what is
    available). With SFR, which has no limit of its own here, the lakes get one: only their stages are lowered,
    against `available`. It doesn't count SFR's own reach losses from the same channel, so where both take water
    from one channel the excess over `available` is a shortfall, booked in the ledger. That should be rare: a lake's
    outlet channel carries the lake's outflow. The ledger records the depth as
    `exchange:<channel>:lake_stage_lowered_m` (with RIV, the same as `stage_lowered_m`).
  - Refused at start-up: a package that isn't RIV, the exchange's own package or the drains', `[rivers]` setting its
    stages, a mapped channel the exchange doesn't couple, and transport with lakes (their losses would need their own
    concentrations; not built).
  - In the synthetic example (`two_way.toml` with a four-cell lake east of the well at 88.5 m, lakebed 2 m deep, 10
    m2 d-1 per cell [A]), over three years the lake gains on 741 days and loses on 355 with RIV (4,900 m3 to the
    stream, 3,000 m3 from it), and 570 and 526 with SFR. The limit lowers the lake on 87 days with RIV and 224 with
    SFR. The shortfall, the one-day lag, is 518 m3 with RIV and 248 m3 with SFR, of about 65,000 m3 exchanged; no
    day's lake and reaches together took more than the channel had under SFR. A lake of zero conductance gives the
    run without lakes bit for bit, and restarts match the unbroken run bit for bit.
- **Stage.** With `[rivers] stage = "surface_depth"`, MODFLOW 6's river stage on day *d* is the river bottom plus
  the bed thickness plus SWAT+'s channel depth for day *d* (`channel_water__depth`, on SWAT+'s rating curve,
  swatplus-bmi 0.2.0). Otherwise it is fixed in the input. Under a NextGen replay the same setting uses t-route's
  depth.
- **Stage floor** (built Nov 27, 2026, for the Virginia study; opt-in, off by default). A river whose bottom was
  set below a surveyed water surface sits that far below it when SWAT+'s channel is shallow. In the Rappahannock
  build the main stems' bottom is the DEM's water surface less 2 m (1 m of channel, 1 m of bed), and SWAT+'s
  low-flow depth on those 100-190 m wide stems is about 0.4 m, so their stage was about 0.6 m below the surveyed
  surface: an extra drain through the lower basin. With `stage_floor = "<column>"`, a column of the rivers' map
  gives each cell's floor (m, the model's datum; blank for none), and the day's stage is the larger of the floor and
  bottom + bed thickness + depth, cell by cell.
  - It needs a stage set from the depth (`surface_depth` or `nextgen_depth`) and a map file (not `all:`), and the
    column; a value that isn't a finite number is refused. Cells no channel sets keep their own stage.
  - The exchange's supply limit still lowers a floored stage when the channel can't supply the loss, so a dry
    channel doesn't leak water it doesn't have. The floor is the stage before the limit.
  - Logged: a start note gives how many cells have a floor and how many of those floors are above the bed's top;
    the ledger counts the cells held at their floor each day (`rivers:stage_floor:cells`) and the raise summed over
    them (`rivers:stage_floor:raise_m`, m); `outputs/stage_floor.csv` lists each floored cell with its channel, its
    floor, the days held at it and the raise summed over those days (m days).
  - Tested (`tests/test_stage_floor.py`, the synthetic two-way example's river): a floor below every day's stage
    gives the run without it bit for bit (the new counts aside); a floor above it holds on the cells that have it and
    not elsewhere, as the ledger and the CSV say; the whole basin closes over two years with half the river floored;
    bad settings are refused. Without `stage_floor` the ledger is the same bit for bit as 0.18.0's, with and without
    a stage from the depth.
- **Timing.** The exchange is explicit with a one-day lag. Volume is conserved (ledger pairs): each day's exchange is
  due in SWAT+ the next day, and the last day's is reported as in transit. Timing shifts by a day, which is small
  against the time scales of groundwater response. Iterating within the day by branching SWAT+
  is possible but costs a SWAT+ day per iteration. It will be tried if the lag matters for an example's question.
- **Mapping.** Channels map to cells by length (`object, layer, row, col, length_m`, or `object, layer, icell2d,
  length_m` on a vertex grid), or to SFR reaches by number (`object, reach`). Conductance comes from the MODFLOW 6
  model, not from the map.

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
  - Built: intake concentrations as an input to user components
    (`water_user_intake_water_<species>__mass_concentration`), set each day from the aquifer at a user's well, and a
    return's heat booked at its point source (`route:<point source>:heat_m3C`). Nitrate and soluble P cross every
    SWAT+-MODFLOW 6 link with their mass (Water quality across links, stage 2). Still needed: a return's loads to its
    point source (swatplus-bmi has `point_source_water_<sp>__mass_flow_rate`; the routing isn't built).
- **Whole basin.** For each example, a basin-wide check closes over the run:
  - precipitation = evapotranspiration + streamflow out + change in soil, aquifer and channel storage + net user
    consumption + transfers out of the basin;
  - the same for each constituent's mass, from SWAT+'s and MODFLOW 6's own budgets plus the ledger.

## Conservation across links (Andrew, Nov 6, 2026; written before the fixes)

The models must never fight. Water may enter from outside (a supply, recharge where no surface model runs, a fixed
head), but every crossing of a model boundary is an account in the ledger, and no link makes or loses water. Sent =
received is necessary, not sufficient: the first coupled run of a study basin (600 km2 of glacial drift) balanced
every pair and still had four leaks, each of which the pairs hid:

| Leak | Size in that run | Why the pairs passed |
|---|---|---|
| Recharge mapped onto more area than the surface model makes it on (lakes in the catchments; a catchment outside SWAT+ given the mean rate) | 28 million m3, 9% of recharge | "sent" was figured on the mapped area, not SWAT+'s |
| Losing reaches taking water the channel didn't have (headwater channels left nearly dry once SWAT+'s aquifers pass recharge through) | 84 million m3, 30% of recharge | booked as a shortfall, while MODFLOW 6 kept the water |
| A starting state not consistent with the coupling (a steady state solved with every reach able to lose) | 49 million m3 of storage drained | storage is inside MODFLOW 6's budget, not a link |
| River cells outside the coupled channels at a fixed stage | an unlimited source or sink | not mapped, so not in the ledger |

Rules, each to be built and tested:

1. **One owner per process.** Each flux is computed by one model. The other models' copy of it is switched off or
   masked: set through BMI each day to the owner's value or to zero (swatplus-bmi's settable inputs, as the soil
   temperature in a study's soil hold). SWAT+'s aquifer return flow, revap and deep recharge are MODFLOW 6's once it is
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
   booked as `exchange:<channel>:seepage`. Without it a drain is a boundary flow, water leaving the basin. Lakes at
   a fixed stage are coupled the same way, through their outlet channels (`[exchange] lakes`, Nov 13, 2026), and
   booked as `exchange:<channel>:lakes`.
6. **The start is consistent with the coupling (built Nov 6, 2026).** `Basin.spin_up` runs the first year again and
   again, SWAT+ from its own start and MODFLOW 6 from the heads the last cycle ended with, until the year's net
   release from storage is within a tolerance of its recharge; `Basin(..., heads=...)` starts the run there. The
   ledger's `groundwater:storage_release` shows any drift that remains.
7. **Each model is checked for the settings coupling needs, at start-up (built Nov 6, 2026).** Below; a run with a
   wrong setting stops before its first day. A start that can drift (soils starting empty, a warm-up whose soils
   haven't settled) is a warning, not a stop (Nov 7, 2026).
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

Warnings, not stops (built Nov 7, 2026). A run starts anyway, but its first years may drift:

- **SWAT+ soils that start empty.** `sw_frac` is 0 in `soil_plant.ini` for the HRUs that use it. Basin Ledger warns
  at start-up (`start_warnings`, a `StartWarning`). A restart doesn't warn, since its soils come from the saved
  state.
- **A warm-up whose soils haven't settled.** `scenarios.warm_up` writes `warm_up.json` beside the saved state. It
  holds the basin's soil water change over the warm-up's last year (all of it, if shorter) as a share of that
  time's precipitation, from SWAT+'s `basin_wb_day` (`balance.soil_water_trend`). Above 12% it warns [A]. In
  settled runs of a study basin, weather alone moves a year by up to 10.4%. Its one-year warm-up from empty soils
  moved 15%, and the synthetic example's 200 days moved 30%.

This is a coarse screen. It catches a warm-up still filling from empty, but not a slow drift of a few millimetres a
year. That shows only when the same years run again.

A study basin's model is the case behind both checks. Its 20-year base case wetted for decades: mean soil water rose from 210 to
470 mm, and flow at the gage rose with it. Stock SWAT+ 62.0.1 alone does the same, and swatplus-bmi with nothing set
matches it byte for byte, so the cause is the setup, not the tools. The project's soils start empty, and its
parameters were calibrated over the two years after a one-year warm-up, while the soils were still filling. With
those parameters the soils settle far above field capacity, and the gage's flow settles at about twice what was
observed. The fix is the study's: settle the soils before calibrating, then calibrate.

**Checked before a run (built Nov 11, 2026; `basinledger check basin.toml`, `Basin.check`).** A run stops at its
first problem, and some problems show only on the first day, once MODFLOW 6 has read its packages. `check` starts
the models as a run would, has MODFLOW 6 read its first day, builds every link and lists every problem it finds,
without solving a day: SWAT+'s settings, parameters for models not run, the recharge map and its areas, the wells',
head outputs' and rules' cells, the routes, users and rules, the exchange and its maps (SFR's checks included), the
rivers, groundwater into the soil and transport. A part with a problem is listed and left out, and what depends on it
isn't listed again (a user on a well whose cell is wrong); a problem that leaves a model unable to start ends the
check there. The start-up's warnings come back as notes. A study runs it after each change to a model or a map,
before a calibration's runs. The checks are the run's own code, so a check and a run can't disagree. Since Nov 13,
2026 `check` takes the run's start-up hooks (`prepare`, `parameters`, `wrap_surface`, `heads`; the builder's
arguments are in the config), so a study checks what it runs. Before, a study whose prepare hook rewrites SWAT+'s
`calibration.cal` was told that the project's stored file turned `flo_min` back on.

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
  aquifer would replace, SSM's AUXMIXED and file sources, and decay (not yet in the ledger; built in stage 2).
- **Restarts** save and restore each species' concentrations with the heads, bit for bit.
- **Not yet:** the steady-state start (the tests start from a stated concentration), concentrations at the river and
  drain cells reported per channel for stage 2, and the stage 2 links below (all built in stage 2, Nov 22).

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

**Stage 2, built (Nov 22, 2026): nitrate and soluble P through every link.** `coupling/surface_mass.py`, on
swatplus-bmi 0.13.1 as released (no new variable was needed), `[transport.<species>] surface/recharge/stream`, and
`examples/synthetic/nitrate.toml` (both species, HRU recharge, rivers that gain and lose, drains, the soil link, a
well). Phosphorus and nitrate take the same wiring; P sorbs (linear, Kd [A], by zone if given), nitrate decays (first
order [A]).
- **Recharge, HRU by HRU.** Each HRU's mass leached below its profile on the day just done, read through
  swatplus-bmi's generated reads of SWAT+'s own daily records: nitrate `hru_pw:percn` (`percn(j)`, the bottom
  layer's percolate in nut_nlch.f90, which hru_hyds.f90 also puts in the HRU's recharge hydrograph, `hd(2)%no3`) and
  soluble P `hru_ls:lchlabp` (`hls_d%lchlabp`, nut_solp.f90). Both kg ha-1. The mass goes over the cells by the same
  conserving map as the water, so a cell's recharge concentration is the mix of the mass and water its HRUs send it.
  Ledger: `recharge:made` (SWAT+'s leached mass) = `recharge:sent` (on cells with water) + `recharge:no_water` (a
  recharge package can't carry mass without water; zero on the synthetic run), sent = received; nitrate's made is
  held to SWAT+'s own aquifer recharge mass (`aquifer_water_no3__recharge_mass_rate`) every day, as the water is.
  With `[recharge] source = "aquifer"` nitrate's mass is per aquifer, every cell under it at one concentration, and a
  warning says so (`AquiferRechargeMass`). Soluble P has no aquifer path: SWAT+ 62.0.1 leaves leached P out of the
  recharge hydrograph (it leaves the model in SWAT+ alone), so it needs HRU recharge and is refused otherwise.
- **The reads are NaN in print.prt's skipped years** (swatplus-bmi's `out_made`: SWAT+ makes its daily records only
  after `nyskip`). A run with a species SWAT+ carries sets nyskip to 0 in the run's copy
  (`swatplus.records_every_year`): only what SWAT+ prints and sums changes, not what it simulates.
- **The aquifer to the streams.** A coupled channel's rivers' gains and drains discharge at their cells' concentrations
  (GWT's outflow). Where the channel's net water for tomorrow is a gain, the gain's mass set in SWAT+ tomorrow
  (`channel_water_<sp>__groundwater_mass_exchange_rate`, with the water) is that discharge less what the channel's
  losing entries took the same day at SWAT+'s outflow concentration for the day just run, not below zero: SWAT+ gets
  the net, as it gets the net water.
- **The streams to the aquifer.** swatplus-bmi takes a loss's mass itself, at the day's outflow concentration
  (`exch_channel`: `rto * ht2%no3`), and refuses a mass set with a loss. The loss's water is solved in MODFLOW 6 the
  day before SWAT+ applies it, so a concentration set on the river's auxiliary that day can't be what SWAT+ later
  takes, and the link wouldn't close. So the loss's water enters at no concentration, and the next day, once SWAT+
  has applied the exchange, what it reports it received (negative where its loss took mass) is set against what the
  aquifer sent the channel the day before; the difference goes into the aquifer through a mass source on the
  channel's cells (GWT SRC, `src_<package>`, added by `add_transport(..., stream_mass=[...])` on the rivers' and
  drains' cells, entry for entry), on the entries that lost water, by the volume each lost (on those that
  discharged, by their mass, if none lost). Ledger per channel: `stream:<c>:due` (what the aquifer sent yesterday)
  = `in_swat` (what SWAT+ reports it kept) + `returned` (back into the aquifer), and returned = `returned_received`
  (GWT's mass source), every day, exactly. The loss's mass reaches its cells a day after its water [A], the
  exchange's own lag; a loss that takes the whole of a day's outflow takes all its mass. The channel's outflow
  concentration is booked (`stream:<c>:concentration`).
- **Lakes at a fixed stage** (`[exchange] lakes`): with the stream side a lake is its outlet channel's, as its water
  is: its gains carry the cells' concentrations into the channel, its losses the channel's mass back, through the
  lake package's own mass sources. Still refused with a stated stream concentration (stage 1's refusal).
- **Into the soil.** The soil link's water leaves each cell at the cell's concentration (GWT's outflow), is split among
  the HRUs' layers as the water is (`GroundwaterToSoil.mass_split`, the same rows and depths), and is set in SWAT+
  the next day with the water (`soil_layer_water_<sp>__groundwater_inflow_mass_rate`). Ledger: `gw_soil:sent` (GWT)
  = `split`, and `due` = `received` (SWAT+'s read-back).
- **SWAT+'s aquifers.** In two-way mode they pass recharge through. Their nitrate leaves only with return flow and
  seepage (`conc_no3 * flo`), both off, so the recharge nitrate stays in their store as the water does. Their soluble
  P didn't: SWAT+ 62.0.1 sends the channels `aqu_d%minp = hin%flo * aqu_dat%minp / 1000` (`hd(1)%solp`,
  aqu_1d_control.f90) whether or not water flows, so with `sol_p` above zero (aquifer.aqu) P reached the channels with
  no water, in every two-way run so far, transport or not. `aquifer_pass_through` now sets `sol_p` to 0, and a
  project or `calibration.cal` (`minp_init`) that sets it is refused. The ledger holds both outflows at zero every day
  with the stream side on (`swatplus:aquifer_no3_outflow`, `swatplus:aquifer_solp_outflow`, from `aquifer:no3_lat`
  and `aquifer:minp`). Groundwater P therefore comes only from MODFLOW 6's GWT, and no SWAT+ aquifer P path runs beside
  it.
- **Decay** (`add_transport(..., decay=, decay_sorbed=)`, MST first order, d-1 [A]): booked as `transport:<sp>:decay`
  (MODFLOW 6's RATEDCY and RATEDCYS), and the closure includes it. Zero-order decay is refused. Tested: the decay
  booked is the rate times the day's end dissolved and sorbed mass, and the day closes.
- **Sorption** (linear, MST): the sorbed store's release is booked on its own (`sorbed_release`, RATESRB) inside the
  storage release. Tested on steady flow: the change in what the cells hold dissolved (porosity x volume x C) and
  sorbed (bulk density x Kd x volume x C) is the storage release, and the sorbed part the sorbed release, every day.
  With transient flow MODFLOW 6 moves the mobile store's volume with the flow model's storage, so there the closure is
  what's checked.
- **The steady-state start** (`Basin.transport_steady_state(cfg, days, heads)`): a coupled run, its flows between
  cells and its inflows' water and mass per cell averaged over it (`Transport.track`), and each cell's concentration
  at steady state under them: c = mass in / (water in + decay coefficient), the decay coefficient the rate times the
  cell's dissolved and sorbed capacity. Advection, mixing and decay, no dispersion [A], solved by sweeps from the
  highest head down. With SWAT+'s mass the start is under the run's own flow-weighted mean recharge concentration.
  Tested: uniform recharge gives its concentration everywhere, and with decay it is the transport model's own answer
  after 22 turnovers (upstream advection, no dispersion), to 1e-6. Start a run from it with `Basin(...,
  concentrations=...)`. Checking it against public well samples is the study's.
- **The transport solve's residual criterion** is 1e-8 kg d-1 (1e-10 in stage 1): at 1e-10 round-off held the
  residual above it for good on Plaster Creek's 13,215 cells (day 47) and on the synthetic grid with sorption.
- **What it costs** (Nov 22, on a machine running a 10-core calibration beside it, load 15-20 on 14 cores, so read
  these as upper bounds). The synthetic run (100 cells, one HRU, one channel; 730 days, best of two): 4.9 ms a day
  without transport, 9.1 with nitrate, 9.0 with soluble P, 12.6 with both; there the added time is Basin Ledger's
  bookkeeping and the BMI reads (about 2 ms per species a day: the budget's 60-odd reads of MODFLOW 6's memory, the
  day's concentrations and mass sources, SWAT+'s reads and sets), not the solve. Plaster Creek's screening
  groundwater model (13,215 active cells in 3 layers, NextGen replayed, stage 1's stated concentrations; 365 days):
  50 ms a day without transport, 272 with nitrate (decay), 251 with soluble P (sorption), 449 with both: about 200
  ms per species a day, all but 1.3 ms of it in MODFLOW 6's transport solve (TVD advection, the solver's settings in
  `add_transport`). Sorption costs no more than decay. On a large basin each species is a second groundwater solve
  a day; the solver's settings (`ims`) and the advection scheme are the levers.
- **Not in stage 2:** a user's return loads to SWAT+'s point sources (built in 0.17.0 with the water systems'
  outfalls: Mass through the human-system routes); nitrate or P in irrigation pumped from wells (its mass leaves the aquifer with the
  well and isn't put on the HRU: swatplus-bmi has no way to set it); SFR, several groundwater models and
  multi-aquifer wells with transport (refused, as in stage 1); SWAT+'s reservoirs coupled to MODFLOW 6's lakes
  (`[[lakes]]`: LAK's exchange with the aquifer is booked as a boundary, its mass not carried to the reservoir); the salt ions (SWAT+'s salt module can't run with aquifers, above).

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

## Salt: chloride and weathering's ions through land, aquifer and streams (designed and built Nov 22, 2026)

For the Virginia study (the Rappahannock above Fredericksburg), which calibrates chloride and specific conductance at
many sites. Sources there: road salt (winter, on roads and developed land; in the stream fast through runoff and
slowly through groundwater, a legacy of decades: Kelly and others 2008, Kaushal and others 2005), weathering in
groundwater by geology, sewer plants and atmospheric deposition. SWAT+ 62.0.1's salt module can't run with aquifers
or point sources (Water quality across links, above), so Basin Ledger carries chloride itself, as it carries stream
temperature, and MODFLOW 6's GWT carries it in the aquifer (stage 1). Unlike the temperature, the balance runs inside
the coupled run, each day: the recharge's concentration comes from the land the same day, and the streams' losses
carry the streams' concentration into the aquifer, so neither can wait for the run's end.

Chloride is a conservative mass (kg): no reactions, no plant uptake, no sorption [A]. It is `salt.py`, asked for with
`[salt.<species>]` beside `[transport.<species>]`. Any number of conservative species run together, each with its
own GWT model, balance, ledger and closure, sharing the day's SWAT+ reads.

**Two species for the Virginia study (Andrew, Nov 22).** Late-summer low-flow conductance in the Rappahannock follows
geology: Blue Ridge about 79 uS/cm (lower quartile 60), Piedmont 103 (81), diabase 88 (80), hornfels 107 (91), the
Culpeper basin's sedimentary rock 467 (372). That background is Ca, Mg, HCO3 and SO4 from weathering, not chloride.
So weathering's ions are a second species (`[salt.weathering]`, their mass as one total, no rain or road salt), with
GWT's weathering source by zone, and conductance is a linear combination of the species (4, below). Chloride's own
weathering source stays available but is expected near zero.

**1. Land: a well-mixed store per HRU.** Each day, after SWAT+'s day:

- **Inputs.** Atmospheric deposition, a concentration in precipitation times the HRU's precipitation (snow included)
  [A]; road salt (below); land sources (below); and groundwater up into the soil, yesterday's water from each cell
  at the cell's concentration, split among the HRUs over it exactly as the water is (`gwsoil.after_solve`'s rows).
- **Road salt,** one of three ways [A]. A rate (kg per m2 of road a day) times the HRU's road area (a share, a number
  or a table by HRU) on snow or ice days: SWAT+'s snowfall of at least `snowfall_mm`, or precipitation of at least
  `precipitation_mm` with the HRU's minimum temperature at or below `tmin_c`. A series of dated rates (a DOT's
  records, say). Or, the study's form, an annual amount per m2 of the HRU by calendar year (`hru, year, kg_m2`; the
  USGS 1 km grid of Bock and others 2018, whose years are calendar years, 1992-2019), spread over that calendar
  year's snow and ice days in proportion to each day's precipitation (`spread = "precipitation"`; `"equal"` the
  option). A winter's two halves come from two calendar years' amounts: December's salt from the year it falls in,
  January's from the next. Precipitation weighting is chosen because salt applied per storm scales with the storm's
  accumulation and length; spread equally, a trace day would take as much as a blizzard. To spread a year's amount
  whole, its days must be known before it starts, so they come from the weather SWAT+ reads, not from its daily
  output: the HRU's station's `.pcp` and `.tmp` records, a day's snowfall being its precipitation at a mean
  temperature at or below the HRU's `fall_tmp` (`snow.sno`, as SWAT+ divides it without elevation bands). A year
  with no snow or ice day takes its amount equally over its frost days (minimum at or below `tmin_c`); with none
  either it isn't applied, with a warning. Years after the table's last take each HRU's mean of its last
  `repeat_last` years, and before its first the mean of its first `repeat_first` [A]; without them such a year is
  refused. (An `hru, kg_m2` table times a `factor` table by year is the other form.) The calendar for 33,525 HRUs
  over 28 years takes 8 s at start-up and reading the 24 MB table 1 s.
- **Land sources,** each named, with its own ledger line (`land:source:<name>`), into the soil store [A]: an amount
  in kg per m2 of the HRU a year, as a number, a table by HRU (`hru, kg_m2_yr`), by SWAT+ land use (`landuse,
  kg_m2_yr` with `by = "landuse"`), or by calendar year (`hru, year, kg_m2`, the form the study's potash, livestock
  and septic estimates come in, with `repeat_last` and `repeat_first`). Its timing: `uniform` (septic effluent, a
  constant input), 12 monthly weights, or `fertilizer`: the year's amount on the HRU's first fertilizer day of the
  calendar year, read through swatplus-bmi (`hru_nb`'s `fertn` and `fertp`, kg/ha applied that day), or on 31
  December if SWAT+ applied none that year; potash goes on with the year's first fertilizer. Rural low-flow chloride
  in the Rappahannock is 6-8 mg/L, several times rain concentrated by evapotranspiration (about 1-2 mg/L), which is
  why these sources are needed.
- **Multipliers.** Each input's mass belongs to a group (`precipitation`, `road_salt`, `point_sources` for the stated
  concentrations, and each land source's `group`, its name unless given), and `multipliers = {group = factor}`
  scales it where it enters, so the ledger books the scaled mass. A calibration fits them as parameters of model
  `salt`, target `<species>:<group>` (the study's four: road salt, agricultural, septic, sewer plants), set at
  start-up over the config's.
- **Outputs** leave in proportion to water at the store's concentration: the water SWAT+ sends to the stream that day
  (surface runoff, lateral flow and tile flow: `hru_wb`'s `surq_gen`, `latq`, `qtile`) and percolation
  (`soil_water__percolation_rate`, the recharge link's own variable, so the mass and the water are read alike).
- **The store's water** is the profile's water at the end of the day, the wilting point's included (chloride is
  dissolved in all of it) [A], plus the water SWAT+ holds back in its surface and lateral lags and an HRU wetland
  (`sw_final`, `wpmm` summed, `lagsurf`, `laglatq`, `wet_stor`). The concentration is the day's mass (yesterday's,
  plus the inputs) over that water plus the day's outflows: an implicit, fully mixed step, so the mass closes by
  construction and evapotranspiration concentrates it (the chloride mass balance's evapoconcentration).
- **Why SWAT+'s printed columns, read through BMI.** swatplus-bmi reads every column of `hru_wb` (and `hru_pw`,
  `channel_sd`) for the day just done, as SWAT+ prints it, without printing (its "every printed daily column"); the
  temperature's after-the-run reading of the files can't serve a balance that the aquifer needs the same day. The run's
  copy of `print.prt` starts printing at year 0, since the records are made only after `nyskip`'s years.
- **A surface store, optional** (`surface_store = true`; off by default; on in the synthetic example as the way a
  snowy basin runs, and for Virginia). One store per HRU
  mutes the winter flush: road salt is diluted into the profile's 200 to 400 mm of water before any reaches the
  stream. With it, deposition and road salt land on a surface store (the snowpack and the road), well mixed with the
  snowpack's water; each day the share of it that leaves as liquid (rain plus melt over that plus the snowpack left)
  is released, and of that the share that ran off (`surq_cont` over rain plus melt) goes to the stream that day and
  the rest into the soil store [A]. Surface runoff then carries the surface's mass and the soil store gives only
  lateral flow, tile flow and percolation. Land sources go to the soil store either way. INCA-Cl (Jin and others 2011, to check) keeps a fast and a slow store the
  same way.

**2. Groundwater.** The existing GWT species (`coupling/transport.py`), its recharge concentration each day the land's:
each HRU's percolation mass placed on the cells under its ground by the same conserving map as its water, so a cell
under two HRUs takes their flow-weighted mean, and the mass sent is the land's exactly (the recharge map needs
`source = "hru"`). Weathering: a mass source per cell by geologic zone (GWT's SRC package, `salt.add_weathering`,
kg per m3 of aquifer a year) [A], chosen over a background recharge concentration because weathering adds solute
with the water's time in the rock, which a source in the aquifer's volume carries and a recharge concentration
doesn't: a recharge background would put as much into a wet year's recharge as a dry year's, and none into slow, old
groundwater. For chloride it is usually small (crystalline Piedmont and Blue Ridge rock holds little chloride); the
major ions that weathering does add are the second species, `[salt.weathering]`, with their own source by zone. The
rate is per m3 of saturated aquifer: the builder writes each cell's rate at full saturation, and each day Basin Ledger
sets each entry to that times its cell's saturated share at the start of the day [A] (MODFLOW 6's NPF saturation at
the end of the day before, kept in the salt state for restarts; on a run's first day, from the starting heads:
`(min(h, top) - bot) / (top - bot)` in a convertible cell, 1 in a confined one). A dry cell weathers nothing. The
static rate put the full source into the Virginia study's dry layer-1 cells (24,842 of 111,897 dry at its steady
start, kept active with no water by Newton) and the weathering solve diverged on day 0 (Nov 22, 2026). The share is
the day before's because the flow is solved after the transport's sources are set, the same one-day lag as the
streams' losses. Wells as built (the user's intake follows the aquifer at its well).

**3. Streams: a mass balance by water body, upstream first.** Each SWAT+ channel and reservoir is a well-mixed store
(kg), in the order `chandeg.con` and `reservoir.con` route them:

- **In:** the upstream bodies' outflow mass, by their routing fractions; the land's runoff, lateral and tile mass
  through the routing units (`rout_unit.def`, `.ele`, `.con`; each HRU's mass in full, split as SWAT+ splits its
  water); rain on the channel or reservoir at the precipitation's concentration; point sources, a direct user's
  returned mass where its component reports the species, the rest of the point source's water at a stated
  concentration [A] (a sewer plant's effluent, say); and the groundwater's gains (rivers, drains and lakes at a fixed
  stage), each at its cell's concentration from GWT, arriving with their water the next day.
- **Out,** at the mixed concentration over SWAT+'s end-of-day volume plus the outflows: the outflow downstream (or
  out of the basin at the outlet), withdrawals (SWAT+ places of use), a reservoir's own seepage, and the water the
  aquifer took the day before through the channel's rivers and lakes (which SWAT+ takes from the channel today, the
  exchange's one-day lag). Evaporation takes no mass.
- **The losses' mass crosses the next day, exactly.** MODFLOW 6 takes a losing reach's water on day d, before the
  stream's day d+1 is known; SWAT+ takes that water from the channel on day d+1. A concentration set on the river
  package on day d (the channel's that day) would put mass into the aquifer that the channel must give up on day
  d+1 at whatever concentration it has then: in the synthetic run, where a pumped well makes the river lose most of
  its flow, the channel was short of that mass by 46% of the losses' mass. So the water enters the aquifer with no
  chloride (the river's and lakes' auxiliary set to zero), and its mass, each entry's water at its channel's mixed
  concentration on day d+1, enters the same cells on day d+1 through a GWT mass source with one entry per river and
  lake entry (`salt.add_stream_sources`, `stream_losses`, built by the example's builder). The stream gives up and
  the aquifer takes the same mass the same day; the aquifer's cells get it a day after the water, the same lag the
  gains have the other way. A lake at a fixed stage takes its outlet channel's concentration [A].
- A dry body keeps its mass until water returns. Conservative, no reactions.

**4. Conductance.** Specific conductance (uS/cm) = a constant + each species' slope x its mixed concentration (mg/L),
fitted per study [A] (`[salt.conductance] base, slopes = {chloride, weathering}`; chloride alone gave r2 0.87 on
Mountain Run). An output only (`conductance_us_cm:<channel>`), never fed back.

**5. Ledger and checks (kg).** Every source, store, link (sent and received) and outflow, as basin totals:
`salt:<sp>:land:*` (precipitation, road salt, each source, from groundwater, to streams, recharge, storage), `salt:<sp>:streams:*`
(from land, precipitation, point sources, gains, losses, withdrawn, outlet, reservoir seepage, storage),
`salt:<sp>:aquifer:from_streams` (the stream-loss source), the links in transit overnight (`transit:gains`,
`transit:soil`; `transit:losses_m3`, the lost water whose mass follows), and GWT's own accounts
(`transport:<sp>:*`, whose recharge sent is now the land's mass, and whose leakage carries none).
Held every day:
- the land's and the streams' closure, each a share of the day's gross mass, within 1e-9;
- land to streams sent = received, recharge sent = GWT's received (the existing pair), the streams' losses = the
  stream-loss source's mass, the soil link's mass out of GWT = the HRUs' share of it;
- the whole system: land + streams + in transit + the aquifer (from GWT's storage release) change by the sources
  (deposition, road salt, land sources, point sources, GWT's boundaries and weathering) less the sinks (the outlet, withdrawals,
  wells, reservoir seepage), within 1e-6 of the day's gross mass flow, as water is held.

Per channel, each species' concentration (`<species>_mg_l:<channel>`) and the conductance are series in
`series.csv` for `[outputs] channels`.

**Not built,** each refused at start-up: reservoirs coupled to LAK lakes (`[[lakes]]`), irrigation from wells
(`[irrigation]`), several groundwater models, an SFR exchange (needs SFT), NextGen, a MODFLOW 6 model not in meters
and days, SWAT+ water allocation other than Basin Ledger's places of use, and routing to anything but channels,
reservoirs and routing units (runon to HRUs). Not refused but not built: a user's intake concentration from a
channel (its withdrawal leaves the stream at the stream's concentration, booked as withdrawn; the user isn't told
it).

## Replay: the salt again from a recorded flow run (designed and built Nov 22, 2026)

A salt calibration changes only salt: input multipliers, sources, weathering rates, porosity, dispersivity, starting
concentrations. None of them changes the flow: MODFLOW 6 solves flow before transport and never reads it back, and
salt sends nothing to SWAT+ (Salt, above). So the flow side of a long, heavy run can be recorded once and the salt run
again from the record as often as a calibration needs, cheaply and in parallel, without SWAT+ or MODFLOW 6's flow
model. `replay.py` does it.

**Recording.** A normal run with `[replay] record = "<folder>"` saves, after each day, what transport and the salt
balances read from the flow side, and nothing else:

| What | Per day | Why |
|---|---|---|
| MODFLOW 6's face flows (FLOWJA) | the upper triangle, one value per connection; the lower is its negative, checked every day | GWT's advection, and which cells rewet |
| each cell's saturation (NPF SAT) | one per cell | GWT's water volume, dispersion's dry cells, weathering's share |
| specific discharge (NPF SPDIS) | three per cell, unless `specific_discharge = false` | dispersion's tensor (its only use) |
| storage (STO STRGSS + STRGSY) | their sum, one per cell (both in float64) | GWT's old water volume reads only the sum |
| the cells MODFLOW 6 made dry | their numbers | GWT inactivates a cell at HDRY |
| every flow package's flow by entry | SIMVALS (and SIMTOMVR with a mover); each package's entry count and a CRC of its cells, checked in a replay | GWT's sources and sinks: recharge, the coupled rivers, drains and lakes, the soil link, wells, the model's own |
| SWAT+'s day as salt reads it | the `hru_wb` columns, `hru_pw:tmn`, percolation, channel and reservoir volume, outflow, rain and seepage, point sources, withdrawals, fertilizer | the land's and the streams' balances, and the road salt's day rule |
| the soil link's plan by map row | gwflow's rule only | the soil link's split |
| ledger lines and controls | `recharge:sent`, every `user:` line, each user's controls | the salt's point sources, a stated recharge concentration, the users run again |
| the flow, base flow and head series | | outputs and observations |

Once: the run's identity (below), the flow model's memory before the first day, the starting heads, the recharge
and soil maps, the exchange's entries by channel, the users and what they take in, and the SWAT+ files the salt reads
at start-up (routing, HRU data, snow, the stations' precipitation and temperature records). A restarted run carries
on the record it was saved from (the run's identity checked, its flow model's files rebuilt for the record's whole
period and compared), so a long record can be made in pieces.

**Storage.** Each array of a chunk of days (about 256 MB raw) is XORed with the day before, its bytes shuffled (each
value's first bytes together), and compressed with zlib. MODFLOW 6's values are float32 by default, float64 with
`precision = "float64"` (bit for bit), or float32 with `mantissa_bits` fewer than its 23 [A]; SWAT+'s are float32 as
SWAT+ holds them, so exact. Compressors were compared on the Rappahannock's flows: lzma, bz2 and a bit shuffle gave
nothing over zlib on shuffled bytes (within 2%); without the XOR, the face flows and specific discharge took
about 24% more.

Measured on the Virginia study's Rappahannock grid (3 layers, 335,691 active cells, Newton, upstream advection; 8
transient days under recharge varied each day and cell, MODFLOW 6 only), bytes per active cell per simulated day:

| precision | with specific discharge | without | 11 years (4,018 days) |
|---|---|---|---|
| float32 (default) | 20.0 | 12.1 | 27.0 GB / 16.3 GB |
| float32, 16-bit mantissa (within 7.6e-6) | 14.1 | 8.6 | 19.0 GB / 11.6 GB |
| float32, 12-bit mantissa (within 1.2e-4) | 10.5 | 6.5 | 14.1 GB / 8.7 GB |
| float64 | 46 | | 62 GB |

Of float32's 20.0: face flows 6.9, specific discharge 7.9, storage about 3.1, recharge 1.1, saturation 0.7, the
seepage drains 0.2. SWAT+'s part is an estimate, not measured [A]: 16 float32 columns for 33,525 HRUs are 2.1 MB a
day raw, 8.6 GB over 11 years at most, likely 3 to 4.5 GB compressed (MODFLOW 6's float32 flows compressed 1.9
times). MODFLOW 6's own budget file at daily resolution would be about 80 GB (estimated for this grid). A study on this disk records with the
16-bit mantissa, or without specific discharge if it runs without dispersion: 15 to 23 GB. Reading back costs 44 ms a
day on this grid; recording, about 0.5 s a day (the compression).

**Replaying** (`Basin(cfg, replay="<folder>")`, `basinledger run --replay`, `basinledger forward --replay` for
PEST++). The run's MODFLOW 6 simulation is built as the recorded run's was, its flow model's input checked against
the record's, and the flow model left in it, unsolved: each day, once MODFLOW 6 has prepared the day, the flow
model's memory is set from the record (face flows, saturation, specific discharge, storage, dry cells, every
package's SIMVALS), and only the transport solutions are solved. The transport models read those arrays through the
GWF-GWT exchange exactly as in the coupled run. The flow model's output control and observations are left out of
the replay's copy.

Two other ways were looked at and set aside (MODFLOW 6 6.8.2's source):

- **A transport-only simulation reading a budget file (FMI).** The whole daily budget file is the size problem. A
  file of one time step is reused for every step (`advance_bfr`), and the arrays could then be set through the API;
  but reusing sets `IFLOWSUPDATED` to 0, so DSP's coefficients would never be recomputed. Writing the file as the
  replay goes, or streaming it through a pipe, fails: the reader rewinds and peeks ahead. And a file is read as the
  step is prepared, so the transport model would take the day's flows for its dry cells and dispersion where the
  coupled run takes the day before's (below): not the same run.
- **The transport models' FMI arrays set directly.** With the GWF-GWT exchange they are the flow model's own memory,
  so setting the flow model's is the same thing, with the flow model's own checks and package structure kept.

**Bit for bit.** In the coupled run the transport model prepares its day before the flow is solved: it inactivates
dry cells and works out dispersion from the flow model's memory as the day before left it (HDRY heads, saturation,
specific discharge, face flows), and solves with the day's. A replay sets the day's arrays after MODFLOW 6 prepares
the day, so the same holds; before its first day it sets the memory the record holds for the day before (or the
start). So a replay from a float64 record reproduces the coupled run's concentrations, masses and ledger lines bit
for bit, every day (`test_float64_replay_is_the_full_run_bit_for_bit`, 365 days of the synthetic salt example, both
species, every salt, transport and user line). From float32, concentrations agree to 3e-9 relative on the synthetic
example and 1.1e-7 on the Rappahannock grid (the test's tolerance is 1e-6: float32 keeps each flow within 6e-8, and
transport is linear in the flows); masses to 1e-6 of each account's largest day. Every account closes every day in a
replay as in a run: the transport's mass budget is the solved equation's, so float32 flows don't open it.

**Speed.** On the Rappahannock grid (one species, upstream advection, 6 days, the machine loaded by a calibration on
10 of its 14 cores), MODFLOW 6's flow took 4.1 to 4.2 s a day:

| dispersion | full run (MODFLOW 6: flow + transport) | replay |
|---|---|---|
| XT3D (MODFLOW 6's default) | 7.68 s | 3.64 s |
| XT3D off (`add_transport(xt3d=False)`) | 4.52 s | 0.37 s |
| none (`dispersivity=None`) | 4.39 s | 0.24 s |

The full run also steps SWAT+ (33,525 HRUs), which a replay doesn't, so a replay without XT3D is 12 to 18 times
cheaper than the coupled day before SWAT+ is counted. With XT3D, half of MODFLOW 6's time is the transport models'
formulate and budget, and a replay saves only the flow's half. On the synthetic salt example (one HRU, 100 cells, two
species) a day takes 17.7 ms coupled and 6.9 ms replayed.

Transport's convergence on this grid (the study's coupled model timed over 10 days, Nov 22): TVD for nitrate and soluble P
didn't reach `outer_dvclose` 1e-9 kg m-3 in 50 outer iterations on day 2; upstream with XT3D took 12.4 s a day,
without XT3D 5.75, with no DSP 5.08 (flow alone 3.9). `add_transport`'s default `outer_dvclose`, 1e-10 kg m-3, is
1e-7 mg/L, absolute: far below any lab's reporting limit and far tighter, relative to the species, for a trace
species than for chloride. Recommended, not yet changed [A]: a closure relative to the species' concentrations
(about 1e-6 of the largest starting or source concentration, 5e-8 kg m-3 for chloride at 50 mg/L, 5e-9 for nitrate
at 5 mg/L), keeping `rcloserecord` strict, since the mass closure checked every day follows the residual, not the
concentration change; and upstream advection without XT3D on coarse grids, where upstream's numerical dispersion
(about half a 200 m cell) is as large as the physical one [A].

**What may change, and what is refused.** The run's identity is hashed when it is recorded and compared in a
replay (`identity`, `gwf_inputs`):

- the config's flow side: every field but `[salt]`, `[transport]`, `[outputs]`, `[replay]`, the run's folder and
  length, and the builder's arguments (which are checked through the files the builder writes);
- every file it names, by content: the SWAT+ project, maps, users' configs, weather series, zones, pilot points, the
  builder's module;
- the parameters of every model but salt's;
- whether the users run, and the starting heads;
- MODFLOW 6's flow model's inputs: the simulation's options and time, the flow model's name file and every package
  file (with what they read), its solution's settings, any GWF-GWF exchange, comments left out. The transport models,
  their solutions and their exchanges are left out, so they may change.

So a replay may change `[salt.<species>]` (multipliers, sources, road salt, point sources, starting stores), salt's
calibration parameters, `[transport.<species>]` (a stated concentration), and the transport models: porosity,
sorption, dispersivity, XT3D, advection scheme, solver, starting concentrations, weathering rates; and it may add a
species. Refused: anything in the list above that differs; more days than the record holds; `heads` and
`wrap_surface`; a record of another format; a day on which a flow package's entries aren't where they were; a
transport model with dispersion from a record without specific discharge; a transport model correcting for the
flow's imbalance (FMI's FLOW_IMBALANCE_CORRECTION, which reads FLOWJA's diagonal) unless the record has it. A record
is refused for several groundwater models, a NextGen surface, and an exchange other than one RIV package (SFR,
several parts).

**Users.** A user who takes in the aquifer's water through a configured well is run again in a replay, with the
recorded controls: its intake follows the replayed aquifer, so its returned chloride does too. Its water is checked
against the record every day (to 1e-12) and a difference stops the replay: that user's water depends on what it
takes in, so changing the salt would change the flow. A user with forcing who takes in the aquifer's water isn't
built. Every other user's lines are the record's.

**Restarts.** A replay saves its concentrations, salt state, ledger and series (`save`) and restarts from them bit
for bit. A replay can also restart from a full run's saved state, and is then the unbroken full run: the flow
model's memory before the restart day comes from the record, where a restarted full run starts its transport from
the flow model's initial memory (a full run's own restart differs from its unbroken run there, by its first day's
dispersion). Weathering's saturation is in the salt state (Salt), so both carry it.

**What a study config needs.** To record: `[replay] record = "<folder>"` (and `mantissa_bits = 16` or
`specific_discharge = false` on a big grid) on the run whose flow is to be reused, run with `basinledger run`; it
must start from day 0 (a restart carries on its record), with `[salt]` or not (the SWAT+ columns salt reads are
printed and recorded either way). `forward` and `Basin.spin_up` never record. To replay: the same config, with only
the salt side changed, `basinledger forward basin.toml salt_parameters.csv observations.csv out.csv --replay <folder>
--run-folder <worker>` for each PEST++ worker (observations of kind `chloride_mg_l`, `weathering_mg_l` or
`conductance_us_cm` at an output channel, and the recorded flow and heads), or `Basin(cfg, replay=folder)`.

[A], in one place: float32 flows (or fewer mantissa bits) stand in for MODFLOW 6's float64; the SWAT+ files the salt
reads are the run's at its start; a user who doesn't take in the aquifer's water returns what it reported in the
recorded run; MODFLOW 6's flow model is prepared each day in a replay but not solved, and nothing it would write is
kept; SWAT+'s part of a record's size is estimated.

**Not built: a nutrient replay.** Nitrate and soluble P do feed back into SWAT+ (the stream and soil sides of the
nitrate-transport branch), so SWAT+ can't be replayed for them. What could be replayed is MODFLOW 6's flow: SWAT+
live, MODFLOW 6's flow from the record, its transport solved. It would take: (1) the record's flow side as now, plus
what MODFLOW 6 sends SWAT+ each day (each channel's exchange, the soil link by layer, the lakes' leakage), set in
SWAT+ from the record in place of the solved flow's; (2) SWAT+'s water checked against the record every day
(percolation by HRU, channel flows): nutrient parameters can change the water through plant growth (nutrient stress
limits growth, so leaf area, evapotranspiration and soil water), and a replay that changes the water must stop, or
run with SWAT+'s nutrient stress off; (3) the nitrate branch's surface-mass links driven from the replayed flows. It
saves only MODFLOW 6's flow solve (3.9 to 4.2 s a day on the Rappahannock grid; SWAT+ and the transport remain),
against salt's nearly all. Not before the nitrate branch is merged.

## Mass through the human-system routes (designed and built Nov 25, 2026, for the Virginia study's water quality; 0.17.0)

Basin Ledger 0.16.0 moves the water of the basin's existing water-use system: water systems with customers, returns
to the ground over service areas, wells over many cells. None of that water carried mass: `ground:` routes and
`cells = [...]` routes were refused with transport and salt, and a system's water carried no constituent. The
Virginia study's water quality stage (nitrate, soluble P and chloride, calibrated at several sites) needs the mass to
move with the water through every one of these routes before its calibration, so that a town's septic chloride,
a plant's phosphorus and a well field's nitrate are where the water is. Andrew approved the work on Nov 25, 2026.
Three parts: mass on the new routes, constituents through water systems, and no double count with salt's land
sources.

**1. Mass on the new routes (GWT, and salt with it).**

- **A return to the ground carries a concentration into its cells,** on the WEL entry that injects its water (the
  plan the `ground:` bullet above stated). Each transport species gets an auxiliary column on the ground package,
  named after the species, and the package is an AUX source in each transport model's SSM. With `add = true` (the
  default) Basin Ledger writes both: the auxiliary columns in the WEL file it adds, and a `<package> AUX <species>`
  line in each transport model's SSM SOURCES block in the run's copy (`existing.add_ssm_sources`). With `add = false`
  the package must already be an AUX source in every transport model, or the start-up refuses it.
  - Each day, after MODFLOW 6 prepares it, each ground entry's concentration is set on the transport model's copy of
    the package's auxiliary values (as recharge's is): the route's mass for the day over its water. A route over
    several cells gives every cell the same concentration, so its mass splits by the same weights as its water,
    the same day.
  - Chosen over a GWT mass source (SRC) on the cells: on the WEL entry the mass goes only with water MODFLOW 6
    delivers, so a return can't put mass into a cell without its water. A route given mass and no water books it
    (`ground:<name>:no_water`), as recharge does; none is expected, since a route's mass is always a
    concentration times its water.
  - Ledger (kg): `transport:<sp>:ground:<name>:due` = `sent` (the concentration times the water, set on the entries)
    + `no_water`, and `sent` = `received` (the transport model's budget on those entries), every day. The budget
    books the ground entries as coupled, not as a boundary.
- **What a route's mass is.** The sum of what is sent to it: a water system's septic water and leaks (part 2), and
  a user's return. A user's return carries what its component reports for the species
  (`water_user_return_water_<species>__mass_flux`); a user that reports nothing for it (a ReportedUser) carries the
  `[[ground]]` entry's stated `concentration = {<species> = kg m-3}` [A]. A route that a non-reporting user sends
  water to with no stated concentration for a transported species is refused at start-up: Basin Ledger doesn't
  guess what a drainfield or a quarry's recharge basin carries.
- **A well over several cells takes each cell's concentration, weighted by what each cell gave.** Each cell's entry
  is in the transport budget (`well:<name>#i:withdrawn`), and the route's `well:<name>:withdrawn` is their sum; its
  `concentration` is that mass over the water the cells gave (`route:well:<name>:delivered`). A user taking in the
  aquifer's water from such a well is given that concentration for the next day (on a run's first day, the cells'
  concentrations by the route's weights). One-cell wells are as before.
- **Salt.** A salted species' ground returns enter the aquifer the same way (its inflows are salt.py's, and the
  ground entries' concentrations are set beside them). In the salt balance's whole-system closure the returns' mass
  is a source from outside the land, streams and aquifer, as a point source's mass is, since what the water systems
  and users withdrew already leaves it as a sink (wells, places of use): `salt:<sp>:aquifer:from_ground`. The water
  systems' own mass closes in their own accounts (part 2), so together the two close the basin.

**2. Constituents through water systems.** A system's water carries its sources' concentrations, mixed by what each
delivered, to its customers; each class's returns carry a stated concentration, or the mixed water's.

- **Which species.** Every transport species (`[transport.<sp>]`): the SWAT+-carried nitrate and soluble P, the
  salted species, and stage 1's stated ones. Concentrations are kg m-3 (= g L-1; mg/L / 1000).
- **The sources' concentrations, for the day just run** (known once the models' day is done, as the delivery is):
  - `well:<w>`: the aquifer's, as the transport budget took it out (over several cells, weighted by what each
    gave);
  - `swat:place_of_use:<p>`: for a salted species, the mixed concentration of the channel or reservoir it draws
    from (salt.py's water bodies, which take the withdrawal out at it); for nitrate and soluble P, SWAT+'s outflow
    concentration of its channel that day (`channel_water_<x>__mass_concentration`) [A], since what SWAT+'s
    withdrawal takes isn't readable; a reservoir's nitrate and P aren't readable through swatplus-bmi, so a
    reservoir intake needs a stated concentration for them; a stage 1 species takes its stated stream
    concentration;
  - `effluent:<system>`: the other system's plant outflow concentration that day;
  - any other (`offsite:`, an import): a stated concentration, required;
  - `concentration = {<species> = kg m-3}` on a source states it [A], over any of the above.
  Mixed: the system's concentration = the sources' mass over the water they delivered (`utility:<s>:<sp>:
  concentration`). Everything the system delivers, loses and leaks that day carries it.
- **The returns, by class:**
  - **sewered:** to the plant. The plant's outflow carries its effluent concentration, `effluent = {<species> =
    value}` on the system [A]: a constant, 12 monthly values, or a CSV of dated values (`date`, then a column per
    species; each row held until the next, a day before the first refused), so a plant's DMR quality by month
    reproduces. `"sewer"` passes the sewer's mass through (a plant that removes nothing: chloride, a water plant's
    backwash). The sewer's mass is the tap water's in the classes' sewered water, plus the users' returns at what
    they report (the tap water's where they report nothing [A]); inflow and infiltration, water from outside,
    carry none [A]. What the plant changes is booked (`plant:added`,
    negative where it removes). The outflow's mass leaves the next day with its water, at the effluent
    concentration of the day it leaves: to the outfall, and to the systems reusing it.
  - **septic:** to the ground through its `septic_to` route, at `septic_effluent = {<species> = value}` (on the
    system, or on a class, which wins) [A]: what reaches the water table below the drainfield, after the tank and the
    soil under it (most soluble P sorbs there; nitrate and chloride mostly pass). The same forms, and `"tap"` for
    the mixed water's. What use adds (`septic_added:<class>`) is booked. It leaves the next day with its water.
  - **leaks:** the mixed source water, to the ground the next day.
  - **consumed** (outdoor use) and unaccounted losses: **the mass is dropped, booked as leaving** (`consumed:kg`,
    `unaccounted:kg`). Chosen over leaving it on the land because the water itself leaves the model: Basin Ledger
    books outdoor use as consumed and puts no water on any HRU, so the mass would land where the model has no
    water carrying it; and swatplus-bmi has no way to add nitrate or P to an HRU's soil (the gap irrigation from
    wells has). For chloride, the mass that leaves is the tap water's on the outdoor share (about 10-15% of a
    town's production at a few tens of mg/L), small beside septic and road salt; a study that wants it on the land
    adds it as a land source.
  - Required: `effluent` for each species if the plant can have an outflow (a sewered class, a user returning to
    it, or inflow and infiltration), and `septic_effluent` if a class has septic water. Missing ones are refused at
    start-up, by system and species; a species that isn't transported is refused (a typo).
- **The users on a system.** Each takes the system's concentration of the day before as its intake
  (`water_user_intake_water_<sp>__mass_concentration`), if its component takes one, the one-day lag the wells'
  intakes have.
- **The outfall.** A plant outflow to `swat:point_source:<p>` carries its mass into the stream:
  - for a salted species, salt.py's water bodies take it (the system reports its mass for the point source, as a
    direct user does; the rest of the point source's water at the stated `point_sources` concentration);
  - for nitrate and soluble P, Basin Ledger sets SWAT+'s point source load before SWAT+'s day
    (`point_source_water_<x>__mass_flow_rate`): the systems' outfall mass, plus the direct users' reported mass,
    plus `[transport.<sp>] point_sources = {<p> = kg m-3}` [A] times the rest of the point source's water (the
    salt balance's rule). Water no one reports and no concentration is stated for carries none, as before 0.17.0,
    and is booked (`transport:<sp>:point_source:<p>:unstated_m3`) with a warning at start-up. Booked:
    `transport:<sp>:point_source:<p>:sent` (kg). This is the "user's return loads to SWAT+'s point sources" stage
    2 left out.
  - to `outfall:` or `offsite:`: booked as leaving.
- **Timing,** as the water's. Day d's sources' mass is known after the models' day; the classes' returns are made
  then and held overnight: septic and leaks with their mass, the sewer with its mass. Day d+1's start releases the
  plant's outflow at that day's effluent concentration, and the ground returns, before SWAT+ (the point sources'
  loads) and MODFLOW 6 (the ground entries' concentrations). A restart carries the held mass with the water.
- **Ledger** (kg, `utility:<s>:<sp>:`), every day: `source:<route>:kg` and `produced:kg`; `concentration`;
  `delivered:kg`, `losses:kg`, `leaks:kg`, `unaccounted:kg`; by class `delivered:<c>:kg`, `sewered:<c>:kg`,
  `septic:<c>:kg`, `septic_tap:<c>:kg`, `septic_added:<c>`, `consumed:<c>:kg`; `users_delivered:kg`,
  `users_return:kg`; `sewer:kg` and
  `held_to_ground:kg` (held); at the release, `plant_outflow:kg`, `plant:added`, `outfall:kg`, `reuse:kg`,
  `reuse:<n>:kg` and `to_ground:<route>:kg`. Pairs: produced = the sources'; produced = delivered + losses;
  delivered = the classes' + the users'; each class's delivered = sewered + the tap water in its septic water +
  consumed, and its septic = that tap water + what use added;
  losses = leaks + unaccounted; plant outflow = outfall + reuse. And the system's closure: what it holds overnight
  changes by the sources' mass, the users' returns and what use and the plant add, less what it delivers to users,
  drops (consumed, unaccounted) and releases (to the ground, the plant's outflow), `closure_share` within 1e-12.
  Across the links: the ground routes' due = what the systems and users sent them, and each point source's mass =
  its systems', users' and stated parts.

**3. No double count with salt's land sources.** Salt's land sources can carry septic chloride by HRU (the study's
septic estimates by HRU, group `septic`). Where a town's septic water now reaches the aquifer through `ground:`
routes with its own chloride, the same homes' chloride would be counted twice. **The rule: septic water that reaches
the ground through a route carries its own mass; a land source in group `septic` carries only the septic systems no
route carries.** Basin Ledger can't see which homes a land source covers, so a salted species with both a land
source in group `septic` and septic water on a route (a system's septic class, or a user's return to the ground) is
refused unless that land source says `outside_routes = true`: the study states it built the source without the
routes' customers. The routed septic concentrations of a salted species scale by the same `septic` multiplier, and
the plants' stated effluent concentrations by `point_sources` (`"tap"` and `"sewer"` aren't scaled), so a calibration's
two multipliers scale every septic and sewer plant source alike. For nitrate and soluble P the same rule holds against
SWAT+'s own septic systems: a project with septic in `landuse.lum` and septic water on a route carrying them is
refused.

**4. Reported irrigation applied on the land model's fields (designed Nov 25, 2026, approved by Andrew the same day;
built Nov 26 on swatplus-bmi 0.15.0, its request 9; 0.18.0).** The Virginia study's irrigators are golf courses, a turf farm, a nursery and a
school's fields. They report monthly withdrawals to DEQ, from channels, ponds and wells, and today the study counts
that water as consumed at the source. It should land on the HRUs it irrigates, so that SWAT+ makes its ET,
percolation and runoff, and carries the nutrients and salt in it. Unreported farm irrigation, from USGS's monthly
HUC-12 estimates, takes the same path.

- **A route onto the land: `land:<name>`** (`[[land]] name, hrus = [{hru, weight}, ...]` or a CSV of `hru, weight`,
  weights normalized; `runoff = 0.0`, the share of the water that runs off the field [A]). An irrigator is a user
  like any other, a ReportedUser from its DEQ months:
  - its withdrawal is any withdrawal route: `well:` (over several cells for HUC-12 groups of farm wells),
    `swat:place_of_use:` on a channel or a reservoir (a pond SWAT+ models as one), or a water system's tap;
  - its return goes to `land:<name>`, and its consumption is the application's loss (spray evaporation and drift)
    [A].

  Nothing in the engine is irrigation-specific: the route is a place water goes, as `ground:` is. A system's
  outdoor class can later send its water there too (`outdoor_to = "land:<town>"`), which would put outdoor use's
  water and mass on the land, now dropped (part 2).
- **Applied the day after it is withdrawn, as a return is held.** The water a place of use delivers is known only
  after SWAT+'s day, and a well's water's concentration only after MODFLOW 6's, but SWAT+ applies irrigation inside
  its HRUs' day, before either is known. So day d's delivered water, and its mass, go on the HRUs on day d+1, split
  by weight. The route then closes exactly for water and mass. For monthly reported volumes the lag is immaterial.
  Andrew chose same-day irrigation for SWAT+'s own decision tables (Irrigation from surface water above), but there
  SWAT+ takes the water itself; here Basin Ledger moves it, and a same-day split couldn't close.
- **What SWAT+ is asked to do (request 9).** Add each HRU's water to its day's irrigation: to the soil (`irrig(j)%
  applied`) and to surface runoff (`irrig(j)%runoff`) by the route's runoff share, after the HRU's management
  operations, so a decision table's own irrigation that day adds rather than overwrites; and add its nitrate and
  soluble P to the top soil layer the same day. SWAT+ 62.0.1 has no way in: its water allocation's delivery to HRUs
  (`wallo_pou_deliv`, case `irr`) reads an HRU index it never sets, an `irrigate` action applies a constant depth
  from its decision table, and upland irrigation's nitrate is never added to the soil (only a paddy's, in
  wetland_control.f90). swatplus-bmi 0.14.0's settable variables are point sources, places of use, the exchanges,
  the soil link's water and mass by layer, and the weather; the soil link's variables are groundwater from below,
  not water on the field, so they aren't borrowed.
- **Mass.** The whole withdrawal's mass goes on the land: its consumption is evaporation and carries none. It is the
  user's reported return mass if its component reports one (scaled by the share delivered); else the `[[land]]`
  entry's stated `concentration` [A]; else its source's water's as part 2 works it out (a well's cells by what they
  gave, a channel's or a reservoir's for salt, SWAT+'s channel outflow for nitrate and P, a water system's mix),
  times the water delivered. A source with none of these (a reservoir intake's nitrate and P) is refused at
  start-up. The mass is fixed at the day's end (`user:<u>:land:<sp>:kg_held`, its `concentration` beside it) and
  goes on with the water the next day. It all goes into the top layer, the runoff's share too [A]: swatplus-bmi
  0.15.0 puts set mass there and refuses mass without water into the soil, so `runoff` stays below 1. Nitrate and soluble P go through
  request 9's mass variables. A salted species goes into the HRUs' soil stores in salt.py (`salt:<sp>:land:
  irrigation`), an input from outside the land like a land source, while its withdrawal leaves the streams or the
  aquifer as a sink, so the whole-system closure holds as with the ground returns.
- **Booked:** `route:land:<name>:requested` (what was set) = `:delivered` (swatplus-bmi's read-back of what went
  in, an HRU's shared among the routes on it by what each set) = `:into_soil` + `:runoff`; the held water and mass
  overnight (`user:<u>:return_held`, `land_held` in a restart); `transport:<sp>:land:<name>:sent` = `:received`
  (SWAT+'s read-back) for nitrate and P; `salt:<sp>:land:source:irrigation` for a salted species. In the whole-basin
  balance the water into the soil is in SWAT+'s `irr`, and the runoff share, which SWAT+ prints in its runoff
  (`surq_cont`) rather than `irr`, is a land input of its own (`irrigation_runoff_from_routes`).
- **What 0.15.0 does with it, documented there and not changed:** on a day with rain, SWAT+ 62.0.1 takes an
  irrigation's runoff out of the rain's infiltration (`precip_eff - surfq`), so the runoff adds water to the HRU only
  by what exceeds that; a sprinkler's `runoff` is best left at 0.
- **Refused:** a route's HRU that SWAT+ doesn't have; a wetland or paddy HRU (`surf_stor` in `hru-data.hru`: SWAT+
  puts its irrigation into the ponded water); weights not positive, `runoff` outside [0, 1); the same HRU irrigated
  by `[irrigation]` (wells pumped for SWAT+'s decision tables) and a `land:` route at once, since the two would mix
  two rules for one field's water; a land user drawing from anything but a well, an intake or a water system; and
  a swatplus-bmi before 0.15.0.
- **In a replay** the water is the record's; the mass is worked out again from the replay's own concentrations and
  goes to salt's land stores, and a `[[land]]` entry's stated concentration may change. A user taking in the aquifer's
  or a system's water and putting it on the land is refused there.
- **Tested** (tests/test_land.py), on Osu_1hru made upland:
  - a golf course on an intake (its delivery short on some days) and a farm on a well over two cells;
  - each day's delivered water lands the next day, split into the soil and the runoff, all of it taken;
  - the mass is what the stream and the aquifer gave up, reaching salt's land the next day, with every closure
    held;
  - nitrate and soluble P reach the top layer as set;
  - a restart and a float64 replay are bit for bit.

**5. In a replay (built Nov 25, 2026; the study's salt calibration runs by replay).**
- **The record** holds each day's route lines (`route:<r>:requested`, `:delivered`, `:shortfall`) beside the users',
  and which routes the run used. The water systems are accounting, so the replay runs them again from those lines:
  their water comes out the same, and their mass is worked out from the replay's own concentrations.
- **The replay's MODFLOW 6 copy** gets the same ground packages, wells over several cells, auxiliaries and SSM
  sources as the run's (`driver.add_route_wells`, shared by both).
- **Run again:** the users taking in a well's water (over several cells too) or a water system's. A user on a system
  returning its water anywhere but the sewer (held overnight) is refused there.
- **What only moves mass may change in a replay:** a return to the ground's stated concentration, a system's
  `effluent`, `septic_effluent` and its sources' `concentration`, as well as salt's settings and multipliers. These
  are left out of the record's identity. The routes' septic and plant masses follow the septic and point_sources
  multipliers.
- **A restart** carries the systems' held water and mass.
- **A record made before 0.17.0** holds no route lines, so it is refused for a config with these routes.
- **Tested** (tests/test_replay_routes.py):
  - a float64 record replays the town's water and mass, the ground routes, the salt and the transport bit for bit;
  - a float32 record replays them to its rounding;
  - a changed septic concentration and multiplier scale the routed septic mass, and every account closes;
  - a replay restarted is the unbroken replay.

**Not built, refused at start-up:** a place of use's concentration
under NextGen, or a reservoir intake's nitrate or P, with none stated. A
system's concentrations are mixed daily; the water in its mains is not a store.

**Tests** (tests/test_ground_mass.py, tests/test_chain.py, tests/test_systems_mass.py): a drainfield's chloride over
three cells, each cell's share of the mass its weight, the transport model's mass closing every day; a well over two
cells withdrawing each cell's concentration by what it gave; the salt balance's whole system closing with ground
returns; a system mixing an intake and a well, its classes' septic and sewered returns at their stated
concentrations, its leaks at the mix, the plant's DMR-like monthly effluent reaching the point source (chloride
through salt, nitrate through SWAT+'s load), every pair and the closure every day; the double count refused, and
allowed with `outside_routes = true`.

## Calibration

Using Basin Ledger depends on calibrated models. Calibration is part of every example, not something a study adds
later. `calibration.py` has:

- the record and its status;
- parameter and observation tables;
- simulated values at observations;
- the statistics;
- a forward model that runs the driver (`python -m basinledger forward`);
- parameters written into a run's copy of the models' files (built Nov 11, 2026);
- PEST++ IES files (built Nov 11, 2026).

Pilot points and USGS groundwater level downloads are still to build. This section is the whole design.

**Parameters (built Nov 11, 2026).** The engine exposes each model's adjustable parameters as one named list, with
bounds (`parameters.csv`: name, model, target, value, lower, upper, transform). They are written into the run's copy
of each model's files before it starts (`apply_parameters`, through the driver's `prepare` step): the example's own
set (`[calibration] parameters`, the calibrated values every run uses), with a calibration's values in place of any
it names.

- **SWAT+:** through its own `calibration.cal`, read once before the first day (swatplus-bmi 0.9.0). A target is
  `<name>:<change>`, with optional objects, soil layers and an HRU's soil group, texture or land use, e.g.
  `cn2:pctchg`, `awc:absval:layers=1-2`, `cn2:pctchg:hsg=D`. The rows go after the project's own. swatplus-bmi
  refuses a change SWAT+ 62.0.1 reads but wouldn't make as asked (its audit of every parameter in the editor's
  `cal_parms.cal`: 91 made, 50 refused, 22 with no code, 111 not checked), so a calibration can't move a parameter
  that moves nothing. A run with changed parameters is an ordinary run.
- **MODFLOW 6:** each value multiplies a package value as built, everywhere, in one layer or in one zone
  (`[calibration] zones`: a CSV of cells and zones, such as a study's hydrostratigraphic units): `npf.k`, `npf.k22`,
  `npf.k33`, `sto.sy`, `sto.ss`, and the columns `riv.cond`, `drn.cond`, `ghb.cond`, `sfr.rhk` and `sfr.man`. Only
  the packages changed are written again, with flopy, at full precision, so a multiplier of one gives the plain run
  bit for bit (tested). Recharge isn't one: Basin Ledger sets it each day from the surface model.
- **MODFLOW 6 by pilot point** (built Nov 15, 2026): a parameter targeting `<package>.<array>:pp=<set>` is one pilot
  point, its location a row of `[calibration] pilot_points`' file (`name, set, layer, x, y`, in the grid's world
  coordinates). A set's values are interpolated over each layer they're on by ordinary kriging, an exponential
  covariance of the given `range` with an optional `nugget` (`calibration.kriging_weights`), in log10 where the
  set's parameters are log-transformed, and the field multiplies the array there; a layer without points is
  unchanged. Written in numpy rather than with pyEMU: the system is small and solved once per layer, and Basin
  Ledger writes its own PEST++ files, where each point is an ordinary parameter row. Tested: the weights honour
  each point and add to 1, far from all points they give the points' declustered mean; on the small model each
  point's value lands at its cell and the field is the kriged one.
- **Water users:** parameters in each component's config file, e.g. `CoolingWaterBmi`'s cycles of concentration.
  These are set from the applicant's data, not calibrated; a user parameter in the list is refused.
- A parameter for a model the example doesn't run is refused.
- A restart writes the parameters into the new run's files again (MODFLOW 6's; SWAT+'s saved state already has
  them), and matches the unbroken run bit for bit (tested).

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
    `python -m basinledger forward basin.toml parameters.csv observations.csv simulated.csv --run-folder runs`
    (`calibration.run_forward_model`). It applies the parameter set, runs every model in the worker's own folder,
    and writes the simulated value at every observation for PEST++ to read.
  - `calibration.write_pest_files(basin.toml, parameters.csv, observations.csv, folder)` writes the control file
    (`basin.pst`), a template of the parameter table and the instructions for the simulated values (built Nov 11,
    2026). Parameters start at their table values within their bounds, log-transformed where the table says;
    observations keep their weights, grouped by kind and site; `ies_num_reals` sets the ensemble. The files are
    checked against each other in the tests; PEST++ itself isn't run there.
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
- **MODFLOW 6:** its state is its heads, plus the storages of advanced packages (MAW, UZF; SFR with its storage
  option, which the exchange refuses) if used. SFR without storage holds no water from day to day.
  At the end of the warm-up the driver saves heads (`X`), and each lake's stage (LAK `XNEWPAK`) and multi-aquifer
  well's head and wellbore storage head (MAW `HEAD`, `XSTO`) beside them, set back before the first day (built Nov
  15, 2026: a model's own pond and an idle farm well restart to the last digit, 1.47 m off without; UZF and CSUB
  still start from their files, and the start-up says so). Each scenario starts a new MODFLOW 6 simulation covering
  only the days left: the example's builder is called for them, or a copied model is cut to them
  (`remaining_simulation`, `timesteps.remaining_days`, built Nov 14, 2026): its steps rewritten to days first, the
  periods already over left out, the one the restart falls in shortened, every input file's period blocks renumbered
  (the block in effect that day becoming period 1), and its time series and time-array series moved back by the
  days run. The files are changed as text, so nothing else in them is rewritten or rounded. The saved heads are then
  put straight into memory after initialize, not into the IC file, which would round them. Restarted runs match the
  unbroken run bit for bit (built and tested, a copied model of monthly periods included). Advanced packages (MAW,
  UZF, LAK, SFR with storage) will need their own state saved.
- **Water users:** a component's state, if any, is saved through `get_state`/`set_state` where the component offers
  them (CSDMS bmi PR #112). `CoolingWaterBmi` keeps no state between hours.
- **Weather:** SWAT+'s gage records are part of its state, so a case keeps the warm-up's weather up to its restart
  and writes its own (`[weather]`) from the first day SWAT+ hasn't read.
- **Driver:** the ledger, the series, the rules' memory and the exchange due are saved with the rest
  (`Basin.save`, `Basin(..., restart=...)`, `--save-at`/`--restart` on the command line). The users are stepped
  forward to the restart day, which reproduces a deterministic component exactly.

Scenarios then run in parallel, each in its own process. `branch()` (fork) isn't used by the driver yet, because
libmf6 lives in the driver's process. A branch would need MODFLOW 6 restarted from heads in the child process.

## Scenarios of water users, and what each one does to the basin (designed and built Nov 25, 2026, for the Virginia study)

Basin Ledger exists to model and investigate the basin-scale effects of individual and multiple water users. The
Virginia study (the Rappahannock above Fredericksburg) has the basin's existing water use as 25 systems and about 80
users. Its scenarios add the proposed data centers to that system: Culpeper's technology zone, Remington and
Wilderness Crossing, each in phases. It reports each campus's effect and their combined effect. The warm-up scenarios
above change the weather. These change users. Andrew asked for them on Nov 25, 2026.

**A scenario is a set of changes to users, on top of the existing system.** The baseline is the example as it is:
its users, systems and rules. A scenario says what it adds and changes. A set of scenarios is one TOML file beside
the example (`scenarios.py`, `basinledger scenarios <file>`):

    name = "data-centers"
    example = "basin.toml"                 # the existing system, the baseline
    # folder = "../runs/data-centers"      # default: the example's run folder, then the set's name

    [warm_up]
    days = 5479                            # the cases carry on from the end of these days
    users = true                           # the existing users run in the warm-up (default)

    [run]
    days = 0                               # the run's days in all, warm-up included (0: the example's)
    workers = 4                            # cases at once, each in its own process
    nice = 10                              # each worker's niceness (default 10; POSIX, not Windows)
    checkpoint_days = 365                  # a case saves its state this often; a stopped case resumes there (0: never)

    [compare]
    channels = ["cha045", "cha112"]        # flow out of these channels
    heads = [{name = "tz_well", cell = [0, 120, 88]}]
    # from = "2027-01-01"                  # the period compared (default: the warm-up's end to the run's end)
    # to = "2045-12-31"
    low_flow = {days = 7, return_period = 10, year_start = "04-01"}   # the default: 7Q10, climatic years from April 1
    exceedance = [0.10, 0.50, 0.90, 0.95, 0.99]                       # flow duration points (the default)
    units = "SI"                           # or "US": cfs, ft, million gallons

    [[scenario]]
    name = "culpeper"
    description = "The technology zone's campuses, on the town's system"

      [[scenario.users]]                   # added: a [[users]] table as in the example
      name = "culpeper_tz1"
      component = "waterledger_models.datacenter.bmi:CoolingWaterBmi"
      config = "users/culpeper_tz1.yaml"
      withdrawal = "utility:culpeper"      # a customer of an existing system
      return = "utility:culpeper"
      forcing = {source = "swatplus", station = "sta012"}
      phases = [{from = "2027-01-01", mw = 60}, {from = "2030-01-01", mw = 300}]

      [[scenario.scale]]                   # an existing user or system, scaled
      utility = "culpeper"                 # or user = "<name>"
      multiplier = 1.15                    # from the scenario's start; or from = "<date>"; or phases = [...]
      # customers = ["homes"]              # a system's classes only (default: all of them)

    [[scenario]]
    name = "remington"
    remove = ["user:old_quarry"]           # an existing user or system, scaled to zero from the scenario's start

      [[scenario.users]]                   # a user with its own source
      name = "remington_dc"
      withdrawal = "well:rem1"
      ...

      [[scenario.wells]]                   # its well: a [[wells]] table, added to every case's model
      name = "rem1"
      cell = [0, 64, 40]
      add = true

      [[scenario.utilities]]               # a system's sources changed: the keys given replace the example's
      name = "fauquier"
      sources = [{route = "swat:place_of_use:rem_intake", priority = 1, capacity = 0.05}, {route = "well:rem1", priority = 2}]

    [[combination]]                        # the scenarios together, compared with each alone
    name = "all"
    scenarios = ["culpeper", "remington", "wilderness"]

What a scenario may hold:

- `users`: users added, `[[users]]` tables (their paths relative to the set's file). A user is a customer of a
  system (`utility:<name>`) or has its own source: a well, an intake on a channel or reservoir
  (`swat:place_of_use:<name>`), or reuse (a system of its own whose source is `effluent:<system>`).
- `wells`, `ground`, `places_of_use`, `point_sources`: the models' new parts, as in the example. A restart carries on
  the same models, so these are added to every case and to the warm-up, at zero where a case doesn't use them. The
  baseline has them too, and nothing reaches them.
- `utilities`: a name the example doesn't have adds a system. A name it has changes that system: the keys given
  replace its own (`sources`, `shares`, `losses`, `leaks`, `customers`, `outfall`). New `sources` drop a dated
  `shares` table unless the scenario gives one. A change of shares or priority is a change of `sources`.
- `scale`: `{user = ... or utility = ..., multiplier}`, from the scenario's start, from `from`, or by `phases`. A
  user's water, mass and power scale (below). A system's customers' demand scales, all classes or the ones named.
  Its users scale only when they are named themselves.
- `remove`: `"user:<name>"` or `"utility:<name>"`, scaled to zero from the scenario's start. A user is kept and
  stepped, not deleted, so whatever it held overnight at the warm-up's end still leaves the next day, and the ledger
  closes. A system removed keeps its sources and plant, which carry yesterday's water out and then sit at zero. Its
  users must be removed with it, or the start refuses.
- `rules`: `[[rules]]` tables added (a curtailment on a new user).

A combination's changes are its scenarios' together. Two scenarios that change the same user or system in the
combination are refused, as is a name added twice.

**Phased growth.** A user's size follows a dated schedule:

    phases = [{from = "2027-01-01", mw = 60}, {from = "2030-01-01", mw = 300, config = {design = "hybrid"}}]

Each phase starts on its date and lasts until the next. Before the first, the user is off. A phase may give:

- `mw`: the IT capacity, CoolingWaterBmi's `it_capacity_mw`;
- `config`: any keys of the component's config for the phase (`mw = 60` is `config = {it_capacity_mw = 60}`);
- `multiplier`: a factor on everything the user takes, returns and uses (default 1).

A phase without `mw` or `config` keeps the config before it. CoolingWaterBmi takes its capacity only at initialize,
and its interface isn't Basin Ledger's to change. So a phased user is one component per distinct config, each built
from the user's config file with the phase's keys put in (written to the run's `users/<name>/`, in the file's own
format). They are all stepped together from the run's start, with the same weather, and each is dropped once its last
phase has passed. Every phase is therefore exactly the component at that size, whatever state or calendar it keeps.
The active phase gets the limits and dry mode from the rules; the others run without them, and no one reads them.
The user is routed as one user, and its accounts are one user's. A schedule of multipliers alone uses one component
built from the file as given. A ReportedUser grows that way. A multiplier is put on each step, before the day's
sum, so a power of two scales the day to the bit, and any other to the sum's rounding (1e-15).

The multiplier scales the outputs whose units are a volume, a volume flux, a mass flux or a power. A temperature isn't
scaled, and neither is a running total of time (CoolingWaterBmi's hours with heat not removed), which counts the
active phase's hours. A withdrawal limit reaches the component divided by the multiplier, so the limit binds on what
the user takes. The user's balance holds step by step, since every term of it is multiplied together. A
`scale` on a user with phases of its own multiplies the two schedules.

**What is compared.** Every case is the baseline with the changes, run from one warm-up through the same weather.
Over the compared period, each case against the baseline:

- **Flow** at each channel in `[compare]`: the daily change, and the change in these statistics. The mean. The
  flow exceeded 10, 50, 90, 95 and 99% of days, from the flow duration curve. The lowest 7-day mean flow in the
  period. The mean of each climatic year's lowest 7-day mean, and the 7Q10: the lowest 7-day mean expected once in
  10 years, by log-Pearson type III on the complete climatic years' minima (USGS practice; the frequency factor by
  Wilson and Hilferty's approximation, the skew the sample's), with the number of years beside it. A 7Q10 from fewer
  than 10 years extrapolates, and the table says how many it had; below 3 it is left out. If a year's minimum is zero,
  log-Pearson can't be fitted, so the empirical quantile (Weibull plotting positions) is given instead, marked. Each
  case's flow duration curve is a table of its own.
- **Heads** at each cell in `[compare]`: the daily change; the change in the mean; the largest drawdown and its
  date; the change on the last day.
- **Capture** of the case's extra pumping (`analysis.capture`, every well's pumping, a system's included): from
  storage, streams, the boundary, recharge, the soil and returns to the ground, and the residual, which is zero.
- **Use:** each user the case adds or changes, its withdrawal, delivery, shortfall, consumption and return, and its
  energy, as totals and mean rates. Each system's production by source, its unmet demand and its outfall: how a new
  customer's water is found and where its wastewater goes.

**Do effects add up?** Each scenario runs alone, and each combination runs its scenarios together. For any
statistic m, a case's effect is E(case) = m(case) - m(baseline). A combination's interaction is

    I = E(together) - (E(first alone) + E(second alone) + ...)

It is zero when the users don't affect each other: their cones of depression don't overlap, a shared source doesn't
run short, a stream doesn't go dry, a rule doesn't fire. The table gives each scenario's effect, their sum, the
effect together, and I, both as a value and as a share of the effect together. For the daily flow and heads it also
gives the share of the combination's day-by-day change that isn't the sum of the parts (the sum of |I| over the sum of
|E|) and the largest |I| on a day. The users' own water adds exactly when nothing feeds back on them (no rule, and no
limited source), and the test holds it to the last bit.

**The runner.** `basinledger scenarios data-centers.toml [--workers 4] [--cases culpeper,all] [--fresh]
[--tables]` runs the
warm-up once, then the baseline, each scenario and each combination, in parallel processes (spawned, one case per
process, since libmf6 and SWAT+ hold one model a process), each at the worker's niceness. Then it writes the tables.
It is restartable. The warm-up and each case record a fingerprint when done: the example's file, its parameter file,
the set's model parts and comparison points, and the case's own changes and the contents of its users' config files.
A case whose fingerprint matches is not run again. A case stopped part way resumes from its last checkpoint. Adding a
scenario runs only the new scenario. Adding a well or a head to compare changes every case's models, so everything
runs again, and the runner says why. `--fresh` runs everything. A model file the example reads that changes under
the same names isn't seen: use `--fresh`. `--tables` writes the tables again from the cases already run. One case's
failure doesn't stop the others. It is listed in `status.json`, and the runner exits with an error after the
tables.

The set's folder holds `warm_up/` (its run and `state/`), `cases/<case>/` (each run's outputs), `checkpoints/` and
`tables/`:

- `flow.csv`, `flow_duration.csv`, `heads.csv`, `capture.csv`, `use.csv`, `systems.csv`, `interaction.csv`: one row
  per case and point (or user, or system), the baseline's value, the case's and the change, in the units named in
  each column's header.
- `daily_flow.csv`, `daily_heads.csv`: each day, the baseline and each case's change.
- `tables.md`: the main tables, rounded, for a paper; a screening example says so at the top.

**Operating rules for systems and reservoirs (designed and built Nov 25, 2026; Andrew approved adding them that
day).** The study's scenarios also need the basin's operations:

- **Pumped storage.** Motts Run Reservoir (res28) and Lake Mooney (res26) are filled from the river. Hunting Run
  (res27) pumps from the Rapidan at high flow and releases back at low flow, so the intake downstream can keep
  pumping. The refill isn't reported, so the study now pumps each month's draw from the river [A].
- **Drought restrictions** on the systems.

These are rules (`rules.py`) acting on two new things: transfers, and a system's customer classes.

**A transfer** moves water from a withdrawal route to a return route, at the rate the rules set each day (0 when no
rule acts). It is accounting, as a system is:

    [swatplus]
    places_of_use = [{name = "hunting_pump", channel = "cha40"}, {name = "hunting_out", reservoir = "res27"}]
    point_sources = [{name = "hunting_in", reservoir = "res27"}, {name = "hunting_back", channel = "cha40"}]

    [[transfers]]
    name = "hunting_fill"
    from = "swat:place_of_use:hunting_pump"   # a SWAT+ intake on a channel or a reservoir
    to = "swat:point_source:hunting_in"       # a point source into a reservoir or a channel (or outfall:, offsite:)
    capacity = 0.6                            # m3 s-1, the pump's

A point source may now go into a reservoir (`{name, reservoir}`): SWAT+'s recall object is routed to it, and the
reservoir takes it with its other inflow. A transfer's withdrawal is known only after SWAT+'s day, and a point
source is set before it. So the water delivered today goes into `to` tomorrow, as an intake user's return does. It
is held overnight and carried by a restart. A transfer on a shared intake gets its part of what the intake delivered.
The ledger books, under `transfer:<name>:`:

- `requested` = `delivered` + `shortfall`;
- `delivered` = `sent` + `held_gain` (what went in today is yesterday's delivery).

Each route's own pairs hold as before.

**A transfer's mass (built Nov 25, 2026, the same night, for the study's water quality runs).** A transfer's water
carries the intake's concentration of each species salt or SWAT+ carries. The concentrations come from the systems'
sources (Mass through the human-system routes, part 2):

- for chloride and weathering's ions, salt's channel or reservoir, which takes the withdrawal out at that
  concentration;
- for nitrate and soluble P, SWAT+'s outflow concentration of the channel [A];
- a reservoir's nitrate and P aren't readable through swatplus-bmi, so a release from one states them
  (`concentration = {<species> = kg m-3}` on the transfer [A], which may also state any other species). The start-up
  refuses a transfer without them.

The mass is held overnight with the water and goes in with it the next day:

- into salt's water body at the point source, reported as a system's outfall is;
- or as SWAT+'s point source load (`point_source_water_<x>__mass_flow_rate`), with the systems' and the users';
- or, to `outfall:` or `offsite:`, it leaves.

A transfer's water is its part of what SWAT+ withdrew, not capped at what it asked, so salt's water body loses
exactly the mass the transfer carries. Booked under `transfer:<name>:<species>:`: `concentration`, and
`delivered:kg` = `sent:kg` + `held_gain:kg`. A restart carries the held mass.

A replay record holds the transfers' water lines, what the rules did. A replay enters them again and works out
their mass from the replay's concentrations. The record also holds each day's restriction cuts, and the replay puts
them on the systems' classes before it plans them again. A record made before this is refused for a config with
transfers or restrictions: make a new one. A transfer's stated concentration is part of the flow's identity, so a
replay can't change it.

**The rules,** each read at the pause from yesterday's flows and the reservoirs' volumes:

    [[rules]]
    kind = "pumped_storage"
    transfer = "hunting_fill"
    channel = "cha40"                 # the river pumped
    above = 2.0                       # m3 s-1: only the flow above this is pumped
    reservoir = "res27"
    below = 1.0                       # fill up to this share of the reservoir's normal pool
    # pool = 1.2e6                    # m3; default SWAT+'s principal spillway volume (hydrology.res vol_ps, ha-m)

    [[rules]]
    kind = "release"
    transfer = "hunting_release"      # from the reservoir's intake to a point source on the channel
    channel = "cha52"                 # the channel kept up (the downstream intake's)
    below = 1.0                       # m3 s-1: release up to its deficit below this
    reservoir = "res27"
    keep = 0.3                        # never below this share of the pool (optional)

    [[rules]]
    kind = "restriction"
    utility = "culpeper"
    stages = [{when = {channel = "cha12", below = 1.5}, cut = {outdoor = 0.5}},
              {when = {reservoir = "res30", below = 0.6}, cut = {outdoor = 1.0, indoor = 0.1}}]

- **`pumped_storage`** sets the rate when the river is above `above` and the reservoir below `below` x pool. The
  rate is the smallest of the capacity, the river's flow above the threshold, and the room left to the target
  volume over a day.
- **`release`** sets the rate when the channel is below `below`. The rate is the smallest of the capacity, the
  channel's deficit, and what the reservoir holds above `keep` x pool over a day.

  The released water reaches the channel the day after it leaves the reservoir, so the release answers the day
  before yesterday's low flow.
- **`restriction`** cuts a system's customer classes by share. Each day the last stage whose trigger held yesterday
  applies: a channel's flow below a threshold, or a reservoir below a share of its pool. A class asks for its demand
  times (1 - cut). The ledger books it as a restriction, not a shortfall, under `utility:<name>:`:
  - `wanted:<class>` = `demand:<class>` + `restricted:<class>`;
  - `restricted`, the day's total.

  `unmet` stays what the sources didn't deliver. A user on the system isn't cut: its withdrawal is its component's,
  and curtailment and dry mode act on it. The stages have no hysteresis yet: a stage ends the day after its trigger
  stops.

When several rules set one transfer's rate, the smallest wins. When several cut one class, the largest cut wins. A
scenario may add transfers (`[[scenario.transfers]]`) and rules. Their intakes and point sources go into every case,
as wells do.

**What the study supplies.** Each campus as a user: its CoolingWaterBmi config (design, plant, parameters), its
phases from the tracker records (date and MW), its weather station, and its source and return: the system it joins,
or its own well, intake or reuse. Wells and intakes need their cells and channels. The channels and cells to report.
The warm-up's length. The compared period's dates, which are on the run's calendar: a build-out in 2027 to 2035 run on
past weather sets the phases on the past dates that stand for those years.

Built and tested on the synthetic example (`tests/test_user_scenarios.py`). A phased user matches the component
built at each phase's size, day by day, to the bit. One user added alone and with a second, and together, with the
interaction. Every case's ledger closes every day, and the capture's residual is zero. The runner skips what is done
and resumes what isn't. The operating rules are tested by construction and on Osu_1hru with a reservoir below its
channel. A pump fills the reservoir at high flow, and its inflow is the channel's and the pump's. A release keeps
the channel up at low flow, never below its keep. A town's outdoor and indoor use is cut in two stages, booked apart
from a shortfall. A restart matches, bit for bit.

## Soil cover for soil evaporation: a process substitution (designed and built Nov 27, 2026, for the Virginia study)

Approved by Andrew on Nov 27, 2026. It is the first substitution built into Basin Ledger. The Clam study's soil hold under
snow (a study's `wrap_surface`, swatplus-bmi's settable soil temperature) was the pattern.

**Why.** SWAT+ 62.0.1's `et_act.f90` takes an HRU's soil cover as all its above-ground biomass plus its residue
(`pl_mass%ab_gr_com%m + pl_mass%rsd_tot%m`, kg/ha). Its potential soil evaporation is PET x exp(-5e-5 x (cover + 0.1)),
before the cap that keeps soil and plant ET within PET. A forest's wood counts as cover: 50 t/ha shades the soil to
nothing, so its soil hardly evaporates all year. A dormant grass, with LAI 0 and little residue, evaporates near PET in
winter. In the Rappahannock joint calibration, the observed forest-minus-open ET contrast was +93 mm a year and the
model gave +16, with the forest's ET parameters pinned at their bounds (the Virginia study's SWAT+ issue 5).

**The rule** (`[soil_cover]`, `coupling/soilcover.py`):

    [soil_cover]
    rule = "leaves"
    land_uses = ["frsd_lum", "frse_lum", "frst_lum", "wetw_lum"]   # and/or hrus = ["hru00001", ...]

- SWAT+ has no marker for a tree, so the config names the tree HRUs, directly or by land use. Each of the
  Rappahannock's forest land uses has a one-plant community (frsd, frse, frst, wetw).
- On those HRUs, before each day, Basin Ledger sets the cover. `"leaves"`: the leaf mass plus the residue, no wood.
  `"swat"`: the above-ground mass plus the residue, SWAT+'s own cover set back, a check that gives the same run.
- swatplus-bmi's settable `soil_surface_cover__mass-per-area_density` replaces `et_act`'s `cover` for that HRU for
  that day only. Every other HRU is left to SWAT+.
- Off without the section, so a run without it is unchanged.
- Unknown HRUs and land uses are refused. So is a swatplus-bmi without the soil cover.

**What it changes, and what it doesn't.**
- Only soil evaporation sees the set cover. SWAT+'s albedo and its erosion (`albedo.f90`, `ero_cfactor.f90`) compute
  their own cover, so a forest's wood still shapes them.
- Under snow (0.5 mm or more), SWAT+ takes half of PET whatever the cover.
- No new science goes into SWAT+: swatplus-bmi only lets the cover be read and set (Andrew's rule that it stays a
  wrapper). The rule lives in Basin Ledger and is stated in the config.

**When the masses are read, and how they are summed.** At the pause before the day, so they are the end of the
previous day. They are summed as SWAT+ sums them, in kg/ha and 4-byte reals, so `"swat"` is SWAT+'s cover to the bit
(a float64 sum in kg m-2 can round a tie the other way). On a run's first day SWAT+'s own cover is the residue alone,
since `pl_community` sums the plants' masses after `et_act`; the masses read before that day agree. `et_act` uses the
same masses unless the HRU's decision-table actions run earlier that day: a management schedule's auto operations,
with any fertilizer they schedule (`hru_control.f90`: `actions` before `et_act`). Scheduled operations
(`mgt_operatn`), growth and the year-end cut come after `et_act`.
- An HRU with no management schedule has no such actions. That is every Rappahannock forest land use (`mgt = null`
  in `landuse.lum`), so its cover is the same day's, exactly.
- On an HRU with auto operations, the cover is a day behind on the days they change the masses. Asking SWAT+ for the
  same day's masses would mean a pause for each HRU, inside its day. With 33,525 HRUs that costs too much, and it
  isn't needed for forests.

**Tested** (`tests/test_soil_cover.py`), on the synthetic example's HRU made a forest (Osu_1hru's `frst_tems_lum`,
growing from the start, its paddy and management schedule taken away):
- the section and its refusals; the HRUs named by name and land use; the wrapper on a stand-in;
- off: the surface model isn't wrapped. Two years of the synthetic and two-way examples give ledgers identical to
  0.18.0's;
- `"swat"`: the same run as off, to the bit, every day;
- `"leaves"`: the forest's ET rises and the whole basin still closes.

## Extending an existing MODFLOW 6 model (decided Nov 14, 2026; first step built)

**Decided (Andrew).** Basin Ledger gets a second way in: take a MODFLOW 6 model someone already has (a consultant's,
a utility's, USGS's, EGLE's), leave it as it is, and couple SWAT+ to it. SWAT+ supplies the land surface (snow, plants,
crops, soils, runoff, percolation, ET) and keeps the streams, where the in-stream questions are: phosphorus, nitrogen,
sediment and water temperature, through SWAT+'s quality routines. The first way in, where a study builds both models
(the Clam, Plaster Creek), stays and must keep working as well as it does.

**Why.** Groundwater work outside academia is done in MODFLOW; consultants rarely use SWAT. A calibrated MODFLOW model
usually already exists for any basin worth a large study, and a reviewer can open it. Extending that model answers
land-use, climate, water-quality and new-user questions without a second groundwater model to defend, and the review
question becomes what the coupling changed. Pitch: bring your MODFLOW 6 model; Basin Ledger adds SWAT+'s land surface
and in-stream water quality, and a ledger that proves every handoff conserves water.

**What exists.** `[modflow6] workspace` already copies an existing simulation into the run, and the driver couples
whatever packages a config names. The maps, the exchange and the checks were built and tested on simulations the
studies' own builders wrote (one daily stress period, our package names).

**Decisions, with what each needs.** Andrew's answers (Nov 14, 2026) are in italics after each.

- **The model as given.** Its grid (DIS or DISV), packages, stress periods, K, storage and boundaries are left
  untouched. Basin Ledger adds or replaces only what it couples: recharge (a recharge package it owns), the
  stream and lake exchange on the model's own RIV, DRN or SFR entries, and users' wells. Anything else in the model
  (its GHB, CHD, its wells, EVT) stays the model's, and the ledger books it as water from or to outside the coupled
  models, as it does now. A start-up check lists what Basin Ledger took over and what it left. *Yes.*
- **Maps on their grid.** HRUs to cells and SWAT+ channels to the model's stream entries, built from the model's grid
  and its package cell lists, with the area and length checks the studies' maps have. Stream entries no SWAT+ channel
  reaches are reported, not guessed. *Yes.*
- **One owner per water body.** SWAT+ owns the stream water and its quality. The model's RIV or SFR becomes the
  exchange face: its stage set each day from SWAT+'s channel, its flux returned to that channel. An SFR package that
  routes flow has its routing driven or switched off, so a river is never carried twice. Lakes follow the same rule.
  *Probably; Andrew wants to understand it better first.* To work through before building: what SWAT+ owning the
  stream means for a model whose SFR was calibrated as the river (its streambed conductance, its routed flows at
  gages), what "driven" routing looks like in MODFLOW 6's SFR, and a worked example on a published model.
- **Stress periods.** Consultant models often have monthly, seasonal or steady periods. Two modes, chosen per study:
  rewrite the model's time discretization to daily steps (the coupling as built), or keep its periods and pass SWAT+'s
  fluxes summed over each period. The ledger checks either. Restarts across several stress periods (built Nov 14,
  2026, under Restarts) became required. *Daily: Basin Ledger rewrites the model's time discretization to daily steps (built Nov 14,
  2026, below). Keeping the model's own periods stays an idea to think about, not built.*
- **MODFLOW 6 only.** Models in MODFLOW-2005 or NWT come in through USGS's converter (mf5to6) first; Basin Ledger
  doesn't read them. *MODFLOW 6 only, for now.*
- **Their calibration stands.** Groundwater parameters stay fixed by default; calibration covers SWAT+ and the coupling
  terms (streambed and lakebed conductance, the soil link). A study may free groundwater parameters, and says so.
  *Yes, but recalibrating the groundwater model may be needed once SWAT+ supplies its recharge: a setting in the
  config (groundwater parameters fixed, or free with their priors), not a code path.*
- **Solutes across the link.** Nitrate and chloride in the aquifer through the model's GWT (or one Basin Ledger adds),
  passed as loads into SWAT+'s channels with the water (Water quality across links). *Yes: groundwater carries
  nitrate and chloride to streams, and the in-stream questions need those loads.*
- **Confidentiality.** A client's model runs locally; nothing from it goes into public outputs or Water Ledger.
  *Yes.* Outputs, maps and ledgers from a client's model stay in the study's private folder; examples and tests use
  only published models.

**Open questions, answered (design pass, Nov 14, 2026).** Checked on MODFLOW 6.8.2 through modflowapi 1.0.1, on a
two-layer test model with a lake, a multi-aquifer well and three SFR reaches; `tests/test_existing_api.py` keeps
each answer as a test, so a new MODFLOW 6 that changes one fails there first.

1. **Can a coupler set LAK, MAW and SFR inputs each step?** Yes, all of them, between `prepare_time_step` and the
   solve, the same way the driver sets RCH and WEL now.
   - LAK: `RAINFALL` and `EVAPORATION` (m d-1 over the lake's surface), `RUNOFF`, `INFLOW` and `WITHDRAWAL` (m3 d-1),
     per lake. Each reaches the lake's budget that day (rainfall 0.01 m d-1 on a 4 ha lake books 400 m3 d-1) and
     the stage (`XNEWPAK`; `STAGE` is the input stage of a constant-stage lake) responds. So a lake can be coupled as
     its own package: SWAT+'s rain, evaporation and runoff in, the lake's stage and outflow from MODFLOW 6.
   - MAW: `RATE` per well. MODFLOW 6 splits it among the well's connections by their heads (`SIMVALS`, per
     connection) and reports what it pumped (`RATESIM`; less than `RATE` if a shut-off or reduction length is set).
     The connections add to the rate less the change in water held in the wellbore (0.02 of 200 m3 on the test day),
     so the ledger's pair for a MAW well is rate = connections + wellbore storage, read from MODFLOW 6. A utility
     well screened across layers can be a MAW well.
   - SFR: `INFLOW`, `RUNOFF`, `RAIN` and `EVAP` per reach (`INFLOW` is what the SFR exchange sets now), and also the
     reach's status (`IBOUNDPAK`: 1 active, 0 inactive, −1 simple), its specified stage (`SSTAGE`) and each reach's
     upstream fraction (`USTRF`).
   - What the API doesn't allow: adding or removing entries, or changing which reaches connect (`IA`, `JA`). A
     package's size and topology are fixed when it is read.
   - One caution: MODFLOW 6 reads a package's period block at the start of a stress period. A value set from Python
     lasts until the next block the model has for that package, so the coupler sets its values every day, after
     `prepare_time_step`, as it does now.
2. **How SFR computes flow, and whether a reach stores water.**
   - Without the `STORAGE` option (MODFLOW 6's default, and what most published models use), a reach holds no water
     from one step to the next. Within each outer iteration SFR goes downstream reach by reach: a reach's inflow
     (its upstream reaches' outflow times its upstream fraction, plus `INFLOW`, `RUNOFF` and `RAIN`), less `EVAP` and
     its leakage to the aquifer (or plus its gain), is its outflow. Its depth comes from that flow by Manning's
     equation, and its leakage from that depth and the cell's head. A reach can't lose more than reaches it. Each
     step is a steady routing of that step's water.
   - With `STORAGE` (kinematic-wave routing, MODFLOW 6.5 and later), reaches hold water between steps, and a
     restart would need it saved. The exchange refuses it now and still would.
   - **Driving the routing** means three things MODFLOW 6 supports, each tried on the test model:
     - *inflows set*: each SWAT+ channel's water goes in as its top reach's `INFLOW` (what the exchange does now);
     - *connections cut at run time*: a reach's `USTRF` set to 0 takes nothing from its upstream reaches, and their
       outflow is booked as SFR's external outflow. SFR's budget still closes (0.00% on the test), and the model's
       files are untouched. The cut is per reach, though, not per connection: at a confluence where a reach takes
       water from two channels, both are cut;
     - *simple status*: a reach with status simple takes its stage as given (`SSTAGE`, set each day) and no longer
       computes it from flow. It still routes the water it has, and its leakage is held to its inflow (stage 18 m
       over a 15 m head with 100 m3 d-1 in: it lost 100 m3 d-1 and its stage dropped to the streambed top). That is
       RIV with the supply limit built in.
   - **Switching the routing off** means the water SFR routes isn't counted anywhere outside SFR: what leaves a
     channel's last reach goes to SFR's external outflow, and SWAT+ is the model that carries the river on. The
     exchange (each reach's leakage or gain) is what crosses to SWAT+.
3. **One owner per water body, for a model whose SFR was calibrated as the river.** *Andrew's leaning (Nov 14,
   2026), not final: "probably it's A and B"* (worked out below the options). In the model as calibrated, SFR is the river: its reaches carry the flow the modeler set at
   the network's upstream ends (and runoff, often from a soil-water-balance model), each reach's stage comes from that
   flow and the reach's width, slope and roughness, and its streambed conductance was calibrated so heads and the
   routed flows at gages (often low flows, or base flow separated from them) match. Coupling SWAT+ takes over the
   flow at the top of each channel and the recharge. What stays the model's is the streambed: its conductance,
   elevation and geometry. The options:
   - **A. SWAT+ owns the water; SFR active within each channel** (the exchange as built, extended to a foreign
     SFR). Each channel's water goes in at its top reach, SFR routes it down that channel's reaches with its own
     Manning's equation and bed, and connections between channels are cut at run time (`USTRF` = 0), so SWAT+ alone
     carries water from one channel to the next. Keeps the calibrated bed and the model's own stage-flow relation.
     Gage flows come from SWAT+, which must be calibrated to them; SFR's routed flows aren't a result any more.
     Costs: a confluence cut by `USTRF` also cuts the reach's own channel where the two networks disagree; such
     reaches are reported; rewriting those connections in the run's copy of the package (the model's own files
     untouched) is the way around it, not built. Inflows the modeler set at upstream ends, and SFR's own runoff, are replaced
     by SWAT+'s and listed by the start-up check.
   - **B. SWAT+ owns the water; SFR as a fixed-stage face** (status simple). Each reach's stage is set each day from
     SWAT+'s channel depth (its own rating curve) over the reach's streambed top, and its channel's water goes in at
     the top reach so the leakage stays within the flow. No cuts needed beyond those of A, and the stage follows the
     model that carries the water. Gives up SFR's calibrated width, slope and roughness (the stage-flow relation);
     the streambed conductance stays. Closest to how RIV is coupled now.
   - **C. MODFLOW 6 owns the water; SWAT+ supplies it.** SWAT+'s runoff and lateral flow go into SFR's reaches as
     `RUNOFF`, the upstream inflows as `INFLOW`, and SFR routes the river as calibrated. Keeps the model's routed
     flows at gages. But SWAT+'s channels no longer carry the water, so its in-stream phosphorus, nitrogen, sediment
     and temperature routines have nothing to work on, which was the reason for the second way in. Listed so the
     trade is plain, not as a candidate.
   - What A and B share, and what a reviewer should be told: the model's bed stays calibrated; the recharge and the
     water in the streams change, because SWAT+ supplies them now; the groundwater calibration may need revisiting
     (the setting for groundwater parameters fixed or free, above). A RIV or DRN model has no routing to give up, so
     A and B only matter for SFR. Lakes follow the same choice: LAK active with SWAT+'s rain, evaporation and runoff
     set each day (point 1), or the lake as SWAT+'s reservoir with LAK's stage set and its routing (outlets) off.
   - **Decided (Andrew, Nov 14, 2026): A and B, B chosen per study and applied per channel.** Built as below.
   - **The leaning, A and B, worked out.** Extending an existing model is one use of Basin Ledger and a test of its
     design; building both models (the Clam's path) may be the main one. So the choice is made for both paths at once:
     - **A is the main option in both.** It is the structure the exchange already has when Basin Ledger builds SFR
       from SWAT+'s channel lines (`channel_reaches`, `coupling/sfr.py`): each channel its own run of reaches, SWAT+
       carrying water from one channel to the next. Nothing changes in the build-both path.
     - **For an existing model, the networks are made to line up** by building SWAT+ on the model's SFR network: its
       stream lines burned into the DEM, and a SWAT+ outlet at every SFR confluence. That is the study's work, in
       QSWAT+ or its own scripts (principle 8); Basin Ledger lists the confluences and where they are
       (`gridmaps.sfr_confluences`). With the networks lined up, every connection between channels runs from one
       channel's last reach to another's top reach, so setting the top reaches' upstream fraction to 0 cuts exactly
       the water SWAT+ carries, and the per-reach cut at a confluence (point 2) is no longer a problem.
     - **One line-up rule for both paths, checked at start-up** (`sfr_problems`, built): each SFR reach maps to one channel
       (or none, and stays the model's); each channel's reaches form one connected path from one top reach; every
       confluence falls on a channel boundary, which means a connection between two channels goes from the first's
       last reach to the second's top reach, and no reach inside a channel takes water from another channel or from
       a reach no channel has. A coupled channel flowing into a reach no channel has is refused too: SFR would carry
       SWAT+'s water on. A mismatch is refused, by reach and channel name. A package Basin Ledger built from channel
       lines has no connections between channels and passes as it does now.
     - **B for what can't be lined up** (built). Chosen per study, applied per channel. `[exchange] sfr = "route"`
       (A everywhere; a mismatch is refused, the default), `"route_or_stage"` (A where a channel lines up, B for a
       channel that doesn't, each named), or `"stage"` (B everywhere). Per channel rather than per reach, because A
       routes a channel's water down all of its reaches, and a channel half routed and half at a fixed stage would
       have a stage that jumps between reaches. B still takes the channel's water at its top reach (a reach at a fixed
       stage still routes what reaches it and loses no more), so it needs the same cuts as A, and a channel that can't
       be cut at its boundaries can't run as B either. The start-up check lists each coupled channel with A or B and,
       for B, the line-up problem that put it there.
     - **How B runs** (`StreamExchange.set_inflow`): each of the channel's reaches takes status simple, its stage
       its streambed top (`STRTOP`) plus SWAT+'s channel depth for the day just run (`channel_water__depth`, booked
       as `exchange:<channel>:sfr_depth_m`). Its water goes in at its entries, the reaches nothing in the channel
       flows into and the reaches cut from water coming from outside it (a confluence inside the channel), each by
       the share of the channel's reach length from it down to the next entry. Each reach still routes what reaches
       it and loses no more. The connections into the channel are cut wherever they fall. Still refused under B: a
       coupled channel flowing into a reach no channel has, and the package settings (units, storage, diversions,
       the mover). A loss that, split this way, comes to more than SWAT+'s channel has the next day is the
       exchange's shortfall, booked as with A and RIV.
4. **The test model.** No USGS-published MODFLOW 6 model of a Michigan basin turned up (searches of USGS data
   releases and ScienceBase, Nov 14, 2026). Candidates, none downloaded:
   - **Kalamazoo River basin (HUC8), Michigan**: the steady MODFLOW-NWT inset model cut from USGS's Lake Michigan
     Basin regional model (Feinstein and others, SIR 2018–5038, doi:10.3133/sir20185038; files doi:10.5066/F76D5R5V,
     edited versions doi:10.5066/P9HS83JL, which also has Boardman-Charlevoix). A Michigan glacial basin of the right
     kind, but it would come in through mf5to6, and being steady it has no storage properties: daily steps need
     specific yield and storage added, which is the study's choice, not Basin Ledger's.
   - **Wisconsin Central Sands**: MODFLOW 6 inset models around lakes on glacial outwash, transient 2012–2018
     (Fienen, Haserodt and Leaf, doi:10.5066/P9BVFSGJ; report SIR 2022–5046). Native MODFLOW 6, transient, with lakes
     and streams, built by another hand (modflow-setup); the closest to what a consultant would bring. From
     published descriptions only (the release's metadata names no packages, and its full metadata and the report
     didn't load): Leaf and Fienen's modflow-setup paper (Frontiers in Earth Science, 2022,
     doi:10.3389/feart.2022.903965) describes the Pleasant Lake inset with LAK (lake extent, bathymetry, lakebed
     leakance zones), SFR from SFRmaker, the mover (MVR) keeping SFR's flow continuous across the inset's boundary,
     recharge from SWB and irrigation wells. Its first period is steady and the rest are monthly. The refined inset is
     a second GWF model in the same simulation, joined to its parent by a GWF-GWF exchange solved in the same matrix.
     Basin Ledger couples one GWF model per simulation and refuses the mover, so this model needs both before it can
     be coupled whole. That is a reason to choose it as the test (consultants use insets and the mover too), or to
     start with the intermediate inset alone, as a single GWF model. The archive's files still need checking.
   - **CONUS surficial models** (Zell and Sanford, doi:10.5066/P91LFFN1): 75 steady, one-layer MODFLOW 6 models, one
     covering the Lower Peninsula; streams as drains. A real foreign MODFLOW 6 grid over Michigan for the maps, too
     simple for the exchange.
   - **Michigan Hydrologic Framework** (EGLE, MSU and USGS): a platform hosting Michigan groundwater models and
     inset extraction. If Andrew has access, it may hold the Michigan MODFLOW 6 model the search didn't find.
   *Chosen (Andrew, Nov 14, 2026): the best one first, more than one likely.* The first is the Central Sands
   **Pleasant Lake model** (16 MB; `examples/central-sands`, `tests/test_central_sands.py`). The second, *yes
   (Andrew)*, is the **Kalamazoo model** through mf5to6, as the Michigan example (below). What Pleasant Lake holds, from its files:
   - **Two groundwater models in one simulation.** A parent (5 layers, 100 × 100, CHD edges, 244 irrigation wells,
     551 SFR reaches) and an LGR inset of 20 m cells around the lake (5 layers, 100 × 120, 25 SFR reaches, one lake
     of 1,382 connections with no outlets), joined by a GWF-GWF exchange with ghost nodes.
   - **The mover**, in that exchange, carrying the inset's two outlet reaches into the parent's SFR.
   - **Recharge from SWB in both models, Newton, and a steady first day** followed by 84 monthly periods of five
     steps growing by 1.2.
   - **Meters and days**, on a structured grid in Wisconsin Transverse Mercator.

   What Basin Ledger does with it already: reads it as given; lists both models' packages; and solves its steady
   first day. It puts both SFR networks through the line-up rule with SFRmaker's lines (`line_id`, NHD flowlines) as
   SWAT+'s channels, and they line up: 2 lines in the inset, 48 in the parent with 21 connections cut at channel
   boundaries, and no channel problems. So SWAT+ built on the same flowlines would route every channel (A), which
   supports building SWAT+ on a model's own stream network. It lists the parent's 10 confluences. And the mover can
   be cut at run time: its factor (`MVR/VALUE`) set to 0 through the API sends the inset's outlet water out of SFR
   (1,167 and 2,493 m3 d-1 on the first day) instead of into the parent, as a cut upstream fraction does.

   What it needs before a coupled run, in order:
   1. Daily steps (built, below): its 85 periods become 2,558 days, its steady first day still a steady day.
   2. The mover allowed where it only joins coupled reaches, cut each day like the line-up rule's connections
      (built, below): on the parent's side its two moves land on top reaches and are cut; on the inset's they take
      coupled streams to another model, refused until both models are coupled (3).
   3. Several groundwater models in one simulation: each with its own maps, recharge and exchange, and the
      GWF-GWF exchange's flow inside the ledger's groundwater, not a boundary. *First step built (below): one model
      of such a simulation coupled, the exchange booked as a boundary. Coupling them together is proposed below, for
      Andrew to decide.*
   4. The lake coupled, as streams are (below, "Lakes").
   5. The parent's irrigation wells, MODFLOW first (below, "Large users and irrigation").

**The Kalamazoo model, converted (Nov 14, 2026; `examples/kalamazoo`, `tests/test_kalamazoo.py`).** USGS's
Kalamazoo inset of the Lake Michigan Basin model (doi:10.5066/F76D5R5V), downloaded by hand (ScienceBase puts a
captcha in front of files over 1 GB) and converted with mf5to6 from the MODFLOW 6.8.2 release in 13 seconds. The
converter needed a copy of the name file with Unix line endings and without MNWI (MNW2's reporting), whose output
files it can't open; the model's input is untouched. What it holds:

- 8 layers × 730 × 1,230 cells of 500 ft, 7.2 million cells, one steady period with no storage package;
- SFR of 20,010 reaches, RIV (52,969 entries), DRN (29,897), GHB (137,100, its regional boundaries), fixed heads
  from IBOUND (28,424), recharge as arrays, and 710 multi-node wells (MNW2) as MAW;
- **feet and days.** The converted DIS says `LENGTH_UNITS UNDEFINED`, its reference file says feet (UTM zone 16 in
  feet, origin there and not in the DIS), and SFR's Manning's constant is for feet (1.486 × 86,400).

What it needs, in order:

1. **Units (decided, Andrew, Nov 14, 2026: convert at each link; built, below).** Basin Ledger took meters. Most
   consultants' models in the US are in feet, so this matters beyond this model. Of the two ways, Basin Ledger
   converting at each link (SWAT+'s m3 and m to the model's ft3 and ft as they cross, the ledger in m3 throughout)
   was chosen over the study rescaling the model to meters: it keeps the second way in's promise, the model as
   given, and its reviewers' numbers with it. Its DIS leaves the length units undefined, so a study states them
   (`[modflow6] length_units = "feet"`); its SFR's Manning's constant passes the check for feet.
2. **Storage** for daily steps: specific yield and specific storage, which a steady model never needed. They are the
   study's [A], and calibration may free them (the setting for groundwater parameters fixed or free). Built
   (Nov 15, 2026, `existing.add_storage`): `[modflow6] storage = {sy = ..., ss = ...}` gives a model with no storage
   package one of Basin Ledger's own in the run's copy (ICONVERT, 1 unless given, SY, and SS given per meter and
   written in the model's length units, read from its DIS options or the config's), transient from the first period;
   `storage = {}` makes a model's own steady periods transient, its own values kept (values given beside its own are
   refused). Without it, a model with steady days is noted: MODFLOW 6 solves those days without storage, each the
   day's steady state. Tested on the synthetic model handed over steady both ways, coupled 75 days, and on
   Kalamazoo's files (text only: its name file and the STO added, SS per foot).
3. **Its size.** A day of 7.2 million cells is minutes of solve, a year days. A study would cut an inset around its
   question, as USGS cut this one from the regional model.
4. **MAW wells** as users' wells (built Nov 15, 2026): a `[[wells]]` entry names a multi-aquifer well by its number
   in a MAW package (`well = 3`, in place of `cell`). Each day the user's withdrawal is the well's `RATE`; MODFLOW 6
   splits it among the well's connections by their heads, and what it pumped (`RATESIM`, less than asked only with
   a shut-off) is what was delivered, booked with where it came from: the aquifer by the connections
   (`route:well:<name>:from_aquifer`) and the wellbore's water (`from_wellbore`, MAW's `QSTO`), which add to it to
   MAW's own solve (a few parts in a million). The well's connections are Basin Ledger's in the start-up check;
   `well` naming a package that isn't MAW, a number it doesn't have, or a well given both a cell and a number is
   refused; transport with a MAW well isn't built. Tested on the synthetic foreign model with a farm well screened
   through both layers, pumped by a second plant for 75 days (`tests/test_existing.py`).

**Lakes: the stream precedent (decided, Andrew, Nov 14, 2026).** SWAT+ owns a lake's water and quality, as a
reservoir, and MODFLOW 6's LAK is the face it exchanges through, as SFR is for a stream. Two forms, as for channels,
chosen per study and applied per lake:

- **A, LAK active within the day.** Each day the lake starts from SWAT+'s stage (`XOLDPAK`, set after
  `prepare_time_step`), takes SWAT+'s inflow, rain and evaporation for the day, and MODFLOW 6 solves its stage and
  leakage within the day on the lake's own bathymetry and lakebed leakance, as the model's builder calibrated them.
  Its outlets are cut (an outlet into SFR is a connection between water bodies, which SWAT+ routes), and its net
  exchange goes to SWAT+'s reservoir the next day. What LAK's stage did within the day isn't kept: tomorrow starts
  from SWAT+'s stage again, so the lake's water is held once, in SWAT+.
- **B, LAK at a fixed stage.** Status constant, its stage SWAT+'s each day. LAK computes the leakage at that stage,
  and the loss is limited to what the reservoir holds, as RIV's is.
- **Checked through the API** (MODFLOW 6.8.2, `tests/test_existing_api.py`): both work. A lake at constant status
  takes the stage given and leaks at it, and an active lake restarted from a given stage each day takes its
  inflows and leaks within the day from there.
- **Built** (`coupling/lakes.py`, `[[lakes]]`: reservoir, package, lake, form; swatplus-bmi 0.11.0, request 6).
  SWAT+ 62.0.1 keeps a reservoir's volume and area, not its level, so the stage is LAK's own at SWAT+'s volume:
  from the lake's table if it has one, else from its connections as MODFLOW 6 computes a lake's volume
  (`lak_calculate_vol` in 6.8.2: each connection's surface area, smoothed at its top, times the depth over its
  bottom, from `TELEV`, `BELEV` and `SAREA` in memory), inverted by bisection. Over a day the change in that volume
  is MODFLOW 6's own lake storage (`QSTO`) to 1e-9, on the test model and on Pleasant Lake's 1,382 connections
  (`tests/test_central_sands_lake.py`). Each day: yesterday's leakage is set on the reservoir before
  SWAT+'s day; after it, the reservoir's day is read (volume, inflow, outflow, precipitation, evaporation, its own
  seepage, and the exchange as SWAT+ took it); the lake's own rainfall, evaporation, runoff, inflow and withdrawal
  are set to 0, then A sets it active from the stage of the day's starting volume with the day's inflow and
  precipitation as `RUNOFF` and its outflow, evaporation and seepage as `WITHDRAWAL` (m3, so no area enters), and B
  sets it constant at the stage of the day's ending volume. After the solve the lake's leakage (its connections'
  `SIMVALS`, positive into the aquifer) is the reservoir's exchange the next day. The ledger books
  `lake:<reservoir>:due` against `applied` (what SWAT+ took) plus `shortfall` (a loss the reservoir couldn't
  supply), and `to_aquifer`, `stage_m`, `volume` and, for A, LAK's own `mf6_end_stage_m`, which isn't kept. A
  coupled lake's outlets are cut each day (a specified outlet's `OUTRATE` and a weir's or Manning's `OUTWIDTH` set
  to 0 after the day is prepared), since SWAT+ routes the reservoir's release, and the ledger checks nothing left
  through them (`lake:<reservoir>:outlet_flow`, zero). The start-up check lists the lake's connections as Basin
  Ledger's, and refuses an outlet bringing another lake's water into a coupled lake unless that lake is coupled too
  (water SWAT+'s reservoir doesn't see), a reservoir whose own seepage (`k` in hydrology.res) isn't 0, which would
  leak twice, and a package that isn't LAK or a lake it doesn't have. A restart carries the leakage due. Tested on the synthetic model with a pond
  in its top layer's inactive ring (its bed on the bottom layer) and Osu_1hru with a reservoir below its channel,
  both forms, 75 days: the pond, below the aquifer's heads, gains 130 m3 a day falling to 79 as it fills, each
  day's gain reaches the reservoir the next day, and the ledger closes. Pleasant Lake, a seepage lake with no
  outlets, is the first study test; it needs its two models coupled together first. Because SWAT+ represents such
  a lake as a reservoir or a wetland, the study's SWAT+ build decides which object it is; a wetland isn't built.

**Large users and irrigation, MODFLOW first (Andrew, Nov 14, 2026).** Basin Ledger exists to answer what a large
water user does to a basin. On the second way in, groundwater is the model's, so pumping is too: every groundwater
withdrawal, a data center's, a city's or a farm's, comes out of MODFLOW 6's wells, and SWAT+ never pumps its own
aquifer.

- **Large users** already work this way: a user's withdrawal from a well is set on that well's entry in MODFLOW 6
  (taken over by the start-up check), and its return goes where its route says. A new user's well is an entry the
  model must have; adding wells to a model that lacks them is still to build.
- **Irrigation.** SWAT+ decides when and how much a field is irrigated (its decision tables, from the weather and the
  soil), and MODFLOW 6 pumps it. The irrigate action takes its water from an outside source in SWAT+, and Basin
  Ledger pumps the same volume that day from the wells the study maps to that field's HRUs. So the demand responds
  to climate and land-use scenarios, the water is pumped once, from the model that owns groundwater, and it lands on
  SWAT+'s fields, where its return flow and its losses are SWAT+'s.
- **Built** (`coupling/irrigation.py`, `[irrigation] wells`; swatplus-bmi 0.11.0 reads each HRU's irrigation,
  request 7). After SWAT+'s day, the water the day's irrigation put on each mapped HRU (what reached the soil and
  what ran off; an application's losses never enter SWAT+'s balance) is requested from the HRU's wells in their
  shares, and MODFLOW 6 pumps it that day. The ledger books `irrigation:pumped` against the wells' shares, and each
  well's request against what MODFLOW 6 delivered, as for a user's well (a well can serve both). Irrigation on HRUs
  not mapped stays water from outside (`irrigation:from_outside`). The start-up reads each HRU's sources from its
  land use's schedule and decision tables (hru.con, hru-data.hru, landuse.lum, management.sch, lum.dtl): only an
  outside source is pumped (`irrigate` from `null`, `irr_demand` from `unlim`, a scheduled `irrm`); a mapped HRU
  that irrigates from a SWAT+ channel, reservoir or aquifer, or through water allocation (its source unknown here),
  is refused, since that water already comes out of SWAT+. Tested on the synthetic model through Osu_1hru's first
  irrigated days: 60 mm on the 10 ha paddy on 8 May, 6,000 m3, pumped that day from two added wells 60:40, the
  ledger closing every day.
- **Keeping a model's own pumping** (its calibrated irrigation rates, as Pleasant Lake's 244 wells have from SWB) is
  the other way: the wells stay the model's, as the start-up check lists them, and SWAT+ doesn't irrigate those
  fields. It reproduces the model as calibrated, but the water pumped never reaches SWAT+'s fields, and a scenario
  can't change the demand. It suits a first run against the model's own calibration, not the large-user questions.
- **Refused now** (built, `swatplus.aquifer_withdrawals`, in two-way mode): a decision table a schedule uses that
  irrigates from SWAT+'s aquifer, and a water allocation point of diversion from an aquifer. In two-way mode
  SWAT+'s aquifer only holds the recharge it passed to MODFLOW 6, so water drawn from it would be pumped twice.
  Irrigation from a channel or a reservoir is surface water and stays SWAT+'s.

**Irrigation from surface water (designed Nov 24, 2026; way chosen by Andrew the same day; built Nov 25).** The aim: a
field that SWAT+ irrigates takes its water out of the stream or reservoir it is drawn from, as `[irrigation] wells`
takes it out of MODFLOW 6. The ledger has to close, and no water may be counted twice against SWAT+'s own irrigation
accounting.

**Chosen: SWAT+'s own irrigation from a channel or reservoir, booked by Basin Ledger.** A decision table's `irrigate`
action from `cha` (or `sdc`) or `res` (actions.f90) runs the same day. It takes the demand out of the channel's
storage (`ch_stor`) or the reservoir, all or nothing: if the source holds less, the field gets none that day. So
SWAT+ moves the water, and in low flow a dry channel gives nothing. Basin Ledger only has to see it. Today it can't:

- swatplus-bmi reads the irrigation that reached the soil and ran off (request 7), not what left the source;
- the application's loss (the demand less those two, `(1 - eff)` of what isn't runoff) leaves SWAT+ with no account;
- `channel_sd_day` has no column for the withdrawal, so the whole-basin channel balance puts it in its residual;
- with `[irrigation]` on, the driver books it as `irrigation:from_outside`, which is wrong today.

There is also a SWAT+ 62.0.1 bug, found by the swatplus-bmi session on Nov 24, 2026. A channel's `irrigate` lowers
`ch_stor`, but routing (ch_rtmusk.f90) starts from `tot_stor`, carried from the day before, and overwrites `ch_stor`
from it. So the field gets the water and the channel never loses it. Reservoirs do lose theirs. Andrew chose to patch
it in swatplus-bmi 0.14.0 (actions.f90 and wet_irrp.f90, drafted upstream). Basin Ledger therefore needs 0.14.0 for
channel irrigation. On an older swatplus-bmi, channel irrigation is water from outside the basin.

So it waits on a read-only swatplus-bmi variable set (docs/requests.md, request 8; in swatplus-bmi 0.14.0 as
`channel_water__irrigation_withdrawal_volume_flow_rate`, `reservoir_water__irrigation_withdrawal_volume_flow_rate`
and `land_irrigation_water__withdrawal_volume_flow_rate`). By channel and by reservoir, the
day's water SWAT+'s irrigate actions took; by HRU, what its irrigation took from a SWAT+ water body. SWAT+ doesn't
keep either one tagged with its source (`irrig(j)%demand` is the action's demand, whatever the source). So
swatplus-bmi's patch adds them up where actions.f90 removes the water: a diagnostic, with no change to what SWAT+
does.

**Built (Nov 25, 2026, on swatplus-bmi 0.14.0; `coupling/irrigation.py`, `NativeIrrigation`).** Basin Ledger:

- **needs no config.** The HRUs whose decision tables irrigate from a channel or reservoir, and the water bodies they
  name, are read from the project (`irrigation.native`). `[irrigation] wells` still refuses those HRUs, so a field's
  water comes from one source.
- **books each day**, over the HRUs, `irrigation:native:withdrawn` (what they took),
  `irrigation:native:from_channels` and `:from_reservoirs` (what the water bodies gave), `:applied` (what reached
  the soil and ran off, on the HRUs that took water that day) and `:loss`. It also books `irrigation:channel:<c>`
  and `irrigation:reservoir:<r>` for each water body the tables name. Two pairs close each day, to float32: what
  was taken = what the channels and reservoirs gave, and what was taken = applied + loss. `:excess` (more applied
  than taken) must be zero. With `[irrigation] wells`, `irrigation:from_outside` leaves out the HRUs that took
  from a SWAT+ water body that day.
- **in the whole-basin balance** (`balance.py`), the channels' and lakes' outflow side gains
  `irrigation_from_streams` and `irrigation_from_lakes`, so their residual closes, and `irrigation_loss` is
  reported. The land keeps SWAT+'s `irr` (only `irrig%applied`, hru_control.f90). Of it, what the streams and lakes
  gave less the loss is the basin's own water; the rest came from outside or MODFLOW 6's wells. Water taken from a
  channel after that channel has routed its day shows in its storage the next day. The balance sums over the run,
  so only the run's last day can be off by that day's withdrawal.
- **at start-up**, a project with channel or reservoir irrigation on an older swatplus-bmi gets a warning that it isn't
  booked. The warning becomes a refusal where `[outputs] balance` is asked for.
- **not checked yet:** the constituents. `irrig(j)%water` is `rto * ch_stor`, nitrate and phosphorus included.
  Whether SWAT+ puts them on the field, and how that meets the stream transport's mass, needs checking before a
  study runs transport with channel irrigation.

**Not chosen.**

- An outside source in SWAT+ with a Basin Ledger place of use paying it back the next day. It needs no swatplus-bmi
  change, but the stream pays a day late, and in low flow the field gets water the stream never had.
- Water allocation's `irr` place of use (wallo_start_day, wallo_pou_deliv). It doesn't work in 62.0.1:
  wallo_pou_deliv tests `irrig(j)` before setting `j` and indexes its loop by `ihru`, not `ird`, and wallo_start_day
  sums the demand with `hru(ihru)%area_ha`.

Water allocation's places of use, Basin Ledger's intakes, have their own SWAT+ bug, found the same day. wallo_withdraw.f90
reads a place's demand as `pou(ipou)%pod(ipod)%duty`, where `ipou` is the place's position in its point of
diversion's list and `ipod` is the point of diversion's number. So it delivers nothing, or a stray value, to a place
of use on any channel but the first, and on a shared channel each place gets another's demand (request 5, now
confirmed). The chosen way doesn't go through water allocation.

**Tested** (`tests/test_irrigation_surface.py`). Osu_1hru's paddy table is switched to irrigate 1 mm (100 m3)
from its channel with sprinkler_med (efficiency 0.85), and the channel's storage coefficient is halved so it keeps
some water. It is coupled both ways on the synthetic example for three years, with the whole-basin balance.

- Thirteen actions found water, 100 m3 each.
- The channel gave what the HRU took every day, and what reached the soil was 85% of it.
- SWAT+'s land budget's irrigation was that 85% (1,105 m3), with none from outside.
- The channels closed to -113 m3 against 611,000 m3 of outflow. Without the irrigation term they would be off by
  about 1,190 m3.

swatplus-bmi 0.14.0's own tests cover:

- the channel's outflow and storage against a run without irrigation;
- a day the channel is short, which gives nothing;
- a reservoir as the source.

**Built (Nov 14, 2026): the start-up check and maps on the model's own grid.**

- **What is taken over** (`coupling/existing.py`; `Basin.takeover`, `basinledger check`). On the first day, once
  MODFLOW 6's packages have their entries, each package with flows is listed with how many of its entries a link
  sets, and which link: recharge from SWAT+ (every entry of `[recharge] package`), the stream exchange (the entries
  of `[exchange] package` the map gives a channel), seepage (its drains), groundwater up into the soil, and users'
  wells. The rest of each package is the model's, booked each day as `groundwater:boundary:<package>`, as before.
  The recharge Basin Ledger sets goes into the model's own recharge package, named in the config, so its calibrated
  recharge is replaced rather than added to; a model without one gets a recharge package of Basin Ledger's own on
  the map's cells with `[recharge] add = true` (built Nov 14, 2026). A user's well is an entry the model has (a city well already in it), or, with `add = true`, a
  well Basin Ledger adds in a WEL package of its own (built Nov 14, 2026; below, "wells the model doesn't have").
- **What would count water twice is refused:** another recharge package left as the model's where it recharges ground
  SWAT+'s recharge covers (on each cell, its own share beyond what SWAT+'s map leaves; 0.8.1), a UZF package (its own
  infiltration), and EVT with groundwater up into the soil. Noted, not refused: EVT otherwise (groundwater ET stays
  the model's), and SFR or LAK left as the model's while SWAT+'s channels are coupled (the same river carried by
  both models, if they stand for the same streams).
- **Maps on the model's grid** (`gridmaps.py`). `hru_map` places each HRU or aquifer polygon on the grid's columns by
  area, on each column's top active layer (inactive cells are common in such models), and reports ground outside the
  grid, over columns with no active cell, and slivers under 0.001 m2 (rounding along cell edges on a rotated grid).
  `stream_entries` reads a package's entries: an SFR package's reaches, or the distinct cells of a RIV or DRN
  package over all its stress periods. `channel_map` gives each entry to the channel whose line runs furthest
  through its cell. An entry no line crosses is reported and left as the model's; a study may set `reach_m`, a
  distance within which a line reaches an entry, and those are reported apart. A cell two channels cross goes to the
  longer and is reported. The length check compares each channel's line inside the model's active area with its
  length through the cells it was given, and, for SFR, with its reaches' own lengths. Both write the CSVs the coupling
  reads, in the grid's own cell form.
- **Tested on a synthetic model in someone else's style** (`examples/synthetic/existing.py`, `existing.toml`): a
  12 × 12 two-layer grid, rotated 15° at UTM-like coordinates, the top layer's outer ring inactive (10 ha active,
  Osu_1hru's aquifer), monthly stress periods of daily steps, its own recharge, a creek with a tributary no channel
  reaches, tile drains, a general-head and a fixed-head boundary, two city wells and riparian ET. Basin Ledger takes
  over the recharge, 10 of the creek's 12 entries and one city well. Coupled both ways for 75 days, the ledger closes
  every day across three stress periods, and DISV gives DIS's heads to 1e-6 m. With a month of one step, or of
  five growing, the run rewrites the time discretization to daily steps first and gives the same run (below,
  "daily steps").
- **The SFR line-up rule** (`sfr_problems`, the same in both paths): connections between channels are allowed from
  one channel's last reach into another's top reach, or from the model's own reaches into a top reach, and are cut
  each day by setting that top reach's upstream fraction to 0. Nothing new enters the ledger: the cut water leaves
  SFR as its external outflow, and SWAT+ carries it. Refused, by reach and channel: a confluence inside a
  channel, a channel leaving SFR from any reach but its last, a coupled channel flowing into a reach no channel has,
  and a channel that isn't one chain. A package built from SWAT+'s channel lines has nothing to cut and runs as
  before. `gridmaps.sfr_confluences` lists the reaches more than one reach flows into, with their cells and places,
  for the outlets a study puts in SWAT+. Tested on a hand-built network with a tributary (`tests/test_sfr_lineup.py`):
  lined up, each channel's reaches carry only that channel's water.
- **B, a fixed stage from SWAT+'s depth** (`[exchange] sfr = "route"`, the default, `"route_or_stage"` or
  `"stage"`): as described under point 3. Tested by hand on the tributary network (a channel with a confluence inside
  it runs as B, its water in at its top and at the confluence by length; every channel at a stage loses no more than
  went in), and coupled on Osu_1hru for a year: the ledger closes, each reach's stage is its streambed top plus
  SWAT+'s depth (or its top, where it lost all it had), and no day's loss is more than the channel's inflow.

**Built (Nov 14, 2026): a model in its own units, converted at each link** (Andrew: "Let's build that option").

- **The model stays as given.** Its units are DIS's `LENGTH_UNITS` (feet, meters or centimeters) and TDIS's
  `TIME_UNITS` (seconds, minutes, hours or days). Where either is undefined, as in Kalamazoo's conversion, the
  study states it (`[modflow6] length_units`, `time_units`) and a run that doesn't is refused at start-up. A unit
  stated otherwise than the model's own is refused, and so are years (a day isn't a whole step of them).
- **One place converts** (`Modflow6.get` and `set`, coupling/modflow6.py). Every link reads and sets MODFLOW 6's
  memory through them, so each value crosses in meters, cubic meters, days and kilograms whatever the model holds:
  recharge (m d-1), the stream exchange's stages, bottoms, conductances and fluxes on RIV, DRN and SFR (A and B
  alike), users' wells, groundwater up into the soil, the cell areas the recharge map is checked against, heads,
  every package's flow and storage booked in the ledger, and a transport model's concentrations, its sources' and
  its mass budget (concentrations as kg per the model's volume unit). Each variable converts by its dimension, a
  (length, time) pair (`DIMENSION`); one name meaning two things (EVT's rate, a depth rate, and MAW's, a flow) is
  told apart by package type. Integers (counts, nodes, status) are never converted. In a model not in meters and
  days, a real value of a dimension not listed is refused rather than passed through; in meters and days nothing
  is converted, so models Basin Ledger builds run exactly as before. Simulated time and the one-day step are
  checked in days.
- **Maps** (`gridmaps.py`): the grid's world coordinates are in its length units, so a grid in feet takes polygons
  and lines projected to a CRS in feet; areas and lengths come out in square meters and meters, by the grid's
  `LENGTH_UNITS` or the `length_units` a study gives where it is undefined.
- **The start-up check reports it.** `basinledger check` prints the model's units and every factor (lengths ×
  0.3048, areas × 0.0929, flows × 0.0283 to m3 d-1, conductances, depth rates, concentrations), with which units
  were stated in the config; a run in other units than meters and days starts with the same line as a warning.
  `Basin.check(..., units=[])` returns it.
- **What isn't converted, because it is the model's:** its solver settings (IMS's closure criteria and the
  backtracking limits of its moderate and complex settings are absolute, in its units, and its builder set them),
  the files MODFLOW 6 writes (heads and budgets in its units), and SFR's Manning's constant (checked for its units,
  above). MODFLOW 6 calibration parameters multiply the model's values, so they carry no units.
- **Tested on twins** (`tests/test_units.py`): the same model in meters and days and in feet and days, feet and
  seconds, or centimeters and hours, every length, area, rate, conductance and tolerance converted, so MODFLOW 6
  solves one problem. MODFLOW 6 alone (recharge, a well and stages set; heads, every package's flows and storage
  read), SFR routed and at SWAT+'s stage, solute mass with sorption, maps on DIS and DISV, and a restart's remaining
  simulation all give the meters twin's values to 1e-6. The synthetic foreign model coupled to Osu_1hru for 75 days
  in feet and days, and in feet and seconds left undefined and stated in the config, gives its meters twin's ledger
  to 1e-6 and closes across each converted link (recharge received = made, each well delivered = requested, the
  model's own well read back at its rate). The twins take IMS's simple settings, with the residual criterion a
  tenth as tight in seconds: the solver can't reach 1e-8 m3 d-1 when it is 4e-12 ft3 s-1.
- **Waiting until a study needs them** (Kalamazoo's points 2 to 4): storage for daily steps, an inset of a model
  too large to run, and MAW wells as users' wells.

**Built (Nov 14, 2026): daily steps, each stress period kept** (Andrew: yes to building it next).

- **Only TDIS changes, in the run's copy** (`coupling/timesteps.py`, `daily_steps`). Each stress period keeps its
  place and length and is cut into one-day steps, its growth factor 1. The number and order of periods are
  unchanged, so every package's period blocks, time series, observations and output control apply to the same
  periods as before. Each rewritten row keeps the original after it (`was ...`), which MODFLOW 6 reads past, and
  the folder a study handed over is never touched. A model in daily steps already is left as it is; Basin Ledger's
  own builders write daily steps, so its own models are unchanged.
- **Periods that aren't whole days** (30.4375-day months, 36 hours) end on the nearest day, each end measured from
  the simulation's start so nothing builds up: no end moves more than 12 hours, and the start-up says how many moved
  and how far. A period too short to keep a whole day of its own is refused. Time units are the model's, or stated
  in `[modflow6] time_units` where TDIS leaves them undefined, as for the conversion.
- **Adaptive time stepping** (ATS6) is commented out: every step is a day.
- **Steady periods** are cut into days like the others: each day is solved without storage, with that day's
  recharge and exchange. A long steady period becomes many steady days; a study that wants it as a start instead
  runs it apart and starts from its heads.
- **Output control** keeps its words: `ALL`, `FIRST` and `LAST` mean what they did, `STEPS` and `FREQUENCY` count
  the new daily steps, and a model that saved every step now writes every day (Pleasant Lake: 2,558 days of heads
  and budgets instead of 421 steps).
- **The calendar.** MODFLOW 6's first day is the run's first day. Where TDIS has a `START_DATE_TIME` that differs
  from SWAT+'s start, the start-up says so and by how many days, since monthly periods carry their dates' inputs.
- **The start-up reports it.** `basinledger check`, and a run's first warnings, say what was rewritten: the periods,
  the days, the steps they had and whether they grew.
- **Restarts across stress periods** (`timesteps.remaining_days`, under Restarts above): a copied model's days left
  from any day, its periods renumbered and its time series moved back.
- **Tested** (`tests/test_timesteps.py`, `tests/test_existing.py`): on TDIS files alone, monthly periods of one
  step and of five growing by 1.2, daily steps left untouched, 365.25-day years in equal months, a model in hours,
  undefined, stated and years' units, ATS, the start date, a one-period restart, and Pleasant Lake's 85 periods (its
  2,558 days, each period its own length). With libmf6, the synthetic foreign model in five growing steps a month,
  rewritten, gives its daily twin's heads to the last digit for 59 days. Coupled to Osu_1hru for 75 days, the same
  model gives the daily model's ledger to 1e-12; handed over as a folder dated two years after the run, the check
  rewrites its steps, notes the 730 days, and leaves the folder as it was.

**Built (Nov 14, 2026): the mover in the line-up rule** (Pleasant Lake's second need).

- **Read from its files** (`coupling/mover.py`). libmf6 holds each mover's types and values (`<model>/MVR` for a
  model's own, `GWF-GWF_<n>/MVR` for an exchange's) but not which entries a move joins, so Basin Ledger reads the
  mover's file: a model's own from its name file's `MVR6`, an exchange's from the exchange's options; with or without
  `MODELNAMES`; each period block in effect until the next, as MODFLOW 6 reads them.
- **A move is a connection to the rule.** A move from one of the coupled SFR package's reaches to another, or into
  one of its reaches from another SFR package or another model's (taken as the model's own reaches), is checked
  as SFR's connections are: allowed from a channel's last reach (or the model's own) into a channel's top reach, a
  channel problem anywhere else (B runs it). Allowed moves are cut each day after the day's period data are read,
  their value set to move nothing: 0 for FACTOR and UPTO, 10^30 for EXCESS and THRESHOLD. The water stays with the
  provider, down its own connections or out of SFR as its external outflow, which SFR's budget books; SWAT+
  carries it on, as for a cut connection.
- **Refused by name**, with the mover and the period block's move: a move within a channel (not built); from a
  coupled reach to a reach no channel has, to another package or to another model (SWAT+ routes the channel's
  water); and into a channel from a lake's outlet or a well (water SWAT+'s channel doesn't see: a lake's outlet
  water isn't groundwater, so it would come into SWAT+'s channel from outside both models; kept refused, Nov 15,
  2026: no model in hand needs it. A lake SWAT+ has a reservoir for is coupled with `[lakes]`, its outlets cut;
  otherwise the move is left out of the run's copy, the outlet's water leaving the model). Moves among the model's own entries stay the model's: a drain's water moved into the model's own
  reaches is the model's boundary flow, as before.
- **A drain's move into a coupled channel** (built Nov 15, 2026) is cut each day and the drain's water is that
  channel's the next day, as the exchange's drains' seepage is (`exchange:<channel>:drains_moved`, part of
  `to_stream`; the drain's entry Basin Ledger's in the takeover): it leaves the aquifer either way, and SWAT+ carries
  the channel's water. A drain the exchange's drains map to another channel is refused.
- **Pleasant Lake:** the exchange's two moves carry the inset's outlet reaches into two of the parent's top reaches.
  Coupling the parent, they line up and are cut; coupling the inset, they take its coupled streams to the parent and
  are refused, until both models are coupled together (need 3).
- **Tested** (`tests/test_mover.py`, `tests/test_central_sands.py`): the line-up network with its tributary joined
  through the mover instead of SFR's connections. Lined up, the move is cut and the main stem carries only its own
  water while the tributary's leaves SFR at its last reach; the model's own tributary moved into a top reach is cut
  and keeps its inflow; each of the four types is cut each day across period blocks; a move into the middle of a
  channel is refused for A and runs as B; and the refusals and the drain's move into the model's own reaches. On
  Pleasant Lake, the mover read from the exchange's file and cut through it.

**Built (Nov 14, 2026): restarts across stress periods** (What's next, 4).

- **The days left, as text** (`timesteps.remaining_days`, called by `remaining_simulation` on the run's copy after
  its steps are rewritten to days). TDIS keeps the periods not yet over, the one the restart day falls in shortened
  to its days left. Every input file the simulation reads (each model's packages, each exchange, and the movers,
  observations, time series and time-array series their options read) has its PERIOD blocks renumbered: blocks from
  that period on move down, the last block before it becomes period 1 where that period has none of its own (its
  data are in effect that day), and earlier blocks are left out. Time series and time-array series move back by the
  days run, so each value falls on the same day. Nothing else in the files changes.
- **Tested:** renumbering for each place a restart can fall (`tests/test_timesteps.py`), time series and
  time-array series moved, and the synthetic foreign model cut at February 10th (19 and 16 days left, its packages
  renumbered, flopy reading the result). MODFLOW 6 alone, the model in five growing steps a month with a city
  well on a time series, restarted at day 40 from the unbroken run's heads, gives the unbroken run's heads and well
  rate every day to the last digit. Coupled to Osu_1hru, handed over as a folder and saved on February 10th, the
  restarted run's ledger and heads are the unbroken run's to the last digit (`tests/test_existing.py`).
- **Still to save:** the state of packages that hold water from day to day (MAW, UZF, LAK, SFR with storage); a
  model with them restarts with their files' starting state, and the restart says so for each such package.

**Built (Nov 14, 2026): one model of a simulation with several.** A simulation can hold more than one groundwater
model, joined by GWF-GWF exchanges and solved together (an LGR inset and its parent, or neighbors). Basin Ledger
couples one (`[modflow6] model`). What crosses its GWF-GWF exchanges to the others is now booked as water from or to
outside, `groundwater:boundary:gwf-gwf_<n>` (`Modflow6.exchanges` and `exchange_flow`, from libmf6's SIMVALS, positive
into the exchange's first model), so the coupled model's daily closure holds; before, a coupled run of such a model
would have stopped on its closure. The start-up check notes each exchange and the model on its other side. Tested
on the synthetic foreign model with a neighbor east of it (`neighbor=True`), coupled 75 days: the ledger closes
every day and the exchange brings the regional water in (`tests/test_existing.py`).

**Built (Nov 15, 2026), on the config recommended and chosen by Andrew: several groundwater models coupled
together** (Pleasant Lake's need 3). libmf6 already solves every model of a simulation together; Basin Ledger now
couples more than one of them.

- **The config.** `[modflow6] models = ["kzoo", "kzoo_east"]` (the first is the one libmf6 is started for; `model`
  and `models` together are refused). Maps stay one file each, with a `model` column naming each row's model (a row
  without one is the first model's); a well takes `model` beside its `cell`, as do a lake, the exchange
  (`[exchange] model`) and a head output. The other way, a `[modflow6.<model>]` table per model with its own maps,
  repeats more and checks less (the area check needs every model's cells at once); Andrew kept the list (Nov 15).
- **The cells in one sequence** (`modflow6.Models`): each model has a view through the one libmf6 (`Modflow6.view`:
  its own grid and packages, the simulation's time and units), and their cells are numbered one model after the
  other. A map, the recharge link and the wells take a cell of any of them without knowing there are several: a
  package named in several models (the recharge package, `rch_swb` in each) is one package whose entries are each
  model's in turn. Models in different length units are refused (the conversion is one per run).
- **Recharge.** One map over every model; each HRU's or aquifer's area is checked against all its cells at once, so an
  overlap where an inset's parent cells aren't inactive is caught. The ledger adds each model's share
  (`recharge:received:<model>`).
- **The ledger.** Each model's own accounts under `<model>/` (`groundwater:boundary:kzoo_east/chd_regional`), its
  storage and its closure (`groundwater:<model>:closure_share`), checked every day with the whole's. A GWF-GWF
  exchange between two coupled models is booked once, inside the groundwater (`groundwater:exchange:gwf-gwf_<n>`,
  positive into its first model), and cancels in the whole's closure; one to a model not coupled stays a boundary,
  as before. The start-up check lists each model's packages as `<model>/<package>` and notes each exchange.
- **Streams in several models.** An exchange map whose rows name more than one model runs the streams through each
  model's package (`exchange.StreamExchanges`): one exchange per model, each with its own channels (a channel's
  reaches are in one model, or it is refused), the day's inflow, supply limit and exchange asked of each and put
  together by channel. The SFR line-up rule takes the other models' maps (`sfr_problems(..., others)`): a coupled
  reach giving water through the mover to another model's coupled reach lines up as a connection between channels
  does, from its channel's last reach into the other's top reach, and is cut; into the same channel there it is a
  move within a channel, refused. On Pleasant Lake, each model's SFR to its own lines with a SWAT+ outlet where the
  inset meets its parent: the mover from the inset's two outlet reaches into the parent's top reaches lines up on both
  sides and is cut; the lines as SFRmaker numbered them, the same line on both sides of the edge, are refused as a
  move within a channel (`tests/test_central_sands.py`).
- **Drains and lakes at a fixed stage** (built Nov 15, 2026): each model's own (`[exchange] drains` or `lakes`, a
  package named alike in each), mapped onto its own channels by the map's rows (`exchange.part_map`); a drain or
  lake whose channel's reaches are in another model is refused.
- **One model's links.** A lake, the head outputs and an exchange in one model name their model (`[exchange]
  model`).
- **Across the models** (built Nov 15, 2026): `[recharge] add` adds a package in each model the map has cells in;
  `[rivers]` sets the stage of each model's river package of that name on its mapped cells; groundwater up into the
  soil takes the soil package in each model with mapped cells, set and read across them (`Models.node_values`,
  `entry_values`, `set_entries`; `Models.set` without a model refuses a package that's in several). Transport is
  refused with several models: it needs a transport model per groundwater model and species, the exchanges between
  them and a mass budget per model.
- **Restarts** save every model's heads and set them back.
- **Tested** (`tests/test_models.py`) on the synthetic foreign model and its neighbor, both coupled: Osu_1hru's
  10 ha aquifer recharging 64 of kzoo's top cells and 36 of its neighbor's through one map (each model's share as the
  map's, 36:64, every day), a second plant on a well added in the neighbor, the exchange between them inside the
  groundwater, each model's closure and the whole's holding over 75 days; a restart continuing the unbroken run to
  the last digit; the models' cells in one sequence; the settings read and refused.
  With streams in both models, MODFLOW 6 alone: the synthetic creek in kzoo and a losing river in its neighbor, one
  exchange through both, the supply limit holding the neighbor's river to the water its channel carries; with each
  model's tile drains and the neighbor's lake, each model's seepage and lake flow back to its own channel. Coupled
  with SWAT+: a recharge package added in both models (75 days), both river stages from SWAT+'s depth (30 days), and
  groundwater up into the soil from both (75 days), every closure holding.
- **Next for Pleasant Lake:** the study's SWAT+ build, with an outlet where the inset meets its parent and the lake as
  a reservoir.

**Built (Nov 14, 2026): wells the model doesn't have** (a new large user's, the question Basin Ledger exists for).
A `[[wells]]` entry with `add = true` is a well the model lacks: Basin Ledger writes a WEL package of its own, named
by the well's `package`, into the run's copy (`existing.add_wells`), one entry per such well at its cell with a rate
of 0 in the file, and lists it in the model's name file; the users' routes set its rate every day, as for the
model's own wells. The package is all Basin Ledger's in the start-up check's table, and the check notes each well
added. Refused: a package name the model already has (a well added never mixes with the model's own entries), a cell
in the wrong form for the grid, and two added wells at one cell of a package. The folder handed over is untouched; a
restart re-keys the added package with the rest. Tested on the synthetic foreign model (`tests/test_existing.py`): a
second plant drawing from a well added in the bottom layer, coupled 75 days, gets what it asked for every day and the
ledger closes; the refusals. MAW wells as users' wells are built too (Kalamazoo's point 4, above).

The same for recharge: `[recharge] add = true` gives a model without a recharge package one of Basin Ledger's own
(`existing.add_recharge`), an RCH list package named by `[recharge] package` with one entry per cell of the recharge
map (a CSV of cells), its rates set every day from SWAT+. A model that has a recharge package of its own is refused
at the takeover, as before (recharge would be counted twice); a map given as `all:` or by nodes is refused (it names
no cells to put the package on). Tested on the synthetic foreign model with its recharge taken out of its name file:
the check lists the package as Basin Ledger's, and coupled 75 days the recharge received is the recharge made and
the ledger closes.

**Plan.** A design pass on the open questions (done above); then the start-up check and maps on a foreign grid, tested
on the synthetic foreign model and then a published model; then stress periods and restarts; then solutes. Each a
minor release with tests.

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
- **Studies** live in their own repositories, one per basin, never here. A client study is private. A study:
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

Done for 0.2.0 (Nov 5 to 7, 2026), each tested on the synthetic examples, and all but transport run on a study's
coupled model:

- `CoolingWaterBmi`'s limit and dry-mode inputs (waterledger-models 3.11), with the rules on a real cooling user;
- intake concentrations on users, and a return's heat at its point source;
- stage from SWAT+'s channel depth (swatplus-bmi 0.2.0);
- conservation across links, the whole-basin balance, capture, stream temperature, groundwater transport;
- scenarios from one saved warm-up (`scenarios.py`; a study runs its cases in parallel processes), with the warm-up's
  soil check.

Done for 0.3.0 (Nov 8, 2026): recharge HRU by HRU on swatplus-bmi 0.6.0's HRU percolation (Recharge from HRUs), tested
on Osu_1hru and Ames_sub1 and run on a study's coupled model.

Done for 0.4.0 (Nov 9, 2026), on swatplus-bmi 0.7.0: weather into SWAT+'s gage records from a run's first day (Weather
into SWAT+), and groundwater up into the soil, with the drain and no cap chosen (Groundwater up into the soil). Not
built: the start-up check refusing EVT or UZF alongside groundwater into the soil.

Done for 0.5.0 (Nov 10, 2026): MODFLOW 6 vertex grids (DISV), cells as layer, icell2d throughout (Cells,
under NextGen geography), for a study's quadtree grid refined around its site.

Done for 0.6.0 (Nov 11, 2026): the two-way exchange through SFR (Two-way: stream-aquifer exchange);
calibration parameters written into SWAT+'s and MODFLOW 6's files, and PEST++ IES files (Calibration).

Done for 0.7.0 (Nov 11, 2026): SFR reaches built from SWAT+'s channel lines (Through SFR), and
`basinledger check` (Conservation across links, checked at start-up).

Done for 0.8.0 (Nov 15, 2026): extending an existing MODFLOW 6 model (the second way in, next item 7 below): the
start-up's takeover, maps onto a grid Basin Ledger didn't write, a model's own units and stress periods, the SFR
line-up rule and the mover, restarts of a model with several stress periods, wells and recharge and storage added
to the run's copy, MAW wells, several models coupled together, LAK lakes on SWAT+'s reservoirs and irrigation from
wells (swatplus-bmi 0.11.0); and lakes at a fixed stage through their outlet channels, check's start-up hooks and the
recharge precision bound, what the Clam basin study ran on.

Next:

1. Water quality's stage 2: loads through SWAT+'s streams, on swatplus-bmi 0.5.0's species variables (Water quality
   across links).
2. Point source loads (asked of swatplus-bmi).
3. Two-way refinements: leaks and infiltration in the chain exchanged with MODFLOW 6; SFT for transport with an SFR
   exchange.
4. Restarts for MODFLOW 6 models with several stress periods (built Nov 14, 2026), and the advanced packages' own
   state, LAK's and MAW's (built Nov 15; UZF's and CSUB's still to save).
5. MODFLOW 6 parameters by pilot point (built Nov 15, 2026: ordinary kriging in numpy, Calibration); USGS groundwater
   levels as observations.
6. AORC forcing for SWAT+ and the users; the SWAT+ project builder on the hydrofabric, with outputs at every nexus.
7. Extending an existing MODFLOW 6 model (decided Nov 14, 2026; the section above): the start-up check and maps on a
   model Basin Ledger didn't write, then stress periods and restarts, then solutes.
8. Plaster Creek, in the calibration order above:
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
