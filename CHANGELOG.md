# Changelog

## Unreleased

- SWAT+'s own irrigation from its channels and reservoirs, booked (docs/design.md, "Irrigation from surface water";
  `coupling/irrigation.py`, `NativeIrrigation`; swatplus-bmi 0.14.0). Where a project's decision tables irrigate
  from a channel or reservoir, the driver reads each day what each water body gave and what each HRU took, and books
  those amounts, what reached the field and the application's loss. Two pairs close each day.
  - The whole-basin balance's channels and lakes gain `irrigation_from_streams` and `irrigation_from_lakes`, and the
    balance reports `irrigation_loss`.
  - With `[irrigation] wells`, `irrigation:from_outside` no longer counts water SWAT+ took from its own channels and
    reservoirs. That was wrong before.
  - On an older swatplus-bmi, such a project is warned that the irrigation isn't booked, or refused where
    `[outputs] balance` is asked for. SWAT+ 62.0.1 never took a channel's irrigation out of the channel, and
    swatplus-bmi 0.14.0 patches that.
  - Tested on Osu_1hru's paddy irrigating from its channel, coupled both ways for three years. The ledger's pairs
    closed every day, and the channels closed in the whole-basin balance with the 1,300 m3 they gave.
- From 0.10.2: places of use on more than one channel, several to a channel served senior first. These need
  swatplus-bmi 0.14.0, and an older one is refused them (`tests/test_intakes.py`).

## 0.11.0 (2026-11-23)

- Groundwater transport, stage 2 (docs/design.md, Water quality across links): nitrate and soluble phosphorus travel
  with the water through every link between SWAT+ and MODFLOW 6, each link's mass in the ledger as sent and received
  and every day's mass closing (`coupling/surface_mass.py`; `[transport.<species>] surface = "no3"` or `"solp"`,
  `recharge = "swatplus"`, `stream = "swatplus"`; `examples/synthetic/nitrate.toml`):
  - recharge carries SWAT+'s own leached mass HRU by HRU (`hru_pw:percn`, `hru_ls:lchlabp`, swatplus-bmi's reads of
    SWAT+'s daily records), on the cells the HRU's water goes to, in the same shares; nitrate's is held to SWAT+'s own
    aquifer recharge mass every day. Per aquifer only with `[recharge] source = "aquifer"`, with a warning; refused
    for soluble P, which SWAT+ 62.0.1 leaves out of the aquifer's recharge;
  - the aquifer's discharge to a channel (rivers' gains, drains) carries its cells' concentrations into SWAT+ the next
    day with the gain; a loss takes SWAT+'s own outflow concentration, and its mass reaches the loss cells the day
    after its water, through a mass source on the cells (`add_transport(..., stream_mass=[...])`), with whatever SWAT+
    didn't keep of the aquifer's discharge, so the link closes exactly;
  - the soil link's water carries its cells' concentrations into SWAT+'s soil layers, split as the water is;
  - lakes at a fixed stage (`[exchange] lakes`) are their outlet channel's, as their water is (still refused with a
    stated stream concentration);
  - SWAT+ makes its daily records from the first year in a run's copy (`print.prt` nyskip 0) when a species is
    carried, since those reads are NaN in skipped years.
- First-order decay in the transport model (`add_transport(..., decay=, decay_sorbed=)`, d-1 [A], denitrification for
  nitrate), no longer refused: the mass decayed is booked (`transport:<sp>:decay`) and the day still closes. Linear
  sorption's release from the sorbed store is booked on its own (`sorbed_release`); Kd may be given by zone.
- `Basin.transport_steady_state`: each species' steady state under a run's mean flows and inflows (advection, mixing
  and decay; dispersion left out [A]), to start a run from (`Basin(..., concentrations=...)`).
- Fixed: SWAT+'s aquifers sent their channels soluble P in two-way mode: SWAT+ 62.0.1 sends the day's recharge volume
  times the aquifer's `sol_p` with no water when return flow is off. `aquifer_pass_through` sets `sol_p` to 0, and
  a project or calibration that sets it is refused.
- The transport solve's residual criterion is 1e-8 kg d-1 (was 1e-10): at 1e-10 round-off held the residual above it
  on Plaster Creek's 13,215 cells (day 47) and on the synthetic grid with sorption. The ledger's closure check is
  unchanged.
- Every MODFLOW 6 solve is checked for non-finite heads, after each outer iteration of a day with a coupling's callback
  and after a day solved whole, and stops at the first, naming the day and iteration (`Modflow6._check_finite`). A
  linear solver that breaks down can return NaN for every head and still report the iteration as converged: in a
  Rappahannock calibration, COMPLEXITY MODERATE's modified ILU (RELAXATION_FACTOR 0.97) did so for 16 of 40 prior
  realizations, and the NaN surfaced only a day later, in the groundwater-to-soil link. The message names the IMS
  settings that usually solve it (RELAXATION_FACTOR 0, or COMPLEXITY COMPLEX).
- Replay (`replay.py`, `[replay] record`, `Basin(cfg, replay=folder)`, `--replay` on `run` and `forward`;
  docs/design.md "Replay"): a coupled run records, each day, what transport and the salt balances read from the flow
  side (MODFLOW 6's face flows, saturation, specific discharge, storage, dry cells and every package's flows by
  entry; SWAT+'s columns the salt reads; the users' lines and controls), and a replay runs only the transport models
  and the salt from it, with MODFLOW 6's flow model in the simulation but unsolved, its memory set from the record
  each day. From a float64 record the coupled run is reproduced bit for bit; from float32 (the default) to 1.1e-7.
  Salt settings, salt parameters, porosity, dispersion, weathering, advection and solver may change; anything that
  would change the flow is refused (the config's flow side, the files it names and the flow model's input, hashed).
  Users who take in the aquifer's water are run again, their water checked against the record. Record size on the
  Virginia study's Rappahannock grid (335,691 cells): 20.0 bytes per cell-day in float32, 14.1 with
  `mantissa_bits = 16`, 12.1 and 8.6 with `specific_discharge = false`; 11.6 to 27 GB for 11 years, against about
  80 GB for MODFLOW 6's daily budget file. A replayed day there: 3.64 s with XT3D, 0.37 s without, 0.24 s with no
  dispersion, against 7.68, 4.52 and 4.39 s for MODFLOW 6's flow and transport, before SWAT+.
- `add_transport(xt3d=False)` turns XT3D off; `dispersivity=None` with no diffusion leaves DSP out (on the Virginia
  grid XT3D was most of the transport's time).
- Observations of salt: kinds `<species>_mg_l` and `conductance_us_cm` at an output channel.
- Weathering by each cell's water: `[salt]` scales each weathering entry (`salt.add_weathering`, written at full
  saturation) every day by its cell's saturated share at the start of the day (MODFLOW 6's saturation at the end of
  the day before; on the first day from the starting heads, `salt.start_saturation`), so a dry cell weathers nothing
  and a partly saturated one in proportion. The static rate put the full source into the Virginia study's dry
  layer-1 cells and its weathering solve diverged on day 0. The saturation is saved with the salt state, so a
  restart carries it. `build_small` takes `upper` (a thin convertible top layer), `newton` and `strt` for the test.
- Salt (`salt.py`, `[salt.<species>]`, docs/design.md "Salt"): chloride carried through every water path of a coupled
  SWAT+ / MODFLOW 6 run, each day, as a conservative mass, since SWAT+ 62.0.1's salt module can't run with aquifers
  or point sources. A well-mixed store per HRU (deposition in precipitation, road salt on snow and ice days or a
  stated series, groundwater up into the soil; out with runoff, lateral and tile flow and percolation, read through
  swatplus-bmi's printed daily columns), optionally behind a snowpack and road store that sends its melt's chloride
  to the stream fast; GWT's recharge concentration each cell's share of the HRUs' percolation mass, placed by the
  recharge map; weathering as a GWT mass source by zone (`salt.add_weathering`); and a mass balance for every SWAT+
  channel and reservoir, upstream first, through SWAT+'s routing units (the land, rain, point sources at a user's
  returned mass or a stated concentration, the aquifer's gains a day after GWT gives them, withdrawals, seepage and
  the outlet). The streams' losses enter the aquifer with no chloride and their mass follows the next day, at the
  channel's mixed concentration, through a GWT source on the river and lake cells (`salt.add_stream_sources`), so
  the stream gives up and the aquifer takes the same mass: setting the channel's concentration on the river package
  instead left the channel short of 46% of the losses' mass the next day in the synthetic run. The land, the streams
  and the whole system (with GWT's storage) close every day, to about 1e-12. Specific conductance from chloride by a
  study's relation, per output channel (`series.csv`). Lakes at a fixed stage (`[exchange] lakes`) now take
  transport when salt carries the species (a lake loses at its outlet channel's concentration). Refused at start-up:
  `[[lakes]]` (LAK), `[irrigation]`, SFR, several groundwater models, NextGen, a model not in meters and days.
  Example `examples/synthetic/salt.toml`; its builder (`build_small`) takes `weathering` and `stream_sources`.
  Cost on the synthetic example: 2.2-2.7 ms a simulated day on 8.4-8.6 ms (one HRU); the land and stream kernels
  take 5 ms a day at the Virginia study's size (33,525 HRUs, 1,192 channels and reservoirs), before the BMI reads.
- Salt runs any number of conservative species together, each with its own GWT model, ledger and closure, sharing
  the day's SWAT+ reads: weathering's ions (`[salt.weathering]`, no rain or road salt, from GWT's weathering source
  by zone) beside chloride. Specific conductance is a constant plus each species' slope times its mg/L
  (`[salt.conductance] base, slopes`), replacing chloride's own `conductance`.
- Land sources (`sources = [{name, rate, by, timing, group}]`): kg per m2 of the HRU a year as a number, by HRU, by
  SWAT+ land use or by calendar year (`hru, year, kg_m2`), timed evenly, by 12 monthly weights or on the HRU's first
  fertilizer day of each year (swatplus-bmi's `hru_nb:fertn`, `fertp`), each with its own ledger line.
- Road salt by an annual amount per calendar year (`road_salt = {annual = "hru, year, kg_m2" CSV}`, the USGS grid's
  form), spread over that year's snow and ice days, worked out at start-up from the stations' weather records, in
  proportion to each day's precipitation (or equally); years outside the table take the mean of its last or first
  N years (`repeat_last`, `repeat_first`).
- Multipliers on each input group's mass (`multipliers`; calibration parameters of model `salt`, target
  `<species>:<group>`), applied where the mass enters, so the ledger books the scaled mass.
- A point source's water from a user that doesn't report the species takes the stated concentration (it had counted
  as carrying none). The synthetic example runs with the surface store on. Cost on it with two species and the land
  sources: 15 ms a simulated day.

## 0.10.1 (2026-11-22)

- A map from pixels reads in seconds: `AreaMap.from_csv` finds each object by name in a table and each distinct cell's
  node once (a 1,025,254-row HRU map took 173 s, most of it in list lookups), as does the soil link's layer index.
- The supply limit bisects all its channels together (`exchange.limit_to_supply`): the same depths as one channel at a
  time, to 1e-12 m (tested against 0.10.0's loop), at a seventh of the cost on 1,157 channels.
- Fixed: a day whose boundary a coupling kept moving (the supply limit) could run past the solver's `outer_maximum`,
  where MODFLOW 6 fails with a segmentation fault (0.10.0 did on the synthetic example's 25 at a tight stage
  tolerance), or end on a moved boundary it hadn't solved, which the ledger refuses. The boundary now moves for at most
  half of `outer_maximum` (read from MODFLOW 6) or 30 rounds after convergence, and the rest of the day solves it as
  it stands; a day that still doesn't converge stops with an error naming `outer_maximum`.
- `[exchange] limit_tolerance` (m, default 0.001 as before): the supply limit's stages agree with the heads once none
  moves more than this in an outer iteration. On a 1,157-channel basin the stages chased the heads for 43 outer
  iterations a day at 1 mm; at 1 cm it took 11, with gage flows within 1e-4 of 1 mm's, but one reservoir's lake lost
  60% more over Jan-Sep 1991 (a lake's conductance is large), so check channels and lakes, not only gages, before
  loosening it.
- `[exchange] limit_warm_start` (default on): each day's supply limit starts its search from yesterday's lowered
  stages. The answer is the same (the limit works from the day's own stages and puts a start that is too low back
  up; tested), found in fewer outer iterations: on the same basin at 1 mm, 12.9 a day instead of 53.6 (35-54 in
  summer), with gage flows within 2e-5 and that lake within 1% of the cold start's. Where each channel's search starts
  is saved with a restart (`limit_lowered` in basin.json), so restarts still match the unbroken run bit for bit.

## 0.10.0 (2026-11-16)

- `pip install basinledger` brings both models on Windows x86-64, Linux x86-64 (glibc 2.34 or later) and macOS 15 or
  later on Apple silicon: Basin Ledger's wheel for each carries USGS's libmf6 (MODFLOW 6 6.8.2, unchanged, with a note
  on its source and license), and swatplus-bmi 0.13.1, now the minimum, carries SWAT+. Each platform wheel is built
  and run installed alone before it's published (`hatch_build.py`, `tools/lib_tag.py`, `tools/check_wheel.py`).
  libmf6 is found in the wheel after a given path and `$LIBMF6`, before `get-modflow`'s folder. Elsewhere, both come
  separately, as before.

## 0.9.1 (2026-11-16)

- Windows: tested on every push (windows-latest; SWAT+ built with MSYS2's gfortran by `swatplus-bmi build`, which needs
  swatplus-bmi 0.12.0, now the minimum). libmf6 is found where `get-modflow` puts it there (`%LOCALAPPDATA%\flopy\bin`), and the
  example's fetch asks it for `mf6.exe`.

## 0.9.0 (2026-11-15)

- Several coupled models, more links: drains and lakes at a fixed stage with streams in several models, each model's
  on its own channels (`exchange.part_map`); `[recharge] add` in each model the map has cells in; `[rivers]` setting
  each model's river stages; groundwater up into the soil from each model's soil package (`Models.node_values`,
  `entry_values`, `set_entries`). `Models.set` without a model refuses a package that's in several. Transport stays
  refused with several models.
- A drain's move into a coupled channel's reach is cut each day and the drain's water is that channel's
  (`exchange:<channel>:drains_moved`, part of `to_stream`), in place of refusing it. A lake's outlet or a well moving
  water into a channel stays refused.
- Restarts carry the advanced packages' state: lakes' stages (LAK) and multi-aquifer wells' heads and wellbore
  storage (MAW), saved beside the heads (`packages.npz`). UZF and CSUB still start from their files, as noted.
- MODFLOW 6 parameters by pilot point: `<package>.<array>:pp=<set>`, the points' locations in `[calibration]
  pilot_points`' file, their values kriged over each layer (ordinary kriging, exponential covariance;
  `calibration.kriging_weights`), in log10 for log-transformed sets. Each point is a parameter row, so the PEST++
  files take them as they are.
- A failed start-up finalizes the models it started.

## 0.8.1 (2026-11-15)

- A recharge package left as the model's is refused only where it recharges ground SWAT+'s recharge also covers: on
  each cell, the package's own share (its auxiliary multiplier, or 1) beyond what SWAT+'s recharge map leaves of the
  cell, by more than 0.1% of the cell, is recharged twice (`existing.recharge_cover`). Off SWAT+'s recharge, or
  standing for the rest of a cell SWAT+ covers in part, it is the model's, booked as water from outside, and noted.
  0.8.0 refused any such package, which stopped the Clam study's recharge below its gage.

## 0.8.0 (2026-11-15)

Needs swatplus-bmi 0.11.0 (reservoirs and HRU irrigation).


- Groundwater up into the soil splits each day's water among HRUs and layers for all map rows at once (`gwsoil.to_layers`, vectorized; the row-by-row version stays as `to_layers_loop`, the reference it is tested against to 1e-12). On a study basin (2,371 HRUs, 11,400 rows) the loop was half of a coupled day; a year went from 103 s to 68 s with the same results.

- The daily pair of percolation sent and SWAT+'s own aquifer recharge (recharge by HRU) is bounded by the precision the project's files give each HRU's area and share (`hru.con`, `rout_unit.ele`), added to its tolerance each day (`recharge:precision`; a ledger pair's `bound`). On a study basin with groundwater up into the soil, a frozen day's 4.76 m3 of recharge, carried by a few small HRUs whose shares are written to 7 decimals, missed the fixed 1e-5 by 9e-5 m3 and stopped the run.

- `Basin.check` takes the run's start-up hooks (`prepare`, `parameters`, `wrap_surface`, `heads`), so a study checks
  what it runs: a prepare hook that rewrites SWAT+'s `calibration.cal` is checked on the file it writes, not on the
  project's stored one.
- Lakes at a fixed stage coupled to channels (`[exchange] lakes = {package, map}`): a RIV package of the lakes' own,
  each lake cell mapped to its outlet channel, one of the exchange's coupled channels. The lakes' net exchange joins
  that channel's the next day (`exchange:<channel>:lakes`), and the package's coupled entries aren't a boundary
  flow. A lake takes no more than its channel has: with RIV the supply limit lowers a channel's river and lake
  stages together; with SFR the lakes have a limit of their own, and a combined excess is a booked shortfall
  (`exchange:<channel>:lake_stage_lowered_m`). Refused at start-up: a package that isn't RIV, the exchange's or the
  drains' package, `[rivers]` on it, a channel the exchange doesn't couple, and transport with lakes. The synthetic
  examples `lakes.toml` and `lakes_sfr.toml`; restarts match the unbroken run bit for bit.

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
- The mover in the SFR line-up rule: movers read from their files (`coupling/mover.py`, a model's own or an
  exchange's), each move between SFR reaches checked as a connection and, lined up, cut each day (its value set to
  move nothing, for each type). A move within a channel, out of a coupled channel to anything but another channel's
  top reach, or into a channel from a drain, lake or well is refused by name; moves among the model's own entries
  stay the model's. Pleasant Lake's exchange mover lines up on the parent's side.
- Restarts of a copied MODFLOW 6 model with several stress periods (`timesteps.remaining_days`): the periods over
  left out, the current one shortened, every input file's period blocks renumbered and its time series moved back,
  all as text. Restarted runs of the synthetic foreign model, alone and coupled, match the unbroken run bit for bit.
- One groundwater model of a simulation with several coupled: what crosses its GWF-GWF exchanges to the others is
  booked as water from or to outside (`groundwater:boundary:gwf-gwf_<n>`), so its daily closure holds, and the
  start-up check notes each exchange. The synthetic foreign model takes a neighbor (`neighbor=True`). Coupling
  several models together is proposed in the design, not built.
- Wells the model doesn't have: `[[wells]] add = true` puts a user's well in a WEL package Basin Ledger adds to the
  run's copy (`existing.add_wells`), its rate set every day; a package name the model has, a cell in the wrong form
  and two added wells at one cell are refused.
- A recharge package for a model without one: `[recharge] add = true` adds an RCH package of Basin Ledger's own on
  the recharge map's cells to the run's copy (`existing.add_recharge`).
- Storage for a model handed over steady (`[modflow6] storage`, `existing.add_storage`): a model with no storage
  package gets one of Basin Ledger's own in the run's copy (SY, and SS per meter in the model's units), transient;
  `storage = {}` makes a model's own steady periods transient. A steady model without it is noted.
- Multi-aquifer wells as users' wells: `[[wells]] well = <n>` names a MAW well by its number in place of a cell;
  its rate is set each day and what it pumped delivered, booked from the aquifer (its connections) and from the
  wellbore (`QSTO`).
- Several groundwater models coupled together (`[modflow6] models`, `modflow6.Models`; the config chosen by Andrew,
  Nov 15): each model a view through the one libmf6, their cells in one sequence, so one recharge
  map (with a `model` column) and the wells take cells of any; each model's accounts, storage and closure in the
  ledger, and a GWF-GWF exchange between coupled models booked inside the groundwater. The stream exchange, a lake
  and head outputs name their model; transport, rivers, the soil link and an added recharge package are refused
  with several models.
- Streams in several coupled models (`exchange.StreamExchanges`): an exchange map whose rows name more than one model
  runs each model's channels through its own package, and the SFR line-up rule lines up the mover between models
  (a coupled channel's last reach into another model's coupled top reach) and cuts it. Pleasant Lake's inset and
  parent line up with an outlet at the inset's edge.
- Lakes (`coupling/lakes.py`, `[[lakes]]`; needs swatplus-bmi 0.11.0): a SWAT+ reservoir holds a MODFLOW 6 lake's
  water and LAK is its face to the aquifer, as SFR is a stream's. Each day the lake is set from the reservoir, A
  active within the day from the reservoir's starting volume with its day's water in and out, or B at a fixed stage,
  and its leakage is the reservoir's groundwater exchange the next day (`lake:<reservoir>:due`, `applied`,
  `shortfall`, `to_aquifer`). The stage is LAK's own at SWAT+'s volume: its table, or its connections as MODFLOW 6
  computes a lake's volume, inverted. A coupled lake's outlets are cut each day (SWAT+ routes its release) and
  checked; an outlet from a lake not coupled into a coupled one, and a reservoir with its own seepage on, are
  refused.
- Irrigation pumped from MODFLOW 6's wells (`coupling/irrigation.py`, `[irrigation] wells`; needs swatplus-bmi
  0.11.0, `land_irrigation_water__volume_flow_rate`): each day the water SWAT+'s irrigation put on each mapped HRU is
  pumped that day from the `[[wells]]` the map names, in their shares (`irrigation:pumped`,
  `irrigation:well:<name>`); irrigation on HRUs not mapped stays water from outside (`irrigation:from_outside`).
  Each HRU's irrigation sources are read from its land use's schedule and decision tables, and a mapped HRU that
  irrigates from a SWAT+ channel, reservoir or aquifer, or through water allocation, is refused.
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
