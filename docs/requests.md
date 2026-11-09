# Requests to other packages

Basin Ledger asks for changes to the packages it depends on, rather than working around them here. Each entry says
what, why, and how Basin Ledger copes until then. Written Nov 4, 2026.

## waterledger-models

1. **Declare pyyaml as a dependency (a bug in 3.9.0).** `waterledger_models/datacenter/bmi.py` imports `yaml`, but
   pyyaml isn't among the package's dependencies, so `from waterledger_models.datacenter.bmi import CoolingWaterBmi`
   fails in a clean environment. Until fixed, Basin Ledger's dev group installs pyyaml. *Status: fixed in 3.10.0
   (yaml is imported only when a YAML config file is read). Basin Ledger still installs pyyaml with `users`, since
   its examples give `CoolingWaterBmi` a YAML config and waterledger-models declares pyyaml only under its `bmi`
   extra.*
2. **A withdrawal limit input on every user component:** `water_user_water__withdrawal_volume_flux_limit` (m3 s-1;
   infinity means no limit). Curtailment and well drawdown rules set it before each day. `CoolingWaterBmi` should
   meet a limit the way a campus would: by running dry or with more assist, using more electricity. Its balance must
   still close, and the electricity must follow. Until then, rules can only be tested on the synthetic user.
   *Status (requests 2 to 5): released in waterledger-models 3.11.0 (Nov 6, 2026), which Basin Ledger 0.2.0
   requires; a study's rules case curtails a real `CoolingWaterBmi` campus. The design follows realistic
   operation:*
   - *The heat to remove is fixed, since the campus must hold its temperature. A limit caps the wet share of the heat
     rejected each step, and the heat balance is solved again.*
   - *Capped evaporative duty moves to dry coolers and chillers, at more electricity. A limit of 0 is dry mode.*
   - *Heat the installed dry and chiller capacity can't take is reported as heat not removed (IT load that would
     have to be shed), never dropped.*
   - *Energy and water balances close every hour.*
   - *Return load = intake concentration × withdrawal for cooling × (1 − drift's share of the dissolved mass). Drift's
     share is drift × circulation × (cycles − 1) / evaporation: about 1–16% across the model's ranges and 3–6% at central
     values, not the 0.1% first written here (drift is 0.005–0.02% of the circulating flow, but blowdown is small beside
     it). Evaporation carries no dissolved solids, so all of the intake's mass leaves in blowdown and drift; blowdown
     concentration is the intake's times the cycles of concentration (waterledger design 0006).*
   - *As built: towers plants, with no dry side, shed the load a limit leaves them and report it as heat not removed.
     Dry mode is `water_user__dry_mode_flag`. The component also reports `water_user_unreturned_water_<name>__mass_flux`
     (drift, leaks, scale), which with the return's mass adds to the intake's, and the common names as aliases
     (request 5).*
3. **A dry-mode input on `CoolingWaterBmi`**, for days a rule calls for no evaporative cooling. It may be the same
   mechanism as the limit, with a limit of zero for cooling water.
4. **Intake concentrations and return loads.** For quality balances, a user needs the concentration of what it takes
   in as an input, and should report the mass in its return. Inputs: `water_user_intake_water_<name>__mass_concentration`,
   in kg m-3. Outputs: `water_user_return_water_<name>__mass_flux`, in kg s-1. Dissolved solids, chloride and
   phosphorus come first. For cooling, blowdown carries the intake's mass, concentrated by the cycles of concentration.
   Return temperature (`water_user_return_water__temperature`, degC) is also needed, for heat.
5. **The common names, as aliases.** Basin Ledger maps any names by role in an example's config, so this is optional.
   Still, reporting `water_user_water__withdrawal_volume_flux` and the other common names alongside
   `data_center_cooling_water__*` would let every user model plug in with no mapping. Each new user model (mines,
   battery and EV plants, fabs, food processing, hydrogen) should report the roles in `docs/design.md`, Water users,
   and close the balance withdrawal + gain = consumption + return + change in storage every step.
6. **State functions** (`get_state`/`set_state`, after CSDMS bmi PR #112), for any user model that keeps state
   between steps (ponds, tailings), so scenarios can start from a shared warm-up.

## swatplus-bmi

*Status: in swatplus-bmi 0.2.0: request 2, channel depth (`channel_water__depth`), and request 4, a channel's groundwater exchange set directly (`channel_water__groundwater_volume_exchange_rate`). Basin Ledger uses both. In swatplus-bmi 0.6.0: request 1, HRU percolation (`soil_water__percolation_rate`), which Basin Ledger's `[recharge] source = "hru"` uses. Still to come: point source loads (3).*

1. **HRU percolation (recharge) per HRU**, for mapping recharge to cells by HRU area as Bailey et al. (2025) do. Today
   only aquifer recharge is exposed, and SWAT+'s aquifers are lumped per subbasin.
2. **Channel water depth (stage)**, for MODFLOW 6's river stage in the two-way exchange.
3. **Loads on point sources:** settable daily mass (and temperature, if SWAT+ 62 carries it) on recall objects, for
   users' outfall loads.
4. **Withdrawal from a channel by channel**, not only through water allocation objects defined in the project. The
   two-way exchange takes water from losing reaches. Today each such reach needs a water allocation place of use set
   up in the SWAT+ project.

5. **To check in SWAT+ 62.0.1's water allocation:** `wallo_withdraw.f90` reads a place of use's duty as
   `pou(ipou)%pod(ipod)%duty`, with `ipod` the point of diversion's global number, though `pou(ipou)%pod` is the
   place of use's own list. With one point of diversion per place of use, as Basin Ledger writes them, the two agree
   for the first channel. A place of use on a second channel may read the wrong entry. Basin Ledger's tests use one
   channel. A test on a two-channel project would settle it.

A note for the swatplus-bmi docs: SWAT+ keeps one point of diversion per channel (`sd_ch%wallo_pod`), so several
places of use on one channel must share one `point_of_diver.wro` row listing them all. Separate rows for the same
channel leave all of them withdrawing nothing, with no error.

A second note for the swatplus-bmi docs: in print.prt, SWAT+ reads object rows by position unless `use_obj_labels` (the
second value under `csvout dbout cdfout`) is `y`. Osu_1hru's rows aren't in the order SWAT+ reads them, so turning on
one object's daily output by editing its row does nothing unless labels are on. swatplus-bmi's tests turn everything
on, which hides this. (An earlier version of this file reported daily outputs missing after an early stop. That was
this issue. A run stopped early writes its daily outputs in full.) *Done in 0.1.1:
`swatplus_bmi.project.set_print` turns labels on and edits rows by name, and Basin Ledger's tests use it.*
