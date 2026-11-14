# Central Sands: Pleasant Lake (a published MODFLOW 6 model, taken as given)

The second way in's first real test model (docs/design.md, "Extending an existing MODFLOW 6 model"): a MODFLOW 6
model Basin Ledger didn't write. It isn't coupled yet; no SWAT+ project is built over it. The tests
(`tests/test_central_sands.py`) run what Basin Ledger can already do with it, and refuse by name what isn't built.

- **Source:** Fienen, M.N., Haserodt, M.J., and Leaf, A.T., 2021, MODFLOW models used to simulate groundwater flow in
  the Wisconsin Central Sands Study Area, 2012-2018: U.S. Geological Survey data release, doi:10.5066/P9BVFSGJ. The
  Pleasant Lake model (`PLEASANT_LAKE_MODEL.zip`, 16 MB; 212 MB unzipped), built with modflow-setup and SFRmaker.
- **What it holds:**
  - two groundwater models in one simulation: a parent (5 layers, 100 × 100 cells, CHD on its edges, 244
    irrigation wells, SFR of 551 reaches) and an LGR inset of 20 m cells around Pleasant Lake (5 layers, 100 × 120,
    SFR of 25 reaches, one lake of 1,382 connections with no outlets). They are joined by a GWF-GWF exchange with
    ghost nodes;
  - the mover, in that exchange, carrying the inset's two outlet reaches into the parent's SFR;
  - recharge from SWB in both models; Newton;
  - a steady first day, then 84 monthly periods (2012–2018) of five steps growing by 1.2;
  - meters and days, a structured grid in Wisconsin Transverse Mercator (origin 554,700 E, 389,300 N).
- **What Basin Ledger does with it now:** reads it as given; lists each model's packages; puts both SFR networks
  through the line-up rule with SFRmaker's lines as SWAT+'s channels. They line up: 2 lines in the inset, 48 in the
  parent (21 connections cut at channel boundaries), no channel problems. SWAT+ built on the same lines would route
  every channel (A). It also lists the parent's 10 confluences, and cuts the mover at run time (its factor set to 0).
- **What stops a coupled run, in order:** time steps longer than a day (the time discretization rewritten to daily
  steps, next); the mover (built as a cut, but the SFR exchange still refuses a package with the mover on); two
  groundwater models in one simulation (Basin Ledger couples one); the lake (LAK's rain, evaporation and runoff from
  SWAT+, possible through the API, not built).

Fetch it:

```bash
uv run python examples/central-sands/fetch.py
uv run pytest tests/test_central_sands.py
```
