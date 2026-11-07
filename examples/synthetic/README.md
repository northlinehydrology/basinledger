# Synthetic example

The engine's test case. It wires real models, but the place and the numbers mean nothing.

- **SWAT+:** the Osu_1hru project that ships with SWAT+ (one 10 ha HRU, one aquifer, 2010–2023), run one day at a time
  through swatplus-bmi. Its aquifer's daily recharge rate is the only thing passed on.
- **MODFLOW 6** (`model.py`, built with flopy): one layer, 41 × 41 cells of 50 m, 100 m thick, K = 5 m/d, specific
  yield 0.15, transmissivity held constant. A river runs along the west edge (stage 88 m); the other edges are
  no-flow. The recharge rate from SWAT+ is applied to every cell (the map `all:aqu001`), so the grid receives about
  4.2 km² × Osu's recharge. A real example maps matching areas from the same geography.
- **Water user:** the synthetic user (`basinledger.testing.SyntheticUser`, `user.toml`), hourly. It withdraws about
  520 m³ a day from well `pw1`, with daily and seasonal cycles. 70% is consumed, and the rest plus a 1% gain goes to
  an outfall account. It also reports a dissolved solids load and electricity.
- **Calibration:** none, and no record, so the engine labels it a screening run. Its objects keep their SWAT+
  names; it has no hydrofabric.

Two more configs use the same SWAT+ project:

- **`surface.toml`:** surface routes, a utility and rules.
  - A mill (`mill.toml`) takes water from Osu_1hru's channel through a SWAT+ water allocation place of use and
    returns 70% of it through a point source. It is curtailed the day after the channel's outflow falls below
    0.004 m³/s.
  - The plant is on a utility, `city`. The utility pumps well `pw1`, adds existing customers (0.002 m³/s) and 10%
    distribution losses, and discharges its sewer flow to the channel through a treatment plant's point source.
  - The plant runs dry the day after low flow, and its withdrawal is capped the day after the head at the well falls
    below 88.3 m.
- **`two_way.toml`:** the stream-aquifer exchange, both ways.
  - The MODFLOW 6 grid is 10 × 10 cells of 1,000 m², exactly Osu_1hru's 10 ha aquifer, so the recharge maps area for
    area (`area_rtol`).
  - SWAT+'s aquifer passes its recharge through. MODFLOW 6's river exchange enters SWAT+'s channel the next day, as a
    point source where the river gains and as a water allocation demand where it loses.
  - A small user pumps a well near the river (`small_user.toml`).
- **`transport.toml`:** `two_way.toml` with chloride through the aquifer (MODFLOW 6 transport, `coupling/transport.py`).
  - Porosity 0.25, longitudinal dispersivity 5 m, a start of 0.01 kg/m³, all [A].
  - Recharge carries 0.03 kg/m³ and the river's losses 0.02 kg/m³ [A], so the aquifer's chloride rises.
  - The user (`intake_user.toml`) takes in the aquifer's chloride at its well, a day late, and returns all of it.

Run it:

```bash
uv run python examples/synthetic/fetch.py --get                     # SWAT+ (2 compile jobs) and libmf6, if missing
uv run pytest                                                       # the tests use these configs
uv run python -m basinledger run examples/synthetic/surface.toml    # outputs in runs/synthetic-surface/base/outputs
uv run --group examples marimo edit examples/synthetic/notebook.py  # the whole period, with and without the user
```

What the notebook shows: drawdown at the well settles at about 1.1–1.5 m with a seasonal swing. Streamflow depletion
reaches the mean pumping rate after about five years, because every cubic meter pumped eventually comes from the
river in this model.

![Recharge and pumping](../../docs/figures/synthetic-recharge.png)
![Drawdown](../../docs/figures/synthetic-drawdown.png)
![Depletion](../../docs/figures/synthetic-depletion.png)
