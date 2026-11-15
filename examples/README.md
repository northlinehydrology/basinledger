# Examples

Each folder is one place: a config (`basin.toml`) that wires real models for one basin and its water users, a script
that fetches its public data (`fetch.py`), and a short write-up (`README.md`). Real examples also keep their maps
between the NextGen hydrofabric and the MODFLOW 6 grid (`maps/`) and their calibration record (`calibration/`). An
example without calibrated models is a screening run, never a study result.

An example is also a study's starting point. Private studies are built the same way in the private northline
repository, never here.

| Example | Status | Calibration | What it is |
|---|---|---|---|
| [synthetic](synthetic/) | built | screening run | SWAT+'s Osu_1hru + a small MODFLOW 6 model + a synthetic water user; the tests run on it |
| [plaster-creek](plaster-creek/) | built | streamflow at USGS 04119055 calibrated (NextGen); nothing else | Plaster Creek, Kent County, Michigan: a calibrated NextGen run replayed, an assumed campus on utility water |
| [central-sands](central-sands/) | test model | the published model's own (USGS); not coupled yet | Pleasant Lake, Wisconsin Central Sands: a published MODFLOW 6 model taken as given, the second way in's first real test |
| [kalamazoo](kalamazoo/) | test model | the published model's own (USGS); not coupled | Kalamazoo River basin, Michigan: USGS's MODFLOW-NWT inset converted to MODFLOW 6, the second way in's Michigan test |
| [michigan-cold-streams](michigan-cold-streams/) | screen | none (a screen of public data) | 78 gaged Michigan basins scored for a first cold-stream study: trout streams, groundwater, temperature data, siting, irrigation, data centers |
