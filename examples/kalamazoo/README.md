# Kalamazoo River basin (Michigan): a published model converted to MODFLOW 6, taken as given

The second way in's Michigan test model (docs/design.md, "Extending an existing MODFLOW 6 model"). It isn't
coupled; no SWAT+ project is built over it. `tests/test_kalamazoo.py` checks what Basin Ledger would meet in it.

- **Source:** the Kalamazoo inset of USGS's Lake Michigan Basin regional model (HUC 04050003): Haserodt, M.J.,
  Kauffman, L.J., and Feinstein, D.T., 2018, U.S. Geological Survey data release, doi:10.5066/F76D5R5V; report
  Feinstein and others, 2018, SIR 2018-5038, doi:10.3133/sir20185038. MODFLOW-NWT.
- **Getting it:** ScienceBase serves `model.zip` (1.04 GB, all five inset models) only after a captcha. Download it by
  hand from the DOI's page (Attached Files) into `~/.cache/basinledger/kalamazoo` or `~/Downloads`, then:

  ```bash
  uv run get-modflow :flopy --repo modflow6 --subset mf5to6    # once: USGS's converter, from the MODFLOW 6 release
  uv run python examples/kalamazoo/fetch.py                     # MD5 checked, kala/ extracted, converted (13 s, 2.8 GB)
  uv run pytest tests/test_kalamazoo.py
  ```

  The converter reads a copy of the name file with Unix line endings, without MNWI (MNW2's output reporting, which
  mf5to6 doesn't convert) and the files only MNWI writes. The model's input is untouched.
- **What it holds, converted:** 8 layers × 730 × 1,230 cells of 500 ft (7.2 million); one steady period with no
  storage; SFR of 20,010 reaches, RIV (52,969 entries), DRN (29,897), GHB (137,100, its regional boundaries), fixed
  heads from IBOUND (28,424), recharge as arrays, and 710 multi-node wells (MNW2) as MAW. **Feet and days**: the
  converted DIS says `LENGTH_UNITS UNDEFINED`, `usgs.model.reference` says feet (UTM zone 16 in feet), and SFR's
  Manning's constant is for feet (1.486 × 86,400). The DIS carries no origin; the reference file has it.
- **What stops a coupled run, in order:** feet (Basin Ledger takes meters); no storage (daily steps need specific
  yield and storage, the study's [A]); its size (a run of 7.2 million cells a day is a long run; a study would cut an
  inset around its question); the MAW wells, which Basin Ledger's wells don't take yet.
