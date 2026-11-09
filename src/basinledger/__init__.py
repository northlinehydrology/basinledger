"""Basin Ledger: a basin's water, quantity and quality, balanced day by day across coupled models.

SWAT+ (through swatplus-bmi), MODFLOW 6 (through libmf6 and modflowapi) and water-user components (from
waterledger-models) are stepped one day at a time by one Python driver, on NextGen geography. See docs/design.md.
"""

__version__ = "0.4.0"
