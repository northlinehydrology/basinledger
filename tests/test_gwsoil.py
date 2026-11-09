"""Groundwater up into the soil, the parts built before the rule is chosen (coupling/gwsoil.py): a cell's water split
among the HRUs over it by their shares of its ground, and among each HRU's layers below the water table by their
thickness below it; the profile's bottom under each cell. Volume is kept exactly."""

from __future__ import annotations

import numpy as np

from basinledger.coupling import gwsoil as G
from basinledger.mapping import AreaMap


def two_hrus():
    # cell 0: HRU a 60 m2, HRU b 40 m2; cell 1: HRU b 100 m2; cell 2: nothing
    return AreaMap(["a", "b"], np.array([0, 1, 1]), np.array([0, 0, 1]), np.array([60.0, 40.0, 100.0]))


def test_a_cells_water_goes_to_its_hrus_by_their_share_of_its_ground():
    s = G.cell_shares(two_hrus(), 3)
    np.testing.assert_allclose(s, [0.6, 0.4, 1.0])


def test_the_layers_below_the_water_table_take_it_by_thickness():
    tops, bottoms = np.array([0.0, 0.3, 1.0]), np.array([0.3, 1.0, 2.0])
    np.testing.assert_allclose(G.layer_split(tops, bottoms, 0.5), np.array([0.0, 0.5, 1.0]) / 1.5)
    np.testing.assert_allclose(G.layer_split(tops, bottoms, 0.0), [0.15, 0.35, 0.5])   # at the surface: all of it
    np.testing.assert_array_equal(G.layer_split(tops, bottoms, 2.5), [0.0, 0.0, 0.0])  # below the profile


def test_volume_is_kept_from_cells_to_layers():
    m = two_hrus()
    layer_hru = np.array([0, 0, 1, 1, 1])
    layer_bottom = np.array([0.4, 1.5, 0.2, 0.8, 1.2])
    cell_volume = np.array([10.0, 5.0, 0.0])
    row_volume = cell_volume[m.node] * G.cell_shares(m, 3)
    head_depth = np.array([0.5, 1.0, 9.0])           # depth to water per cell
    vol, unplaced = G.to_layers(row_volume, head_depth[m.node], m, layer_hru, layer_bottom)
    assert unplaced == 0.0
    assert vol.sum() == np.float64(row_volume.sum())
    # HRU a, 6 m3 at 0.5 m: its first layer is above the water table
    np.testing.assert_allclose(vol[:2], [0.0, 6.0])
    # HRU b: 4 m3 at 0.5 m (layers 0.3 and 0.4 m below it), 5 m3 at 1.0 m (only the 0.2 m of its last layer)
    np.testing.assert_allclose(vol[2:], [0.0, 4.0 * 0.3 / 0.7, 4.0 * 0.4 / 0.7 + 5.0])
    rates = G.layer_rates(vol, layer_hru, np.array([60.0, 140.0]))
    np.testing.assert_allclose((rates * 86400.0 * np.array([60.0, 140.0])[layer_hru]).sum(), 15.0)


def test_water_with_its_table_below_the_profile_is_returned_unplaced():
    m = two_hrus()
    vol, unplaced = G.to_layers(np.array([1.0, 2.0, 3.0]), np.array([5.0, 5.0, 0.1]), m, np.array([0, 1]),
                                np.array([1.0, 1.0]))
    assert unplaced == 3.0 and vol.tolist() == [0.0, 3.0]


def test_the_drain_sits_at_the_profiles_bottom_under_each_cell():
    m = two_hrus()
    depth = G.profile_depth(np.array([0, 0, 1, 1, 1]), np.array([0.4, 1.5, 0.2, 0.8, 1.2]), 2)
    np.testing.assert_allclose(depth, [1.5, 1.2])
    e = G.drain_elevation(m, np.array([100.0, 101.0, 102.0]), depth, 3)
    np.testing.assert_allclose(e[:2], [100.0 - (0.6 * 1.5 + 0.4 * 1.2), 101.0 - 1.2])
    assert np.isnan(e[2])
