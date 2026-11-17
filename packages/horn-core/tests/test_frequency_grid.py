import subprocess
import sys

import numpy as np
import pytest

from horn_core.frequency_grid import band_edges, band_grids, merged_grid, points_per_band


@pytest.mark.parametrize("low,high,bands,intervals", [
    (500 / np.sqrt(2), 4000 * np.sqrt(2), 8, 100),
    (500, 8000, 8, 100),
    (320 / np.sqrt(2), 5000 * np.sqrt(2), 2, 161),
    (100, 1000, 3, 7),
])
def test_merged_grid_has_uniform_log_spacing_across_band_joins(low, high, bands, intervals):
    grid = merged_grid(low, high, bands, intervals)
    steps = np.diff(np.log(grid))
    assert grid[0] == low and grid[-1] == high
    assert np.max(np.abs(steps - steps.mean())) < 1e-9
    assert len(grid) == bands * (points_per_band(low, high, bands, intervals) - 1) + 1


def test_bands_have_equal_log_width_and_share_edges():
    edges = band_edges(354., 5657., 8)
    assert np.allclose(np.diff(np.log(edges)), np.log(5657. / 354.) / 8, rtol=0, atol=1e-12)
    grids = band_grids(354., 5657., 8, 100)
    for (_, high, f), (low, _, g) in zip(grids[:-1], grids[1:]):
        assert high == low == f[-1] == g[0]


def test_low_octave_is_not_undersampled_relative_to_top_octave():
    low, high = 500 / np.sqrt(2), 4000 * np.sqrt(2)
    grid = merged_grid(low, high, 8, 100)
    first = np.sum(grid < 2 * low)
    last = np.sum(grid > high / 2)
    assert abs(first - last) <= 1
    assert first >= 24


def test_minimum_points_per_octave_overrides_a_coarse_request():
    # 3 points per band over 4 octaves would be 0.5 points per octave.
    assert points_per_band(100., 1600., 2, 6, 0) == 3
    assert points_per_band(100., 1600., 2, 6, 24) == 49
    grid = merged_grid(100., 1600., 2, 6, 24)
    assert np.log2(grid[1] / grid[0]) == pytest.approx(1 / 24)
    # An exact product does not round up to an extra point.
    assert points_per_band(354., 354. * 16, 8, 2, 24) == 13


@pytest.mark.parametrize("args", [(0, 10, 2, 10), (10, 5, 2, 10), (10, 100, 0, 10), (10, 100, 2, 1)])
def test_invalid_grid_requests_fail(args):
    with pytest.raises(ValueError):
        points_per_band(*args)


def test_cli_band_reproduces_library_grid_exactly():
    out = subprocess.run([sys.executable, "-m", "horn_core.frequency_grid", "--min-freq", "353.5533905932738",
                          "--max-freq", "5656.854249492381", "--num-bands", "8", "--num-intervals", "100",
                          "--band", "3"], capture_output=True, text=True, check=True).stdout.split()
    low, high, f = band_grids(353.5533905932738, 5656.854249492381, 8, 100)[3]
    assert (float(out[0]), float(out[1]), int(out[2])) == (low, high, len(f))
    assert np.array_equal(np.geomspace(float(out[0]), float(out[1]), int(out[2])), f)
