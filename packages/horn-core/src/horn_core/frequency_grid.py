"""Shared frequency grid for parallel solver bands.

The sweep [min_freq, max_freq] is split into ``num_bands`` bands of equal
width on a logarithmic axis. Every band holds the same number of
logarithmically spaced points and shares its end points with its
neighbours, so the merged grid (after dropping duplicate joins) has uniform
logarithmic spacing across the whole sweep.

The number of points per band is the larger of
``ceil(num_intervals / num_bands)`` and the count needed to give at least
``min_points_per_octave`` intervals per octave. The solver processes,
band merging, refinement and the resolution study all use this module so
that the expected grid is defined in one place.
"""
import argparse
import math

import numpy as np

from horn_core.acoustics import validate_band

DEFAULT_MIN_POINTS_PER_OCTAVE = 24


def _check(min_freq, max_freq, num_bands, num_intervals, min_points_per_octave):
    validate_band(min_freq, max_freq)
    if int(num_bands) != num_bands or num_bands < 1:
        raise ValueError("num_bands must be a positive integer")
    if int(num_intervals) != num_intervals or num_intervals < 2:
        raise ValueError("num_intervals must be an integer of at least 2")
    if not np.isfinite(min_points_per_octave) or min_points_per_octave < 0:
        raise ValueError("min_points_per_octave must be finite and nonnegative")


def band_edges(min_freq, max_freq, num_bands):
    """Return ``num_bands + 1`` geometrically spaced band edges in Hz."""
    validate_band(min_freq, max_freq)
    if int(num_bands) != num_bands or num_bands < 1:
        raise ValueError("num_bands must be a positive integer")
    edges = min_freq * (max_freq / min_freq) ** (np.arange(num_bands + 1) / num_bands)
    edges[0], edges[-1] = min_freq, max_freq
    return edges


def points_per_band(min_freq, max_freq, num_bands, num_intervals,
                    min_points_per_octave=DEFAULT_MIN_POINTS_PER_OCTAVE):
    """Number of solver points in each band, including both band edges."""
    _check(min_freq, max_freq, num_bands, num_intervals, min_points_per_octave)
    points = max(2, math.ceil(num_intervals / num_bands))
    if min_points_per_octave > 0:
        octaves = math.log2(max_freq / min_freq) / num_bands
        # Tolerance stops an exact product such as 24 * 0.5 rounding up.
        intervals = math.ceil(min_points_per_octave * octaves - 1e-9)
        points = max(points, intervals + 1)
    return int(points)


def band_grids(min_freq, max_freq, num_bands, num_intervals,
               min_points_per_octave=DEFAULT_MIN_POINTS_PER_OCTAVE):
    """Return ``[(low, high, frequencies), ...]`` for every solver band."""
    points = points_per_band(min_freq, max_freq, num_bands, num_intervals, min_points_per_octave)
    edges = band_edges(min_freq, max_freq, num_bands)
    return [(float(lo), float(hi), np.geomspace(lo, hi, points)) for lo, hi in zip(edges[:-1], edges[1:])]


def merged_grid(min_freq, max_freq, num_bands, num_intervals,
                min_points_per_octave=DEFAULT_MIN_POINTS_PER_OCTAVE):
    """Unique sorted frequencies of all bands, as the merge step produces them."""
    grids = band_grids(min_freq, max_freq, num_bands, num_intervals, min_points_per_octave)
    return np.unique(np.concatenate([f for _, _, f in grids]))


def main():
    p = argparse.ArgumentParser(description="Print one solver band as 'low_hz high_hz points'.")
    p.add_argument("--min-freq", type=float, required=True)
    p.add_argument("--max-freq", type=float, required=True)
    p.add_argument("--num-bands", type=int, required=True)
    p.add_argument("--num-intervals", type=int, required=True)
    p.add_argument("--min-points-per-octave", type=float, default=DEFAULT_MIN_POINTS_PER_OCTAVE)
    p.add_argument("--band", type=int, required=True, help="Zero-based band index.")
    a = p.parse_args()
    grids = band_grids(a.min_freq, a.max_freq, a.num_bands, a.num_intervals, a.min_points_per_octave)
    if not 0 <= a.band < len(grids):
        p.error(f"band must be in 0..{len(grids) - 1}")
    low, high, frequencies = grids[a.band]
    # repr keeps full float precision so the solver reproduces the same grid.
    print(f"{low!r} {high!r} {len(frequencies)}")


if __name__ == "__main__":
    main()
