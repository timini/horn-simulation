"""Widest flat band: the widest contiguous range where a response stays within +-N dB.

A response is "within +-N dB" over a range when its peak-to-peak level there
is at most 2N dB, so it fits inside a 2N dB corridor centred on the midpoint
of its highest and lowest level. This is the same peak-to-peak measure as the
target-band ripple limit, so a 6 dB ``max_ripple_db`` and a +-3 dB flat band
agree.

The search runs over the whole sweep it is given. The response is linearly
interpolated in log frequency onto a grid of every original sample plus a
fixed step (1/96 octave by default). Every original sample, including a
narrow notch, still counts. Nothing is smoothed. A band edge that falls
between two solver samples comes from the linear interpolation, so its
accuracy is set by the solver point density; the 1/96 octave step is only the
grid quantisation. A two-pointer sweep with running maximum and minimum finds the widest
window in linear time. When ``min_level_db`` is given, points below it break
the band. A band that touches either end of the sweep may extend further, so
it is flagged as a lower bound on the true width.
"""
from collections import deque

import numpy as np

DEFAULT_TOLERANCES_DB = (1.0, 2.0, 3.0)
DEFAULT_POINTS_PER_OCTAVE = 96
MIN_SAMPLES = 3


def _empty(tolerance_db):
    return {"tolerance_db": float(tolerance_db), "f_low_hz": None, "f_high_hz": None, "octaves": 0.0,
            "centre_level_db": None, "mean_level_db": None, "samples": 0,
            "low_is_bound": False, "high_is_bound": False}


def widest_flat_band(frequencies, levels, tolerance_db, min_level_db=None,
                     points_per_octave=DEFAULT_POINTS_PER_OCTAVE, min_samples=MIN_SAMPLES):
    """Widest contiguous band where ``levels`` stay within +-``tolerance_db``.

    Returns a dict with ``f_low_hz``, ``f_high_hz``, ``octaves``, the corridor
    ``centre_level_db``, the log-frequency ``mean_level_db``, the number of
    original ``samples`` inside, and ``low_is_bound``/``high_is_bound`` when
    the band reaches the end of the sweep. Only windows holding at least
    ``min_samples`` original samples are candidates; when none does, the
    frequencies are None and ``octaves`` is 0. Equal widths are resolved towards the higher mean level.
    """
    f = np.asarray(frequencies, dtype=float)
    y0 = np.asarray(levels, dtype=float)
    if f.ndim != 1 or f.shape != y0.shape or len(f) < 2:
        raise ValueError("Frequencies and levels must be matching one-dimensional arrays")
    if not (np.isfinite(f).all() and np.isfinite(y0).all()) or f[0] <= 0 or np.any(np.diff(f) <= 0):
        raise ValueError("Frequencies must be positive and strictly increasing; levels finite")
    if not np.isfinite(tolerance_db) or tolerance_db <= 0 or points_per_octave <= 0:
        raise ValueError("Tolerance and grid step must be positive")
    x0 = np.log2(f)
    x = np.union1d(x0, np.arange(x0[0], x0[-1], 1.0 / points_per_octave))
    y = np.interp(x, x0, y0)
    valid = np.ones(len(x), bool) if min_level_db is None else y >= min_level_db - 1e-12
    span = 2.0 * tolerance_db + 1e-9
    # Cumulative count of original samples, to skip windows holding too few.
    count = np.r_[0, np.cumsum(np.isin(x, x0))].tolist()
    # Cumulative trapezoid area for the window mean.
    area = np.r_[0.0, np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))]
    xs, ys = x.tolist(), y.tolist()
    best = None  # (width, mean, left, right)
    hi, lo = deque(), deque()
    left = 0
    for right in range(len(xs)):
        if not valid[right]:
            hi.clear(); lo.clear()
            left = right + 1
            continue
        while hi and ys[hi[-1]] <= ys[right]:
            hi.pop()
        hi.append(right)
        while lo and ys[lo[-1]] >= ys[right]:
            lo.pop()
        lo.append(right)
        while ys[hi[0]] - ys[lo[0]] > span:
            left += 1
            if hi[0] < left:
                hi.popleft()
            if lo[0] < left:
                lo.popleft()
        width = xs[right] - xs[left]
        if width <= 0 or count[right + 1] - count[left] < min_samples:
            continue
        mean = (area[right] - area[left]) / width
        if best is None or width > best[0] + 1e-12 or (abs(width - best[0]) <= 1e-12 and mean > best[1]):
            best = (width, mean, left, right)
    if best is None:
        return _empty(tolerance_db)
    width, mean, left, right = best
    inside = int(count[right + 1] - count[left])
    window = y[left:right + 1]
    return {"tolerance_db": float(tolerance_db), "f_low_hz": float(2.0 ** x[left]), "f_high_hz": float(2.0 ** x[right]),
            "octaves": float(width), "centre_level_db": float((window.max() + window.min()) / 2),
            "mean_level_db": float(mean), "samples": inside,
            "low_is_bound": bool(left == 0), "high_is_bound": bool(right == len(x) - 1)}


def flat_band_fields(frequencies, levels, min_level_db=None, tolerances_db=DEFAULT_TOLERANCES_DB):
    """Flat ranking-row fields ``flat_band_<N>db_*`` for each tolerance."""
    fields = {}
    for tolerance in tolerances_db:
        band = widest_flat_band(frequencies, levels, tolerance, min_level_db)
        key = f"flat_band_{tolerance:g}db"
        fields.update({
            f"{key}_low_hz": band["f_low_hz"], f"{key}_high_hz": band["f_high_hz"],
            f"{key}_octaves": band["octaves"], f"{key}_mean_level_db": band["mean_level_db"],
            f"{key}_is_lower_bound": bool(band["low_is_bound"] or band["high_is_bound"]),
        })
    fields["flat_band_min_level_db"] = min_level_db
    return fields


def format_flat_bands(row, tolerances_db=DEFAULT_TOLERANCES_DB):
    """One plain-text line, for example '+-3 dB: 0.95 oct, 410-790 Hz'; '' when absent."""
    parts = []
    for tolerance in tolerances_db:
        key = f"flat_band_{tolerance:g}db"
        if f"{key}_octaves" not in row:
            return ""
        low, high = row.get(f"{key}_low_hz"), row.get(f"{key}_high_hz")
        if low is None:
            parts.append(f"+-{tolerance:g} dB: none")
            continue
        mark = ">=" if row.get(f"{key}_is_lower_bound") else ""
        parts.append(f"+-{tolerance:g} dB: {mark}{row[f'{key}_octaves']:.2f} oct, {low:.0f}-{high:.0f} Hz")
    return "; ".join(parts)
