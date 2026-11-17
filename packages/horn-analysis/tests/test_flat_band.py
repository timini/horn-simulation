"""Widest flat band on synthetic responses and on a recorded FEM response."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from horn_analysis.evaluation import evaluate_response
from horn_analysis.flat_band import flat_band_fields, format_flat_bands, widest_flat_band
from horn_analysis.scoring import TargetSpec

ROOT = Path(__file__).resolve().parents[3]
F = np.geomspace(100., 10000., 401)


def _bandpass(f, low, high, order=2):
    """Level of a Butterworth-like band-pass, 0 dB in the middle."""
    return -10 * np.log10(1 + (low / f) ** (2 * order)) - 10 * np.log10(1 + (f / high) ** (2 * order))


def test_flat_response_spans_the_whole_sweep_and_is_flagged_as_a_bound():
    band = widest_flat_band(F, np.full(len(F), 90.), 1.)
    assert band["f_low_hz"] == pytest.approx(100.) and band["f_high_hz"] == pytest.approx(10000.)
    assert band["octaves"] == pytest.approx(np.log2(100.))
    assert band["low_is_bound"] and band["high_is_bound"]
    assert band["mean_level_db"] == pytest.approx(90.)


def test_bandpass_edges_match_the_analytic_crossings():
    levels = _bandpass(F, 300., 3000., order=4)
    band = widest_flat_band(F, levels, 1.5)
    # The middle is ~0 dB, so a 3 dB corridor ends where each skirt is 3 dB down.
    low = 300. * (10 ** 0.3 - 1) ** (1 / 8)
    high = 3000. / (10 ** 0.3 - 1) ** (1 / 8)
    assert band["f_low_hz"] == pytest.approx(low, rel=0.01)
    assert band["f_high_hz"] == pytest.approx(high, rel=0.01)
    assert not band["low_is_bound"] and not band["high_is_bound"]


def test_wider_tolerance_never_gives_a_narrower_band():
    rng = np.random.default_rng(3)
    levels = _bandpass(F, 250., 4000.) + rng.normal(0, 0.4, len(F))
    widths = [widest_flat_band(F, levels, t)["octaves"] for t in (0.5, 1., 2., 3., 6.)]
    assert widths == sorted(widths)


def test_a_narrow_notch_is_not_smoothed_away():
    levels = np.full(len(F), 90.)
    notch = np.argmin(abs(F - 800.))
    levels[notch] = 80.
    band = widest_flat_band(F, levels, 3.)
    # Best band is the longer side of the notch, above 800 Hz.
    assert 800. < band["f_low_hz"] < 820. and band["f_high_hz"] == pytest.approx(10000.)
    assert band["octaves"] == pytest.approx(np.log2(10000. / 800.), abs=0.03)


def test_minimum_level_breaks_the_band():
    levels = np.where(F < 1000., 80., 90.)
    assert widest_flat_band(F, levels, 10.)["f_low_hz"] == pytest.approx(100.)
    assert widest_flat_band(F, levels, 10., min_level_db=85.)["f_low_hz"] > 990.
    assert widest_flat_band(F, levels, 10., min_level_db=95.)["f_low_hz"] is None


def test_ties_go_to_the_louder_band():
    levels = np.where(F < 1000., 80., 90.)  # two equal one-decade halves
    band = widest_flat_band(F, levels, 1.)
    assert band["mean_level_db"] == pytest.approx(90., abs=0.1)


def test_rejects_unordered_input():
    with pytest.raises(ValueError):
        widest_flat_band(F[::-1], np.zeros(len(F)), 1.)
    with pytest.raises(ValueError):
        widest_flat_band(F, np.zeros(len(F)), 0.)


def test_recorded_fem_response_contains_the_target_band():
    """examples/6nmb420-320-5000 has 4.42 dB peak to peak over 320-5000 Hz."""
    data = pd.read_csv(ROOT / "examples/6nmb420-320-5000/response.csv")
    f, levels = data.frequency_hz.to_numpy(), data.spl_db_2p83vrms_1m.to_numpy()
    fields = flat_band_fields(f, levels)
    assert fields["flat_band_3db_low_hz"] <= 320. and fields["flat_band_3db_high_hz"] >= 5000.
    assert fields["flat_band_3db_octaves"] >= np.log2(5000. / 320.)
    assert fields["flat_band_1db_octaves"] < fields["flat_band_2db_octaves"] <= fields["flat_band_3db_octaves"]
    # The band agrees with the peak-to-peak of the raw samples inside it.
    inside = (f >= fields["flat_band_3db_low_hz"]) & (f <= fields["flat_band_3db_high_hz"])
    assert np.ptp(levels[inside]) <= 6. + 1e-9
    assert "+-3 dB:" in format_flat_bands(fields)


def test_ranking_rows_carry_the_flat_band_fields():
    levels = _bandpass(F, 300., 3000., order=4) + 100.
    row = evaluate_response(F, levels, TargetSpec(500, 2000))
    for tolerance in (1, 2, 3):
        assert row[f"flat_band_{tolerance}db_octaves"] > 0
    assert row["flat_band_3db_low_hz"] < 500. and row["flat_band_3db_high_hz"] > 2000.
    assert row["flat_band_min_level_db"] is None
    capped = evaluate_response(F, levels, TargetSpec(500, 2000, min_output_db=99.9))
    assert capped["flat_band_min_level_db"] == 99.9
    assert capped["flat_band_3db_octaves"] < row["flat_band_3db_octaves"]


def test_webster_screen_rows_skip_the_flat_band_fields():
    levels = _bandpass(F, 300., 3000., order=4) + 100.
    row = evaluate_response(F, levels, TargetSpec(500, 2000), flat_bands=False)
    assert not any(key.startswith("flat_band_") for key in row)
    assert format_flat_bands(row) == ""


def test_too_sparse_widest_window_falls_back_to_a_narrower_valid_band():
    # The widest +-3 dB window (100-1000 Hz) holds only two solver samples;
    # the narrower 1000-1030 Hz band holds four and must be returned.
    band = widest_flat_band([100., 1000., 1010., 1020., 1030.], [90., 95., 80., 80., 80.], 3.)
    assert band["samples"] >= 3
    assert band["f_high_hz"] == pytest.approx(1030.)
    assert band["f_low_hz"] < 1010.


def test_edge_between_sparse_samples_is_the_linear_interpolation():
    band = widest_flat_band([1000., 2000., 4000.], [0., 0., -20.], 3., min_samples=2)
    # 6 dB of a 20 dB drop over one octave: 2 ** (1 + 0.3) kHz.
    assert band["f_high_hz"] == pytest.approx(2000. * 2 ** 0.3, rel=0.01)


def test_html_cell_marks_lower_bound_and_missing_band():
    from horn_analysis.html_report import _flat_band_cell
    row = {"flat_band_3db_octaves": 1.45, "flat_band_3db_low_hz": 566., "flat_band_3db_high_hz": 1543.,
           "flat_band_3db_is_lower_bound": True}
    assert _flat_band_cell(row) == "≥1.45 (566 to 1543 Hz)"
    assert _flat_band_cell({"flat_band_3db_octaves": 0., "flat_band_3db_low_hz": None}) == "none"
    assert _flat_band_cell({}) == "n/a"
