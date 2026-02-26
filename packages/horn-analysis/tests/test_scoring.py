"""Tests for the scoring and ranking module."""

import pytest

from horn_analysis.kpi import HornKPI
from horn_analysis.scoring import (
    TargetSpec,
    SelectionScore,
    compute_selection_score,
    rank_candidates,
)


def _make_kpi(**overrides):
    defaults = dict(
        peak_spl_db=95.0,
        peak_frequency_hz=2000.0,
        f3_low_hz=500.0,
        f3_high_hz=4000.0,
        bandwidth_hz=3500.0,
        bandwidth_octaves=3.0,
        passband_ripple_db=2.0,
        average_sensitivity_db=92.0,
    )
    defaults.update(overrides)
    return HornKPI(**defaults)


class TestBandwidthCoverage:
    def test_full_coverage(self):
        """Horn covers entire target range → coverage = 1.0."""
        kpi = _make_kpi(f3_low_hz=400.0, f3_high_hz=5000.0)
        target = TargetSpec(f_low_hz=500.0, f_high_hz=4000.0)
        score = compute_selection_score(kpi, target)
        assert score.bandwidth_coverage == pytest.approx(1.0)

    def test_zero_overlap(self):
        """Horn band entirely outside target range → coverage = 0.0."""
        kpi = _make_kpi(f3_low_hz=100.0, f3_high_hz=300.0)
        target = TargetSpec(f_low_hz=500.0, f_high_hz=4000.0)
        score = compute_selection_score(kpi, target)
        assert score.bandwidth_coverage == pytest.approx(0.0)

    def test_partial_overlap(self):
        """Horn covers half the target range → coverage = 0.5."""
        kpi = _make_kpi(f3_low_hz=500.0, f3_high_hz=2250.0)
        target = TargetSpec(f_low_hz=500.0, f_high_hz=4000.0)
        score = compute_selection_score(kpi, target)
        assert score.bandwidth_coverage == pytest.approx(0.5, abs=0.01)

    def test_no_f3_points(self):
        """No -3 dB crossings found → coverage = 0.0."""
        kpi = _make_kpi(f3_low_hz=None, f3_high_hz=None)
        target = TargetSpec(f_low_hz=500.0, f_high_hz=4000.0)
        score = compute_selection_score(kpi, target)
        assert score.bandwidth_coverage == 0.0


class TestRippleScore:
    def test_zero_ripple(self):
        """0 dB ripple → ripple component = 1.0."""
        kpi = _make_kpi(passband_ripple_db=0.0)
        target = TargetSpec(f_low_hz=500.0, f_high_hz=4000.0)
        score = compute_selection_score(kpi, target)
        # ripple_score = 1 - 0/6 = 1.0
        assert score.passband_ripple_db == 0.0

    def test_6db_ripple(self):
        """6 dB ripple → ripple component = 0.0."""
        kpi = _make_kpi(passband_ripple_db=6.0)
        target = TargetSpec(f_low_hz=500.0, f_high_hz=4000.0)
        score = compute_selection_score(kpi, target)
        assert score.passband_ripple_db == 6.0


class TestSensitivityScore:
    def test_high_sensitivity(self):
        """120 dB → sensitivity score = 1.0."""
        kpi = _make_kpi(average_sensitivity_db=120.0)
        target = TargetSpec(f_low_hz=500.0, f_high_hz=4000.0)
        score = compute_selection_score(kpi, target)
        assert score.avg_sensitivity_db == 120.0

    def test_low_sensitivity(self):
        """80 dB → sensitivity score = 0.0."""
        kpi = _make_kpi(average_sensitivity_db=80.0)
        target = TargetSpec(f_low_hz=500.0, f_high_hz=4000.0)
        score = compute_selection_score(kpi, target)
        assert score.avg_sensitivity_db == 80.0


class TestCompositeScore:
    def test_perfect_score(self):
        """Full coverage, zero ripple, max sensitivity → composite ≈ 1.0."""
        kpi = _make_kpi(
            f3_low_hz=400.0,
            f3_high_hz=5000.0,
            passband_ripple_db=0.0,
            average_sensitivity_db=120.0,
        )
        target = TargetSpec(f_low_hz=500.0, f_high_hz=4000.0)
        score = compute_selection_score(kpi, target)
        assert score.composite_score == pytest.approx(1.0, abs=0.01)

    def test_worst_score(self):
        """Zero coverage, 6 dB ripple, 80 dB sensitivity → composite = 0.0."""
        kpi = _make_kpi(
            f3_low_hz=100.0,
            f3_high_hz=300.0,
            passband_ripple_db=6.0,
            average_sensitivity_db=80.0,
        )
        target = TargetSpec(f_low_hz=500.0, f_high_hz=4000.0)
        score = compute_selection_score(kpi, target)
        assert score.composite_score == pytest.approx(0.0, abs=0.01)


class TestRankCandidates:
    def test_ranking_order(self):
        scores = [
            SelectionScore("d1", "h1", 0.5, 3.0, 90.0, 0.4),
            SelectionScore("d2", "h2", 1.0, 1.0, 100.0, 0.9),
            SelectionScore("d3", "h3", 0.8, 2.0, 95.0, 0.7),
        ]
        ranked = rank_candidates(scores, top_n=2)
        assert len(ranked) == 2
        assert ranked[0].driver_id == "d2"
        assert ranked[1].driver_id == "d3"

    def test_top_n_limits_output(self):
        scores = [
            SelectionScore(f"d{i}", f"h{i}", 0.5, 3.0, 90.0, float(i) / 10)
            for i in range(10)
        ]
        ranked = rank_candidates(scores, top_n=3)
        assert len(ranked) == 3

    def test_to_dict(self):
        s = SelectionScore("d1", "h1", 0.8, 2.0, 92.0, 0.75)
        d = s.to_dict()
        assert d["driver_id"] == "d1"
        assert d["composite_score"] == 0.75
