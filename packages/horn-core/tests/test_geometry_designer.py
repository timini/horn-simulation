"""Tests for analytical geometry derivation (fullauto mode)."""

import json
import math

import pytest

from horn_core.geometry_designer import (
    C0,
    DerivedGeometry,
    derive_length_range,
    derive_mouth_radius,
    derive_mouth_radius_range,
    derive_simulation_freq_range,
    generate_fullauto_candidates,
)


class TestDeriveMouthRadius:
    def test_formula_500hz(self):
        """At 500 Hz: r = 343 / (2*pi*500) = 0.1092 m."""
        expected = C0 / (2 * math.pi * 500)
        assert derive_mouth_radius(500) == pytest.approx(expected, rel=1e-6)
        assert derive_mouth_radius(500) == pytest.approx(0.1092, rel=1e-2)

    def test_formula_1000hz(self):
        """At 1000 Hz: r = 343 / (2*pi*1000) = 0.0546 m."""
        expected = C0 / (2 * math.pi * 1000)
        assert derive_mouth_radius(1000) == pytest.approx(expected, rel=1e-6)

    def test_lower_freq_gives_larger_radius(self):
        """Lower frequency -> larger mouth radius."""
        assert derive_mouth_radius(200) > derive_mouth_radius(500)
        assert derive_mouth_radius(500) > derive_mouth_radius(2000)

    def test_formula_100hz_large_horn(self):
        """At 100 Hz: r = 343 / (2*pi*100) = 0.546 m (large horn)."""
        r = derive_mouth_radius(100)
        assert r == pytest.approx(0.546, rel=1e-2)

    def test_formula_4000hz_small_horn(self):
        """At 4000 Hz: r = 343 / (2*pi*4000) = 0.01365 m (small horn)."""
        r = derive_mouth_radius(4000)
        assert r == pytest.approx(0.01365, rel=1e-2)


class TestDeriveMouthRadiusRange:
    def test_default_spread(self):
        """Default +-30% around ideal."""
        ideal = derive_mouth_radius(500)
        lo, hi = derive_mouth_radius_range(500)
        assert lo == pytest.approx(ideal * 0.7, rel=1e-6)
        assert hi == pytest.approx(ideal * 1.3, rel=1e-6)

    def test_custom_spread(self):
        ideal = derive_mouth_radius(500)
        lo, hi = derive_mouth_radius_range(500, spread=0.5)
        assert lo == pytest.approx(ideal * 0.5, rel=1e-6)
        assert hi == pytest.approx(ideal * 1.5, rel=1e-6)

    def test_range_brackets_ideal(self):
        ideal = derive_mouth_radius(1000)
        lo, hi = derive_mouth_radius_range(1000)
        assert lo < ideal < hi


class TestDeriveLengthRange:
    def test_formula_500hz(self):
        """At 500 Hz: wavelength=0.686m, L_min=0.1715, L_max=0.343."""
        wavelength = C0 / 500
        l_min, l_max = derive_length_range(500)
        assert l_min == pytest.approx(wavelength / 4, rel=1e-6)
        assert l_max == pytest.approx(wavelength / 2, rel=1e-6)

    def test_lower_freq_gives_longer_horn(self):
        lo_min, lo_max = derive_length_range(200)
        hi_min, hi_max = derive_length_range(2000)
        assert lo_min > hi_min
        assert lo_max > hi_max

    def test_max_is_double_min(self):
        """L_max should always be 2x L_min (half-wave vs quarter-wave)."""
        l_min, l_max = derive_length_range(500)
        assert l_max == pytest.approx(2 * l_min, rel=1e-9)


class TestDeriveSimulationFreqRange:
    def test_extends_half_octave(self):
        """Should extend by +/-0.5 octave."""
        sim_min, sim_max = derive_simulation_freq_range(500, 4000)
        assert sim_min == pytest.approx(500 / math.sqrt(2), rel=1e-6)
        assert sim_max == pytest.approx(4000 * math.sqrt(2), rel=1e-6)

    def test_sim_range_brackets_target(self):
        sim_min, sim_max = derive_simulation_freq_range(500, 4000)
        assert sim_min < 500
        assert sim_max > 4000


class TestGenerateFullautoCandidates:
    def test_default_grid_size(self):
        """3 profiles x 1 throat x 3 mouth x 3 lengths = 27 max."""
        candidates, derived = generate_fullauto_candidates(
            target_f_low=500,
            target_f_high=4000,
            throat_radii=[0.025],
        )
        # All 27 should pass since derived mouth radii >> 0.025
        assert len(candidates) == 27
        assert derived.candidate_count == 27

    def test_all_profiles_represented(self):
        candidates, _ = generate_fullauto_candidates(
            target_f_low=500,
            target_f_high=4000,
            throat_radii=[0.025],
        )
        profiles = {c.profile for c in candidates}
        assert profiles == {"conical", "exponential", "hyperbolic"}

    def test_mouth_exceeds_throat(self):
        """Every candidate must have mouth_radius > throat_radius."""
        candidates, _ = generate_fullauto_candidates(
            target_f_low=500,
            target_f_high=4000,
            throat_radii=[0.025],
        )
        for c in candidates:
            assert c.mouth_radius > c.throat_radius, (
                f"{c.candidate_id}: mouth={c.mouth_radius} <= throat={c.throat_radius}"
            )

    def test_unique_candidate_ids(self):
        candidates, _ = generate_fullauto_candidates(
            target_f_low=500,
            target_f_high=4000,
            throat_radii=[0.025],
        )
        ids = [c.candidate_id for c in candidates]
        assert len(ids) == len(set(ids))

    def test_mouth_radius_within_derived_range(self):
        """All mouth radii should fall within the derived range."""
        candidates, derived = generate_fullauto_candidates(
            target_f_low=500,
            target_f_high=4000,
            throat_radii=[0.025],
        )
        lo, hi = derived.mouth_radius_range
        for c in candidates:
            # Allow 1e-6 tolerance for rounding in candidate generation
            assert lo - 1e-6 <= c.mouth_radius <= hi + 1e-6

    def test_length_within_derived_range(self):
        """All lengths should fall within the derived range."""
        candidates, derived = generate_fullauto_candidates(
            target_f_low=500,
            target_f_high=4000,
            throat_radii=[0.025],
        )
        lo, hi = derived.length_range
        for c in candidates:
            assert lo - 1e-9 <= c.length <= hi + 1e-9

    def test_custom_grid_size(self):
        """num_mouth_radii=2, num_lengths=2 -> 3*1*2*2=12 max."""
        candidates, _ = generate_fullauto_candidates(
            target_f_low=500,
            target_f_high=4000,
            throat_radii=[0.025],
            num_mouth_radii=2,
            num_lengths=2,
        )
        assert len(candidates) == 12

    def test_multiple_throat_radii(self):
        """Two throat radii should roughly double the candidates."""
        single, _ = generate_fullauto_candidates(
            target_f_low=500,
            target_f_high=4000,
            throat_radii=[0.025],
        )
        double, _ = generate_fullauto_candidates(
            target_f_low=500,
            target_f_high=4000,
            throat_radii=[0.020, 0.025],
        )
        assert len(double) == 2 * len(single)

    def test_large_throat_filters_small_mouths(self):
        """When throat is large relative to derived mouth radius,
        some candidates should be filtered."""
        # At 4000 Hz, ideal mouth radius ~0.0137m. With throat=0.012,
        # some mouth radii in the range will be <= throat.
        candidates, derived = generate_fullauto_candidates(
            target_f_low=4000,
            target_f_high=8000,
            throat_radii=[0.012],
        )
        # Should still produce some valid candidates
        assert len(candidates) > 0
        # But fewer than the full grid
        assert len(candidates) < 27

    def test_derived_geometry_populated(self):
        _, derived = generate_fullauto_candidates(
            target_f_low=500,
            target_f_high=4000,
            throat_radii=[0.025],
        )
        assert isinstance(derived, DerivedGeometry)
        assert derived.target_f_low == 500
        assert derived.target_f_high == 4000
        assert derived.ideal_mouth_radius > 0
        assert derived.mouth_radius_range[0] < derived.mouth_radius_range[1]
        assert derived.length_range[0] < derived.length_range[1]
        assert derived.sim_freq_range[0] < 500
        assert derived.sim_freq_range[1] > 4000

    def test_low_freq_produces_large_horn(self):
        """100 Hz target should produce mouth radius ~0.55m."""
        _, derived = generate_fullauto_candidates(
            target_f_low=100,
            target_f_high=1000,
            throat_radii=[0.025],
        )
        assert derived.ideal_mouth_radius == pytest.approx(0.546, rel=1e-2)
        assert derived.length_range[1] > 1.0  # half-wave at 100 Hz > 1.7m

    def test_high_freq_produces_small_horn(self):
        """4000 Hz target should produce mouth radius ~0.014m."""
        _, derived = generate_fullauto_candidates(
            target_f_low=4000,
            target_f_high=8000,
            throat_radii=[0.005],
        )
        assert derived.ideal_mouth_radius == pytest.approx(0.01365, rel=1e-2)
        assert derived.length_range[1] < 0.05  # half-wave at 4000 Hz ~0.043m
