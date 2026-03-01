"""BEM vs analytical radiation model validation tests.

Compares BEM radiation coupling against analytical models and exact solutions
to quantify the accuracy improvement from nonlocal radiation conditions.

Test cases:
  V7a: Straight tube — BEM should match analytical ~94 dB constant SPL
  V7b: Conical horn at high ka — BEM and flanged_piston should converge
  V7c: Conical horn at low ka — quantify the accuracy delta
"""

import numpy as np
import pytest

try:
    import bempp.api as bempp_api
    BEMPP_AVAILABLE = True
except ImportError:
    BEMPP_AVAILABLE = False

from .conftest import (
    assert_spl_within_tolerance,
    run_solver_and_get_spl,
    load_reference,
    _generate_cylinder_step,
)

pytestmark = [
    pytest.mark.validation,
    pytest.mark.skipif(not BEMPP_AVAILABLE, reason="bempp-cl not installed"),
]

EXPECTED_SPL_STRAIGHT_TUBE = 20 * np.log10(1.0 / 20e-6)  # ~93.98 dB


@pytest.fixture(scope="module")
def straight_tube_bem_results(tmp_path_factory):
    """Solve V7a: straight tube with BEM radiation BC."""
    ref = load_reference("straight_tube_analytical.json")
    geom = ref["geometry"]
    freq_cfg = ref["frequency_range"]
    tmp = tmp_path_factory.mktemp("v7a")
    step_file = tmp / "straight_tube.step"
    _generate_cylinder_step(step_file, geom["throat_radius_m"], geom["length_m"])
    frequencies, spl = run_solver_and_get_spl(
        step_file=step_file,
        freq_range=(freq_cfg["min_hz"], freq_cfg["max_hz"]),
        num_intervals=freq_cfg["num_points"],
        horn_length=geom["length_m"],
        tmp_dir=tmp,
        radiation_model="bem",
    )
    return frequencies, spl, ref


class TestStraightTubeBEM:
    """V7a: BEM radiation on a straight tube should match analytical."""

    def test_spl_matches_analytical(self, straight_tube_bem_results):
        """BEM SPL should be within 2 dB of analytical 93.98 dB.

        We use a wider tolerance than the Robin BC test (0.5 dB) because
        BEM has additional discretization error from the boundary mesh and
        the iterative coupling scheme.
        """
        frequencies, spl, ref = straight_tube_bem_results
        tolerance = 2.0  # dB — wider than Robin BC due to BEM discretization

        assert_spl_within_tolerance(
            computed=spl,
            reference=EXPECTED_SPL_STRAIGHT_TUBE,
            tolerance_db=tolerance,
            frequencies=frequencies,
            label="V7a straight tube (BEM)",
        )

    def test_spl_is_reasonably_flat(self, straight_tube_bem_results):
        """BEM SPL should be reasonably flat (std < 2 dB)."""
        frequencies, spl, ref = straight_tube_bem_results
        std_spl = np.std(spl)
        assert std_spl < 2.0, (
            f"BEM SPL should be reasonably flat (std < 2 dB) but std={std_spl:.3f} dB.\n"
            f"SPL range: {np.min(spl):.2f} - {np.max(spl):.2f} dB"
        )


@pytest.fixture(scope="module")
def conical_horn_bem_results(tmp_path_factory):
    """Solve V7b/c: conical horn with BEM radiation BC."""
    ref = load_reference("conical_horn_webster.json")
    geom = ref["geometry"]
    freq_cfg = ref["frequency_range"]
    tmp = tmp_path_factory.mktemp("v7bc")
    step_file = tmp / "conical_horn.step"
    from horn_geometry.generator import create_conical_horn
    create_conical_horn(
        throat_radius=geom["throat_radius_m"],
        mouth_radius=geom["mouth_radius_m"],
        length=geom["length_m"],
        output_file=step_file,
    )
    frequencies, spl = run_solver_and_get_spl(
        step_file=step_file,
        freq_range=(freq_cfg["min_hz"], freq_cfg["max_hz"]),
        num_intervals=freq_cfg["num_points"],
        horn_length=geom["length_m"],
        tmp_dir=tmp,
        radiation_model="bem",
    )
    return frequencies, spl, ref


@pytest.fixture(scope="module")
def conical_horn_flanged_results(tmp_path_factory):
    """Solve conical horn with flanged_piston for comparison."""
    ref = load_reference("conical_horn_webster.json")
    geom = ref["geometry"]
    freq_cfg = ref["frequency_range"]
    tmp = tmp_path_factory.mktemp("v7_flanged")
    step_file = tmp / "conical_horn.step"
    from horn_geometry.generator import create_conical_horn
    create_conical_horn(
        throat_radius=geom["throat_radius_m"],
        mouth_radius=geom["mouth_radius_m"],
        length=geom["length_m"],
        output_file=step_file,
    )
    frequencies, spl = run_solver_and_get_spl(
        step_file=step_file,
        freq_range=(freq_cfg["min_hz"], freq_cfg["max_hz"]),
        num_intervals=freq_cfg["num_points"],
        horn_length=geom["length_m"],
        tmp_dir=tmp,
        radiation_model="flanged_piston",
    )
    return frequencies, spl, ref


class TestConicalHornBEMvsAnalytical:
    """V7b/c: Compare BEM vs analytical on a conical horn."""

    def test_bem_produces_finite_spl(self, conical_horn_bem_results):
        """BEM should produce finite SPL at all frequencies."""
        frequencies, spl, ref = conical_horn_bem_results
        assert all(np.isfinite(spl)), "All BEM SPL values should be finite"

    def test_converges_to_flanged_piston_at_high_ka(
        self, conical_horn_bem_results, conical_horn_flanged_results
    ):
        """At high ka, BEM and flanged_piston should give similar results.

        Both approximate the true radiation condition, and at high frequencies
        (large ka) both approach the plane-wave limit. The difference should
        be small (< 3 dB).
        """
        freq_bem, spl_bem, ref = conical_horn_bem_results
        freq_fp, spl_fp, _ = conical_horn_flanged_results

        geom = ref["geometry"]
        a = geom["mouth_radius_m"]
        c = 343.0

        ka_values = 2 * np.pi * freq_bem / c * a
        high_ka_mask = ka_values > 3.0

        if not np.any(high_ka_mask):
            pytest.skip("No frequencies with ka > 3 in the test range")

        diff = np.abs(spl_bem[high_ka_mask] - spl_fp[high_ka_mask])
        max_diff = np.max(diff)

        print(f"High ka (ka > 3): BEM vs flanged_piston max diff = {max_diff:.2f} dB")

        assert max_diff < 3.0, (
            f"At high ka, BEM and flanged_piston should converge. "
            f"Max SPL difference: {max_diff:.2f} dB"
        )

    def test_low_ka_deviation_quantified(
        self, conical_horn_bem_results, conical_horn_flanged_results
    ):
        """At low ka, quantify the BEM vs flanged_piston difference.

        This is the key measurement: at low frequencies where the piston
        approximation breaks down, BEM should give a different (more accurate)
        result. We just verify the comparison runs and report the delta.
        """
        freq_bem, spl_bem, ref = conical_horn_bem_results
        freq_fp, spl_fp, _ = conical_horn_flanged_results

        geom = ref["geometry"]
        a = geom["mouth_radius_m"]
        c = 343.0

        ka_values = 2 * np.pi * freq_bem / c * a
        low_ka_mask = ka_values < 1.0

        if not np.any(low_ka_mask):
            pytest.skip("No frequencies with ka < 1 in the test range")

        diff = spl_bem[low_ka_mask] - spl_fp[low_ka_mask]
        max_abs_diff = np.max(np.abs(diff))
        mean_diff = np.mean(diff)

        print(f"\nLow ka (ka < 1): BEM vs flanged_piston")
        print(f"  Mean difference: {mean_diff:+.2f} dB (positive = BEM higher)")
        print(f"  Max absolute difference: {max_abs_diff:.2f} dB")
        print(f"  This quantifies the accuracy delta from nonlocal BEM radiation.")

        # No hard assertion — this is a characterization test
        # The result tells us how much the BEM differs from the piston model
        assert np.isfinite(max_abs_diff), "Deviation should be finite"
