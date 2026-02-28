"""V1: Straight tube with Robin BC — exact analytical validation.

Physics
-------
A straight tube (cylinder) with Dirichlet p=1 at z=0 and Robin BC
(dp/dn = -jkp) at z=L. The analytical solution is a pure forward-traveling
wave p(z) = exp(-jkz), giving |p(L)| = 1.0 at all frequencies.

Expected SPL = 20*log10(1 / 20e-6) = 93.98 dB, constant across frequency.

Tolerance: 0.5 dB (mesh discretisation error only).
"""

import numpy as np
import pytest

from .conftest import assert_spl_within_tolerance, run_solver_and_get_spl, load_reference, _generate_cylinder_step

EXPECTED_SPL = 20 * np.log10(1.0 / 20e-6)  # ~93.98 dB


@pytest.mark.validation
class TestStraightTube:
    """V1 validation: straight tube should produce flat ~94 dB SPL."""

    def test_spl_matches_analytical(self, straight_tube_results):
        """All SPL values should be within 0.5 dB of the analytical 93.98 dB."""
        frequencies, spl, ref = straight_tube_results
        tolerance = ref["expected"]["tolerance_db"]

        assert len(spl) == ref["frequency_range"]["num_points"], (
            f"Expected {ref['frequency_range']['num_points']} frequency points, got {len(spl)}"
        )

        assert_spl_within_tolerance(
            computed=spl,
            reference=EXPECTED_SPL,
            tolerance_db=tolerance,
            frequencies=frequencies,
            label="V1 straight tube",
        )

    def test_spl_is_flat(self, straight_tube_results):
        """SPL standard deviation should be < 0.5 dB (nearly constant)."""
        frequencies, spl, ref = straight_tube_results

        std_spl = np.std(spl)
        assert std_spl < 0.5, (
            f"SPL should be nearly flat (std < 0.5 dB) but std={std_spl:.3f} dB.\n"
            f"SPL range: {np.min(spl):.2f} - {np.max(spl):.2f} dB"
        )


@pytest.mark.validation
class TestStraightTubeFlangedPiston:
    """Flanged piston BC on a straight tube.

    At high ka the flanged piston impedance approaches Z=rho*c (plane wave),
    so at high frequencies the SPL should converge to the plane_wave result.
    At low ka the impedance is lower (more reflection), so SPL should differ.
    """

    @pytest.fixture(scope="class")
    def flanged_piston_results(self, tmp_path_factory):
        """Solve straight tube with flanged_piston BC."""
        ref = load_reference("straight_tube_analytical.json")
        geom = ref["geometry"]
        freq_cfg = ref["frequency_range"]
        tmp = tmp_path_factory.mktemp("v1_flanged")
        step_file = tmp / "straight_tube.step"
        _generate_cylinder_step(step_file, geom["throat_radius_m"], geom["length_m"])
        frequencies, spl = run_solver_and_get_spl(
            step_file=step_file,
            freq_range=(freq_cfg["min_hz"], freq_cfg["max_hz"]),
            num_intervals=freq_cfg["num_points"],
            horn_length=geom["length_m"],
            tmp_dir=tmp,
            radiation_model="flanged_piston",
        )
        return frequencies, spl, ref

    def test_produces_finite_spl(self, flanged_piston_results):
        """Flanged piston BC should produce finite SPL at all frequencies."""
        frequencies, spl, ref = flanged_piston_results
        assert all(np.isfinite(spl)), "All SPL values should be finite"

    def test_converges_to_plane_wave_at_high_ka(self, flanged_piston_results, straight_tube_results):
        """At high frequencies (large ka), flanged piston should converge to plane wave.

        For the straight tube geometry, we select frequencies where ka > 5
        and verify the SPL difference is small.
        """
        freq_fp, spl_fp, ref = flanged_piston_results
        freq_pw, spl_pw, _ = straight_tube_results

        geom = ref["geometry"]
        a = geom["throat_radius_m"]  # tube radius
        c = 343.0

        # Find indices where ka > 5
        ka_values = 2 * np.pi * freq_fp / c * a
        high_ka_mask = ka_values > 5.0

        if not np.any(high_ka_mask):
            pytest.skip("No frequencies with ka > 5 in the test range")

        spl_fp_high = spl_fp[high_ka_mask]
        spl_pw_high = spl_pw[high_ka_mask]

        max_diff = np.max(np.abs(spl_fp_high - spl_pw_high))
        assert max_diff < 2.0, (
            f"At high ka, flanged piston should converge to plane wave. "
            f"Max SPL difference: {max_diff:.2f} dB"
        )

    def test_lower_spl_at_low_ka(self, flanged_piston_results, straight_tube_results):
        """At low ka, flanged piston has lower radiation resistance -> different SPL.

        The impedance mismatch at low ka causes more reflection, which
        changes the SPL compared to the plane wave (Z=rho*c) case.
        """
        freq_fp, spl_fp, ref = flanged_piston_results
        freq_pw, spl_pw, _ = straight_tube_results

        geom = ref["geometry"]
        a = geom["throat_radius_m"]
        c = 343.0

        ka_values = 2 * np.pi * freq_fp / c * a
        low_ka_mask = ka_values < 0.5

        if not np.any(low_ka_mask):
            pytest.skip("No frequencies with ka < 0.5 in the test range")

        # At low ka, we just verify the results are different (not identical)
        spl_diff = np.abs(spl_fp[low_ka_mask] - spl_pw[low_ka_mask])
        mean_diff = np.mean(spl_diff)
        assert mean_diff > 0.01, (
            f"At low ka, flanged piston should differ from plane wave. "
            f"Mean SPL difference: {mean_diff:.4f} dB"
        )
