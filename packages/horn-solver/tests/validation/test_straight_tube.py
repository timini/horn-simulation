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

from .conftest import (
    assert_spl_within_tolerance,
    run_solver_and_get_spl,
    load_reference,
)

EXPECTED_SPL = 20 * np.log10(1.0 / 20e-6)  # ~93.98 dB


@pytest.mark.validation
class TestStraightTube:
    """V1 validation: straight tube should produce flat ~94 dB SPL."""

    def test_spl_matches_analytical(self, straight_tube_step, tmp_path):
        """All SPL values should be within 0.5 dB of the analytical 93.98 dB."""
        step_file, ref = straight_tube_step
        geom = ref["geometry"]
        freq_cfg = ref["frequency_range"]
        tolerance = ref["expected"]["tolerance_db"]

        frequencies, spl = run_solver_and_get_spl(
            step_file=step_file,
            freq_range=(freq_cfg["min_hz"], freq_cfg["max_hz"]),
            num_intervals=freq_cfg["num_points"],
            horn_length=geom["length_m"],
            mesh_size=0.005,  # 5mm for good accuracy at 4kHz (lambda/6 ~ 14mm)
            tmp_dir=tmp_path,
        )

        assert len(spl) == freq_cfg["num_points"], (
            f"Expected {freq_cfg['num_points']} frequency points, got {len(spl)}"
        )

        assert_spl_within_tolerance(
            computed=spl,
            reference=EXPECTED_SPL,
            tolerance_db=tolerance,
            frequencies=frequencies,
            label="V1 straight tube",
        )

    def test_spl_is_flat(self, straight_tube_step, tmp_path):
        """SPL standard deviation should be < 0.5 dB (nearly constant)."""
        step_file, ref = straight_tube_step
        geom = ref["geometry"]
        freq_cfg = ref["frequency_range"]

        frequencies, spl = run_solver_and_get_spl(
            step_file=step_file,
            freq_range=(freq_cfg["min_hz"], freq_cfg["max_hz"]),
            num_intervals=freq_cfg["num_points"],
            horn_length=geom["length_m"],
            mesh_size=0.005,
            tmp_dir=tmp_path,
        )

        std_spl = np.std(spl)
        assert std_spl < 0.5, (
            f"SPL should be nearly flat (std < 0.5 dB) but std={std_spl:.3f} dB.\n"
            f"SPL range: {np.min(spl):.2f} - {np.max(spl):.2f} dB"
        )
