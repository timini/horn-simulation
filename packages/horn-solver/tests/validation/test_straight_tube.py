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

from .conftest import assert_spl_within_tolerance

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
