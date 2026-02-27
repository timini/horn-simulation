"""V6: Exponential horn community cross-validation — published geometry, independent reference.

Cross-validates the FEM solver against Webster equation predictions for an
exponential horn geometry taken from published literature (IJERT Type B horn:
throat_dia=10mm, mouth_dia=57mm, length=165mm).

This is a Tier 2 validation: FEM (3D) vs Webster (1D analytical) on a
real-world horn geometry. The Webster reference data is pre-computed and
stored in CSV to decouple this test from the V4 Webster implementation.

Source: Choudhari et al., IJERT Vol.3 Issue 2, 2014 — "Theoretical,
Simulation and Experimental Analysis of Sound Frequency and Sound Pressure
Level of Different Air Horn Amplifier" (Type B horn, modeled as exponential).

Tolerance: 6 dB (3D vs 1D model mismatch, small throat → higher-order modes).
"""

import numpy as np
import pytest

from .conftest import (
    assert_spl_within_tolerance,
    has_reference_data,
    load_csv_reference,
    run_solver_and_get_spl,
)

try:
    from horn_geometry.generator import create_exponential_horn

    GEOMETRY_AVAILABLE = True
except ImportError:
    GEOMETRY_AVAILABLE = False


@pytest.mark.validation
class TestExponentialHornCommunity:
    """V6 cross-validation: FEM vs published exponential horn geometry."""

    @pytest.mark.skipif(
        not has_reference_data("exponential_horn_ijert.csv"),
        reason="No reference data in exponential_horn_ijert.csv",
    )
    def test_exponential_horn_vs_published(self, tmp_path):
        """Compare FEM output against Webster reference for IJERT Type B horn."""
        ref_data = load_csv_reference("exponential_horn_ijert.csv")
        assert ref_data is not None, "Reference data should be loaded"
        ref_freq, ref_spl = ref_data

        # Geometry from IJERT Type B horn
        throat_radius = 0.005
        mouth_radius = 0.0285
        horn_length = 0.165

        # Generate matching exponential geometry
        step_file = tmp_path / "community_exponential_horn.step"
        create_exponential_horn(
            throat_radius=throat_radius,
            mouth_radius=mouth_radius,
            length=horn_length,
            output_file=step_file,
        )

        # Run FEM at the reference frequencies
        freq_range = (ref_freq.min(), ref_freq.max())
        frequencies, fem_spl = run_solver_and_get_spl(
            step_file=step_file,
            freq_range=freq_range,
            num_intervals=len(ref_freq),
            horn_length=horn_length,
            mesh_size=0.002,  # Smaller than throat radius (5mm)
            tmp_dir=tmp_path,
        )

        # Interpolate FEM results to reference frequency points if needed
        if not np.allclose(frequencies, ref_freq, rtol=0.01):
            fem_spl = np.interp(ref_freq, frequencies, fem_spl)
            frequencies = ref_freq

        assert_spl_within_tolerance(
            computed=fem_spl,
            reference=ref_spl,
            tolerance_db=6.0,
            frequencies=frequencies,
            label="V6 exponential community horn (FEM vs published)",
        )
