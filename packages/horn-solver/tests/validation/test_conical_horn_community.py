"""V3: Community horn cross-validation — published geometry, independent reference.

Cross-validates the FEM solver against Webster equation predictions for a
horn geometry taken from published literature (IJERT Type A conical horn:
throat_dia=10mm, mouth_dia=57mm, length=250mm).

This is a Tier 2 validation: FEM (3D) vs Webster (1D analytical) on a
real-world horn geometry. The Webster reference data is pre-computed and
stored in CSV to decouple this test from the Webster implementation in V2.

Source: Choudhari et al., IJERT Vol.3 Issue 2, 2014 — "Theoretical,
Simulation and Experimental Analysis of Sound Frequency and Sound Pressure
Level of Different Air Horn Amplifier"

Tolerance: 6 dB (3D vs 1D model mismatch, small throat → higher-order modes).

To add a new reference case:
    1. Digitize published frequency response data into a CSV file
    2. Add metadata as header comments (source, geometry, conditions)
    3. Place CSV in reference_data/ directory
    4. Add a new test method or parametrize the existing one
    See validation/README.md for full instructions.
"""

import csv
import numpy as np
import pytest
from pathlib import Path

from .conftest import (
    REFERENCE_DATA_DIR,
    assert_spl_within_tolerance,
    run_solver_and_get_spl,
)

try:
    from horn_geometry.generator import create_conical_horn

    GEOMETRY_AVAILABLE = True
except ImportError:
    GEOMETRY_AVAILABLE = False


def load_csv_reference(filename: str) -> tuple[np.ndarray, np.ndarray] | None:
    """Load frequency response reference data from a CSV file.

    Returns (frequencies, spl) arrays, or None if the file has no data rows.
    Header comment lines starting with '#' are skipped.
    """
    csv_path = REFERENCE_DATA_DIR / filename
    if not csv_path.exists():
        return None

    frequencies = []
    spl_values = []

    with open(csv_path) as f:
        reader = csv.DictReader(
            (row for row in f if not row.startswith("#")),
        )
        for row in reader:
            try:
                frequencies.append(float(row["frequency_hz"]))
                spl_values.append(float(row["spl_db"]))
            except (ValueError, KeyError):
                continue

    if not frequencies:
        return None

    return np.array(frequencies), np.array(spl_values)


def has_reference_data(filename: str) -> bool:
    """Check if a reference CSV file has actual data rows."""
    result = load_csv_reference(filename)
    return result is not None


@pytest.mark.validation
class TestConicalHornCommunity:
    """V3 cross-validation: FEM vs published horn measurements."""

    @pytest.mark.skipif(
        not has_reference_data("conical_horn_hornresp.csv"),
        reason="No reference data in conical_horn_hornresp.csv",
    )
    def test_conical_horn_vs_published(self, tmp_path):
        """Compare FEM output against Webster reference for IJERT Type A horn."""
        ref_data = load_csv_reference("conical_horn_hornresp.csv")
        assert ref_data is not None, "Reference data should be loaded"
        ref_freq, ref_spl = ref_data

        # Geometry parameters from CSV header comments
        # These should match the published horn being compared against.
        # Update these when populating the CSV with real data.
        throat_radius = 0.005
        mouth_radius = 0.0285
        horn_length = 0.25

        # Generate matching geometry
        step_file = tmp_path / "community_horn.step"
        create_conical_horn(
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
            mesh_size=0.004,  # Smaller than throat radius (5mm)
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
            label="V3 community horn (FEM vs published)",
        )
