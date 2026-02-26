"""Tests for extract_kpis_from_arrays() — array-based KPI extraction."""

import pytest
import numpy as np
import pandas as pd

from horn_analysis.kpi import extract_kpis, extract_kpis_from_arrays, HornKPI


@pytest.fixture
def bandpass_data():
    """Create bandpass frequency/SPL arrays."""
    freq = np.geomspace(100, 10000, 100)
    spl = 90 - 10 * ((np.log10(freq / 2000)) ** 2) * 20
    spl = np.clip(spl, 50, 95)
    return freq, spl


@pytest.fixture
def bandpass_csv(tmp_path, bandpass_data):
    """Write bandpass data to a CSV file."""
    freq, spl = bandpass_data
    csv_path = tmp_path / "bandpass.csv"
    pd.DataFrame({"frequency": freq, "spl": spl}).to_csv(csv_path, index=False)
    return str(csv_path)


class TestExtractKpisFromArrays:
    def test_returns_horn_kpi(self, bandpass_data):
        """Should return a HornKPI dataclass."""
        freq, spl = bandpass_data
        result = extract_kpis_from_arrays(freq, spl)
        assert isinstance(result, HornKPI)

    def test_matches_csv_version(self, bandpass_data, bandpass_csv):
        """Array version should produce identical results to CSV version."""
        freq, spl = bandpass_data
        from_arrays = extract_kpis_from_arrays(freq, spl)
        from_csv = extract_kpis(bandpass_csv)

        assert from_arrays.peak_spl_db == pytest.approx(from_csv.peak_spl_db)
        assert from_arrays.peak_frequency_hz == pytest.approx(from_csv.peak_frequency_hz)

        if from_csv.f3_low_hz is not None:
            assert from_arrays.f3_low_hz == pytest.approx(from_csv.f3_low_hz, rel=1e-6)
        else:
            assert from_arrays.f3_low_hz is None

        if from_csv.f3_high_hz is not None:
            assert from_arrays.f3_high_hz == pytest.approx(from_csv.f3_high_hz, rel=1e-6)
        else:
            assert from_arrays.f3_high_hz is None

        if from_csv.bandwidth_hz is not None:
            assert from_arrays.bandwidth_hz == pytest.approx(from_csv.bandwidth_hz, rel=1e-6)

        if from_csv.bandwidth_octaves is not None:
            assert from_arrays.bandwidth_octaves == pytest.approx(from_csv.bandwidth_octaves, rel=1e-6)

        if from_csv.passband_ripple_db is not None:
            assert from_arrays.passband_ripple_db == pytest.approx(from_csv.passband_ripple_db, rel=1e-6)

        if from_csv.average_sensitivity_db is not None:
            assert from_arrays.average_sensitivity_db == pytest.approx(from_csv.average_sensitivity_db, rel=1e-6)

    def test_peak_detection(self, bandpass_data):
        """Should detect the peak correctly."""
        freq, spl = bandpass_data
        result = extract_kpis_from_arrays(freq, spl)
        assert result.peak_spl_db > 0
        assert result.peak_frequency_hz > 0

    def test_flat_response(self):
        """Flat response should have no -3dB crossings."""
        freq = np.geomspace(100, 10000, 50)
        spl = np.full_like(freq, 90.0)
        result = extract_kpis_from_arrays(freq, spl)
        assert result.peak_spl_db == pytest.approx(90.0)
        assert result.f3_low_hz is None
        assert result.f3_high_hz is None
