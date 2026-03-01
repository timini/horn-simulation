import json

import pytest
import numpy as np
import pandas as pd
from pathlib import Path


@pytest.fixture
def sample_csv(tmp_path):
    """Create a sample frequency response CSV."""
    csv_path = tmp_path / "results.csv"
    df = pd.DataFrame({
        "frequency": [100.0, 200.0, 500.0, 1000.0, 2000.0],
        "spl": [80.0, 85.0, 90.0, 92.0, 88.0],
    })
    df.to_csv(csv_path, index=False)
    return csv_path


class TestPlotter:
    def test_plot_creates_image(self, sample_csv, tmp_path):
        from horn_analysis.plotter import plot_spl_vs_frequency

        output = tmp_path / "plot.png"
        plot_spl_vs_frequency(str(sample_csv), str(output))
        assert output.exists()
        assert output.stat().st_size > 0

    def test_plot_with_single_point(self, tmp_path):
        from horn_analysis.plotter import plot_spl_vs_frequency

        csv_path = tmp_path / "single.csv"
        pd.DataFrame({"frequency": [1000.0], "spl": [90.0]}).to_csv(csv_path, index=False)
        output = tmp_path / "plot.png"
        plot_spl_vs_frequency(str(csv_path), str(output))
        assert output.exists()


class TestCompareHorns:
    def test_comparison_creates_image(self, sample_csv, tmp_path):
        from horn_analysis.compare_horns import plot_comparison

        csv_b = tmp_path / "results_b.csv"
        pd.DataFrame({
            "frequency": [100.0, 200.0, 500.0, 1000.0, 2000.0],
            "spl": [78.0, 83.0, 89.0, 94.0, 91.0],
        }).to_csv(csv_b, index=False)

        output = tmp_path / "comparison.png"
        plot_comparison(str(sample_csv), "Horn A", str(csv_b), "Horn B", str(output))
        assert output.exists()
        assert output.stat().st_size > 0


class TestKPI:
    """Tests for KPI extraction."""

    @pytest.fixture
    def bandpass_csv(self, tmp_path):
        """Create a CSV with a clear bandpass shape for predictable KPIs."""
        csv_path = tmp_path / "bandpass.csv"
        freq = np.geomspace(100, 10000, 100)
        # Simulate a bandpass: peak at 2kHz, rolls off below 500Hz and above 6kHz
        spl = 90 - 10 * ((np.log10(freq / 2000)) ** 2) * 20
        # Clip to a reasonable range
        spl = np.clip(spl, 50, 95)
        df = pd.DataFrame({"frequency": freq, "spl": spl})
        df.to_csv(csv_path, index=False)
        return csv_path

    def test_extract_kpis_returns_dataclass(self, bandpass_csv):
        from horn_analysis.kpi import extract_kpis, HornKPI

        kpis = extract_kpis(str(bandpass_csv))
        assert isinstance(kpis, HornKPI)

    def test_peak_detection(self, bandpass_csv):
        from horn_analysis.kpi import extract_kpis

        kpis = extract_kpis(str(bandpass_csv))
        assert kpis.peak_spl_db > 0
        assert kpis.peak_frequency_hz > 0

    def test_f3_points_bracket_peak(self, bandpass_csv):
        from horn_analysis.kpi import extract_kpis

        kpis = extract_kpis(str(bandpass_csv))
        assert kpis.f3_low_hz is not None
        assert kpis.f3_high_hz is not None
        assert kpis.f3_low_hz < kpis.peak_frequency_hz
        assert kpis.f3_high_hz > kpis.peak_frequency_hz

    def test_bandwidth_positive(self, bandpass_csv):
        from horn_analysis.kpi import extract_kpis

        kpis = extract_kpis(str(bandpass_csv))
        assert kpis.bandwidth_hz is not None
        assert kpis.bandwidth_hz > 0
        assert kpis.bandwidth_octaves is not None
        assert kpis.bandwidth_octaves > 0

    def test_passband_ripple_nonnegative(self, bandpass_csv):
        from horn_analysis.kpi import extract_kpis

        kpis = extract_kpis(str(bandpass_csv))
        assert kpis.passband_ripple_db is not None
        assert kpis.passband_ripple_db >= 0

    def test_to_dict_and_json(self, bandpass_csv):
        from horn_analysis.kpi import extract_kpis

        kpis = extract_kpis(str(bandpass_csv))
        d = kpis.to_dict()
        assert isinstance(d, dict)
        assert "peak_spl_db" in d
        # Should be JSON-serializable
        json_str = json.dumps(d)
        assert len(json_str) > 0

    def test_flat_response_no_f3(self, tmp_path):
        """A flat response at the peak should have no -3dB crossing below the peak."""
        from horn_analysis.kpi import extract_kpis

        csv_path = tmp_path / "flat.csv"
        freq = np.geomspace(100, 10000, 50)
        spl = np.full_like(freq, 90.0)
        pd.DataFrame({"frequency": freq, "spl": spl}).to_csv(csv_path, index=False)

        kpis = extract_kpis(str(csv_path))
        assert kpis.peak_spl_db == pytest.approx(90.0)
        # With perfectly flat response, there's no -3dB crossing
        assert kpis.f3_low_hz is None
        assert kpis.f3_high_hz is None


class TestMultiComparison:
    """Tests for multi-horn comparison plotting."""

    def _make_csv(self, tmp_path, name, spl_offset=0):
        csv_path = tmp_path / name
        freq = np.geomspace(100, 10000, 50)
        spl = 85 + spl_offset - 5 * ((np.log10(freq / 2000)) ** 2) * 10
        pd.DataFrame({"frequency": freq, "spl": spl}).to_csv(csv_path, index=False)
        return csv_path

    def test_multi_comparison_creates_image(self, tmp_path):
        from horn_analysis.compare import plot_multi_comparison

        csv_a = self._make_csv(tmp_path, "a.csv", spl_offset=0)
        csv_b = self._make_csv(tmp_path, "b.csv", spl_offset=3)
        csv_c = self._make_csv(tmp_path, "c.csv", spl_offset=-2)

        output = tmp_path / "multi.png"
        plot_multi_comparison(
            [(str(csv_a), "Horn A"), (str(csv_b), "Horn B"), (str(csv_c), "Horn C")],
            str(output),
        )
        assert output.exists()
        assert output.stat().st_size > 0

    def test_multi_comparison_with_kpi_table(self, tmp_path):
        from horn_analysis.compare import plot_multi_comparison

        csv_a = self._make_csv(tmp_path, "a.csv")
        csv_b = self._make_csv(tmp_path, "b.csv", spl_offset=5)

        output = tmp_path / "multi_kpi.png"
        plot_multi_comparison(
            [(str(csv_a), "Horn A"), (str(csv_b), "Horn B")],
            str(output),
            kpi_table=True,
        )
        assert output.exists()
        assert output.stat().st_size > 0

    def test_parse_file_label(self):
        from horn_analysis.compare import parse_file_label

        assert parse_file_label("data.csv:My Horn") == ("data.csv", "My Horn")
        path, label = parse_file_label("results/final.csv")
        assert path == "results/final.csv"
        assert label == "final"


@pytest.fixture
def solver_csv_with_impedance_phase(tmp_path):
    """Create a CSV with all solver output columns including impedance and phase."""
    csv_path = tmp_path / "full_results.csv"
    freq = np.geomspace(500, 8000, 50)
    spl = 85 - 5 * ((np.log10(freq / 2000)) ** 2) * 10
    phase_deg = -np.cumsum(np.ones_like(freq) * 15)  # Decreasing phase
    # Simulate impedance: real part ~400, imaginary varies
    z_real = 400 + 50 * np.sin(2 * np.pi * np.log10(freq))
    z_imag = 100 * np.cos(2 * np.pi * np.log10(freq))
    df = pd.DataFrame({
        "frequency": freq,
        "spl": spl,
        "phase_deg": phase_deg,
        "z_real": z_real,
        "z_imag": z_imag,
    })
    df.to_csv(csv_path, index=False)
    return csv_path


class TestImpedancePlot:
    def test_impedance_plot_creates_image(self, solver_csv_with_impedance_phase, tmp_path):
        from horn_analysis.impedance_plot import plot_impedance

        output = tmp_path / "impedance.png"
        plot_impedance(str(solver_csv_with_impedance_phase), str(output))
        assert output.exists()
        assert output.stat().st_size > 0


class TestPhasePlot:
    def test_phase_plot_creates_image(self, solver_csv_with_impedance_phase, tmp_path):
        from horn_analysis.phase_plot import plot_phase

        output = tmp_path / "phase.png"
        plot_phase(str(solver_csv_with_impedance_phase), str(output))
        assert output.exists()
        assert output.stat().st_size > 0

    def test_phase_plot_with_group_delay(self, solver_csv_with_impedance_phase, tmp_path):
        from horn_analysis.phase_plot import plot_phase

        output = tmp_path / "phase_gd.png"
        plot_phase(str(solver_csv_with_impedance_phase), str(output), group_delay=True)
        assert output.exists()
        assert output.stat().st_size > 0


class TestDashboard:
    def test_dashboard_creates_image(self, solver_csv_with_impedance_phase, tmp_path):
        from horn_analysis.dashboard import generate_dashboard

        output = tmp_path / "dashboard.png"
        generate_dashboard(str(solver_csv_with_impedance_phase), str(output))
        assert output.exists()
        assert output.stat().st_size > 0


# -- Horn 3D Render Tests ---------------------------------------------------

class TestHornRender:
    def test_render_conical_creates_image(self, tmp_path):
        from horn_analysis.horn_render import render_horn_3d

        output = tmp_path / "conical.png"
        render_horn_3d(0.05, 0.2, 0.5, "conical", str(output))
        assert output.exists()
        assert output.stat().st_size > 0

    def test_render_exponential_creates_image(self, tmp_path):
        from horn_analysis.horn_render import render_horn_3d

        output = tmp_path / "exponential.png"
        render_horn_3d(0.05, 0.2, 0.5, "exponential", str(output))
        assert output.exists()
        assert output.stat().st_size > 0

    def test_render_hyperbolic_creates_image(self, tmp_path):
        from horn_analysis.horn_render import render_horn_3d

        output = tmp_path / "hyperbolic.png"
        render_horn_3d(0.05, 0.2, 0.5, "hyperbolic", str(output))
        assert output.exists()
        assert output.stat().st_size > 0

    def test_render_without_profile_panel(self, tmp_path):
        from horn_analysis.horn_render import render_horn_3d

        output = tmp_path / "no_profile.png"
        render_horn_3d(0.05, 0.2, 0.5, "conical", str(output), show_profile=False)
        assert output.exists()
        assert output.stat().st_size > 0

    def test_radius_profile_values(self):
        from horn_analysis.horn_render import _radius_profile

        r_t, r_m, L = 0.05, 0.2, 0.5
        for profile in ("conical", "exponential", "hyperbolic", "os", "cd"):
            r = _radius_profile(np.array([0.0, L]), r_t, r_m, L, profile)
            assert r[0] == pytest.approx(r_t, rel=1e-10)
            assert r[-1] == pytest.approx(r_m, rel=1e-10)
        # Tractrix and Le Cléac'h use parametric interpolation; endpoints are approximate
        for profile in ("tractrix", "lecleach"):
            r = _radius_profile(np.array([0.0, L]), r_t, r_m, L, profile)
            assert r[0] == pytest.approx(r_t, rel=0.02)
            assert r[-1] == pytest.approx(r_m, rel=0.02)

    def test_radius_profile_unknown_raises(self):
        from horn_analysis.horn_render import _radius_profile

        with pytest.raises(ValueError, match="Unknown profile"):
            _radius_profile(np.array([0.0]), 0.05, 0.2, 0.5, "parabolic")

    def test_fig_to_b64_3d_returns_data_uri(self):
        from horn_analysis.horn_render import fig_to_b64_3d

        result = fig_to_b64_3d(0.05, 0.2, 0.5, "conical")
        assert result.startswith("data:image/png;base64,")
        assert len(result) > 100


# -- Directivity Plot Tests --------------------------------------------------

@pytest.fixture
def directivity_csv(tmp_path):
    """Create a synthetic directivity CSV with cos^2 pattern narrowing with frequency."""
    csv_path = tmp_path / "directivity.csv"
    freqs = [500, 1000, 2000, 4000, 8000]
    angles = np.arange(0, 181, 5)
    rows = []
    for f in freqs:
        # Narrower beam at higher frequency
        n = f / 500  # exponent increases with frequency
        for theta in angles:
            # cos^n pattern gives narrowing beam
            spl = 90 + 10 * np.log10(max(np.cos(np.radians(theta)) ** n, 1e-10))
            rows.append({"frequency": f, "theta_deg": theta, "spl_db": spl})
    df = pd.DataFrame(rows)
    df.to_csv(csv_path, index=False)
    return csv_path


class TestDirectivityPlot:
    def test_load_directivity(self, directivity_csv):
        from horn_analysis.directivity_plot import load_directivity

        df = load_directivity(str(directivity_csv))
        assert "frequency" in df.columns
        assert "theta_deg" in df.columns
        assert "spl_db" in df.columns

    def test_polar_directivity(self, directivity_csv, tmp_path):
        from horn_analysis.directivity_plot import load_directivity, plot_polar_directivity

        df = load_directivity(str(directivity_csv))
        output = tmp_path / "polar.png"
        plot_polar_directivity(df, output_file=str(output))
        assert output.exists()
        assert output.stat().st_size > 0

    def test_directivity_contour(self, directivity_csv, tmp_path):
        from horn_analysis.directivity_plot import load_directivity, plot_directivity_contour

        df = load_directivity(str(directivity_csv))
        output = tmp_path / "contour.png"
        plot_directivity_contour(df, output_file=str(output))
        assert output.exists()
        assert output.stat().st_size > 0

    def test_beamwidth(self, directivity_csv, tmp_path):
        from horn_analysis.directivity_plot import load_directivity, plot_beamwidth

        df = load_directivity(str(directivity_csv))
        output = tmp_path / "beamwidth.png"
        plot_beamwidth(df, output_file=str(output))
        assert output.exists()
        assert output.stat().st_size > 0

    def test_directivity_index(self, directivity_csv, tmp_path):
        from horn_analysis.directivity_plot import load_directivity, plot_directivity_index

        df = load_directivity(str(directivity_csv))
        output = tmp_path / "di.png"
        plot_directivity_index(df, output_file=str(output))
        assert output.exists()
        assert output.stat().st_size > 0

    def test_compute_beamwidth_values(self, directivity_csv):
        from horn_analysis.directivity_plot import load_directivity, compute_beamwidth

        df = load_directivity(str(directivity_csv))
        freqs, bw = compute_beamwidth(df)
        # Higher frequency should have narrower beamwidth
        assert len(freqs) == 5
        # Compare lowest and highest frequency beamwidths
        assert bw[-1] < bw[0], "Higher frequency should have narrower beam"

    def test_compute_directivity_index_values(self, directivity_csv):
        from horn_analysis.directivity_plot import load_directivity, compute_directivity_index

        df = load_directivity(str(directivity_csv))
        freqs, di = compute_directivity_index(df)
        # Horn concentrates sound -> DI should be positive
        assert np.all(di > 0), "Directivity index should be positive for a horn"

    def test_generate_directivity_report(self, directivity_csv, tmp_path):
        from horn_analysis.directivity_plot import generate_directivity_report

        output_dir = tmp_path / "directivity_report"
        generate_directivity_report(str(directivity_csv), str(output_dir))

        expected_files = [
            "polar_directivity.png",
            "directivity_contour.png",
            "beamwidth.png",
            "directivity_index.png",
        ]
        for fname in expected_files:
            fpath = output_dir / fname
            assert fpath.exists(), f"Missing: {fname}"
            assert fpath.stat().st_size > 0, f"Empty: {fname}"
