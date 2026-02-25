import pytest
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
