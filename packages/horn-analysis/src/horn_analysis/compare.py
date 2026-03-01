"""Multi-horn frequency response comparison with optional KPI table.

Generalizes compare_horns.py to N files. Accepts labeled file pairs
via CLI: --files a.csv:Label b.csv:Label ...
"""

import argparse
from pathlib import Path
from typing import List, Tuple, Optional

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from horn_analysis import plot_theme


def plot_multi_comparison(
    file_label_pairs: List[Tuple[str, str]],
    output_file: str,
    kpi_table: bool = False,
) -> Path:
    """Plot N frequency responses on the same axes.

    Args:
        file_label_pairs: List of (csv_path, label) tuples.
        output_file: Output image path.
        kpi_table: If True, add a KPI summary table below the plot.

    Returns:
        Path to the output image.
    """
    if kpi_table:
        fig, (ax_plot, ax_table) = plot_theme.create_figure(
            figsize=(12, 10), nrows=2, ncols=1, gridspec_kw={"height_ratios": [3, 1]},
        )
    else:
        fig, ax_plot = plot_theme.create_figure(figsize=(12, 8))

    all_freq = []
    all_spl = []
    for i, (csv_path, label) in enumerate(file_label_pairs):
        df = pd.read_csv(csv_path)
        color = plot_theme.MULTI_COLORS[i % len(plot_theme.MULTI_COLORS)]
        ax_plot.plot(df["frequency"], df["spl"], label=label, color=color, linewidth=1.4)
        all_freq.extend(df["frequency"].values)
        all_spl.extend(df["spl"].values)

    plot_theme.setup_freq_axis(ax_plot, min(all_freq), max(all_freq))
    plot_theme.setup_spl_axis(ax_plot, np.array(all_spl))
    plot_theme.setup_grid(ax_plot)

    ax_plot.set_title("Horn Frequency Response Comparison")
    ax_plot.legend()

    if kpi_table:
        from horn_analysis.kpi import extract_kpis

        rows = []
        for csv_path, label in file_label_pairs:
            kpis = extract_kpis(csv_path)
            rows.append([
                label,
                f"{kpis.peak_spl_db:.1f}",
                f"{kpis.peak_frequency_hz:.0f}",
                f"{kpis.f3_low_hz:.0f}" if kpis.f3_low_hz else "\u2014",
                f"{kpis.f3_high_hz:.0f}" if kpis.f3_high_hz else "\u2014",
                f"{kpis.bandwidth_octaves:.1f}" if kpis.bandwidth_octaves else "\u2014",
                f"{kpis.passband_ripple_db:.1f}" if kpis.passband_ripple_db is not None else "\u2014",
                f"{kpis.average_sensitivity_db:.1f}" if kpis.average_sensitivity_db is not None else "\u2014",
            ])

        col_labels = [
            "Horn", "Peak (dB)", "Peak Freq (Hz)",
            "f3 Low (Hz)", "f3 High (Hz)", "BW (oct)",
            "Ripple (dB)", "Avg Sens (dB)",
        ]

        ax_table.axis("off")
        table = ax_table.table(
            cellText=rows,
            colLabels=col_labels,
            loc="center",
            cellLoc="center",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 1.4)

    output_path = Path(output_file)
    plot_theme.save_figure(fig, str(output_path))
    print(f"Comparison plot saved to {output_path}")
    return output_path


def parse_file_label(s: str) -> Tuple[str, str]:
    """Parse 'path.csv:Label' into (path, label). If no colon, use filename as label."""
    if ":" in s:
        path, label = s.rsplit(":", 1)
        return path, label
    return s, Path(s).stem


def main():
    """CLI for multi-horn comparison."""
    parser = argparse.ArgumentParser(description="Compare N horn frequency responses.")
    parser.add_argument(
        "--files", nargs="+", required=True,
        help="CSV files with optional labels: file.csv:Label",
    )
    parser.add_argument("--output", type=str, default="comparison.png", help="Output image path.")
    parser.add_argument("--kpi-table", action="store_true", help="Add KPI summary table.")
    args = parser.parse_args()

    pairs = [parse_file_label(f) for f in args.files]
    plot_multi_comparison(pairs, args.output, kpi_table=args.kpi_table)


if __name__ == "__main__":
    main()
