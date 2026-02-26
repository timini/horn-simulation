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
        fig, (ax_plot, ax_table) = plt.subplots(
            2, 1, figsize=(12, 10), gridspec_kw={"height_ratios": [3, 1]},
        )
    else:
        fig, ax_plot = plt.subplots(figsize=(12, 8))

    for csv_path, label in file_label_pairs:
        df = pd.read_csv(csv_path)
        ax_plot.plot(df["frequency"], df["spl"], label=label)

    ax_plot.set_xscale("log")
    ax_plot.set_title("Horn Frequency Response Comparison")
    ax_plot.set_xlabel("Frequency (Hz)")
    ax_plot.set_ylabel("Sound Pressure Level (dB)")
    ax_plot.grid(True, which="both", ls="--")
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
                f"{kpis.f3_low_hz:.0f}" if kpis.f3_low_hz else "—",
                f"{kpis.f3_high_hz:.0f}" if kpis.f3_high_hz else "—",
                f"{kpis.bandwidth_octaves:.1f}" if kpis.bandwidth_octaves else "—",
                f"{kpis.passband_ripple_db:.1f}" if kpis.passband_ripple_db is not None else "—",
                f"{kpis.average_sensitivity_db:.1f}" if kpis.average_sensitivity_db is not None else "—",
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

    plt.tight_layout()
    output_path = Path(output_file)
    plt.savefig(output_path, dpi=150)
    plt.close()
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
