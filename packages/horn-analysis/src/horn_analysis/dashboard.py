"""Combined 3-panel dashboard: SPL + impedance + phase/group-delay on shared x-axis."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from horn_analysis import plot_theme


def generate_dashboard(csv_file: str, output_file: str):
    """Generate a combined 3-panel dashboard from solver CSV output.

    Row 1 (tall): SPL with 5 dB divisions
    Row 2: |Z| magnitude + angle(Z) on dual y-axis
    Row 3: Unwrapped phase + group delay on dual y-axis

    Expects CSV columns: frequency, spl, z_real, z_imag, phase_deg.

    Args:
        csv_file: Input CSV path with all solver output columns.
        output_file: Output image path.
    """
    df = pd.read_csv(csv_file)
    freq = df["frequency"].values

    fig, (ax_spl, ax_z, ax_ph) = plot_theme.create_figure(
        figsize=(11, 14), nrows=3, ncols=1,
        sharex=True, gridspec_kw={"height_ratios": [2, 1.2, 1.2]},
    )

    # -- Row 1: SPL --
    spl = df["spl"].values
    ax_spl.plot(freq, spl, color=plot_theme.COLORS["primary"], linewidth=1.4)
    ax_spl.set_title("Sound Pressure Level (SPL) vs. Frequency")
    plot_theme.setup_spl_axis(ax_spl, spl)
    plot_theme.setup_grid(ax_spl)

    # -- Row 2: Impedance (dual y-axis) --
    z_complex = df["z_real"].values + 1j * df["z_imag"].values
    z_mag = np.abs(z_complex)
    z_angle = np.degrees(np.angle(z_complex))

    color_mag = plot_theme.COLORS["primary"]
    ax_z.plot(freq, z_mag, color=color_mag, linewidth=1.4, label="|Z|")
    ax_z.set_ylabel("|Z| (Pa\u00b7s/m)", color=color_mag)
    ax_z.tick_params(axis="y", labelcolor=color_mag)

    ax_z2 = ax_z.twinx()
    color_angle = plot_theme.COLORS["secondary"]
    ax_z2.plot(freq, z_angle, color=color_angle, linestyle="--", linewidth=1.4, label="\u2220Z")
    ax_z2.set_ylabel("\u2220Z (degrees)", color=color_angle)
    ax_z2.tick_params(axis="y", labelcolor=color_angle)

    ax_z.set_title("Throat Impedance")
    plot_theme.setup_grid(ax_z)

    lines_z1, labels_z1 = ax_z.get_legend_handles_labels()
    lines_z2, labels_z2 = ax_z2.get_legend_handles_labels()
    ax_z.legend(lines_z1 + lines_z2, labels_z1 + labels_z2, loc="upper right", fontsize=9)

    # -- Row 3: Phase + group delay (dual y-axis) --
    phase_deg = df["phase_deg"].values
    phase_rad = np.radians(phase_deg)
    phase_unwrapped = np.degrees(np.unwrap(phase_rad))

    color_phase = plot_theme.COLORS["primary"]
    ax_ph.plot(freq, phase_unwrapped, color=color_phase, linewidth=1.4, label="Phase")
    ax_ph.set_ylabel("Phase (degrees)", color=color_phase)
    ax_ph.tick_params(axis="y", labelcolor=color_phase)

    ax_ph2 = ax_ph.twinx()
    if len(freq) > 1:
        dphi = np.diff(np.unwrap(phase_rad))
        df_vals = np.diff(freq)
        tau_g_ms = -dphi / (2 * np.pi * df_vals) * 1000
        freq_mid = (freq[:-1] + freq[1:]) / 2

        color_gd = plot_theme.COLORS["tertiary"]
        ax_ph2.plot(freq_mid, tau_g_ms, color=color_gd, linestyle="--", linewidth=1.2, label="Group delay")
        ax_ph2.set_ylabel("Group Delay (ms)", color=color_gd)
        ax_ph2.tick_params(axis="y", labelcolor=color_gd)

    ax_ph.set_title("Phase Response & Group Delay")
    plot_theme.setup_grid(ax_ph)

    lines_ph1, labels_ph1 = ax_ph.get_legend_handles_labels()
    lines_ph2, labels_ph2 = ax_ph2.get_legend_handles_labels()
    ax_ph.legend(lines_ph1 + lines_ph2, labels_ph1 + labels_ph2, loc="upper right", fontsize=9)

    # Shared frequency axis on bottom panel only
    plot_theme.setup_freq_axis(ax_ph, freq.min(), freq.max())
    # Apply log scale to upper panels (shared x) but hide their tick labels
    ax_spl.set_xlabel("")
    ax_z.set_xlabel("")

    plot_theme.save_figure(fig, output_file)
    print(f"Dashboard saved to {output_file}")


def main():
    """CLI for combined dashboard generation."""
    parser = argparse.ArgumentParser(
        description="Generate combined SPL/impedance/phase dashboard from solver CSV.",
    )
    parser.add_argument("csv_file", type=str, help="Input CSV with frequency, spl, z_real, z_imag, phase_deg.")
    parser.add_argument("output_file", type=str, help="Output image path.")
    args = parser.parse_args()
    generate_dashboard(args.csv_file, args.output_file)


if __name__ == "__main__":
    main()
