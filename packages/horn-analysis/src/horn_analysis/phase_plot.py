"""Phase response plotting: unwrapped phase and optional group delay vs frequency."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from horn_analysis import plot_theme


def plot_phase(csv_file: str, output_file: str, group_delay: bool = False):
    """Plot phase response from solver CSV output.

    Expects CSV columns: frequency, phase_deg.

    Args:
        csv_file: Input CSV path.
        output_file: Output image path.
        group_delay: If True, add a second subplot with group delay.
    """
    df = pd.read_csv(csv_file)
    freq = df["frequency"].values
    phase_deg = df["phase_deg"].values

    # Unwrap phase for continuous display
    phase_rad = np.radians(phase_deg)
    phase_unwrapped_rad = np.unwrap(phase_rad)
    phase_unwrapped_deg = np.degrees(phase_unwrapped_rad)

    if group_delay:
        fig, (ax1, ax2) = plot_theme.create_figure(figsize=(10, 8), nrows=2, ncols=1, sharex=True)
    else:
        fig, ax1 = plot_theme.create_figure(figsize=(10, 6))

    ax1.plot(freq, phase_unwrapped_deg, color=plot_theme.COLORS["primary"], linewidth=1.4)
    ax1.set_ylabel("Phase (degrees)")
    ax1.set_title("Phase Response vs Frequency")

    f_min, f_max = freq.min(), freq.max()

    if group_delay:
        plot_theme.setup_freq_axis(ax1, f_min, f_max)
        # Remove x-label from top panel since bottom panel will have it
        ax1.set_xlabel("")
        plot_theme.setup_grid(ax1)

        # Group delay: τ_g = -dφ/dω = -(1/(2π)) * dφ_rad/df
        if len(freq) > 1:
            dphi = np.diff(phase_unwrapped_rad)
            df_vals = np.diff(freq)
            tau_g = -dphi / (2 * np.pi * df_vals)
            freq_mid = (freq[:-1] + freq[1:]) / 2
            # Convert to milliseconds
            tau_g_ms = tau_g * 1000

            ax2.plot(freq_mid, tau_g_ms, color=plot_theme.COLORS["tertiary"], linewidth=1.4)
            ax2.set_ylabel("Group Delay (ms)")
            plot_theme.setup_freq_axis(ax2, f_min, f_max)
            plot_theme.setup_grid(ax2)
        else:
            ax2.text(0.5, 0.5, "Insufficient data for group delay",
                     ha="center", va="center", transform=ax2.transAxes)
            ax2.set_xlabel("Frequency (Hz)")
    else:
        plot_theme.setup_freq_axis(ax1, f_min, f_max)
        plot_theme.setup_grid(ax1)

    plot_theme.save_figure(fig, output_file)
    print(f"Phase plot saved to {output_file}")


def main():
    """CLI for phase response plotting."""
    parser = argparse.ArgumentParser(description="Plot phase response from solver CSV.")
    parser.add_argument("csv_file", type=str, help="Input CSV with frequency, phase_deg columns.")
    parser.add_argument("output_file", type=str, help="Output image path.")
    parser.add_argument("--group-delay", action="store_true", help="Include group delay subplot.")
    args = parser.parse_args()
    plot_phase(args.csv_file, args.output_file, group_delay=args.group_delay)


if __name__ == "__main__":
    main()
