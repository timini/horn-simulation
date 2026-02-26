"""Phase response plotting: unwrapped phase and optional group delay vs frequency."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


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
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    else:
        fig, ax1 = plt.subplots(figsize=(10, 6))

    ax1.semilogx(freq, phase_unwrapped_deg, color="tab:blue")
    ax1.set_ylabel("Phase (degrees)")
    ax1.set_title("Phase Response vs Frequency")
    ax1.grid(True, which="both", ls="--", alpha=0.5)

    if group_delay:
        # Group delay: τ_g = -dφ/dω = -(1/(2π)) * dφ_rad/df
        if len(freq) > 1:
            dphi = np.diff(phase_unwrapped_rad)
            df_vals = np.diff(freq)
            tau_g = -dphi / (2 * np.pi * df_vals)
            freq_mid = (freq[:-1] + freq[1:]) / 2
            # Convert to milliseconds
            tau_g_ms = tau_g * 1000

            ax2.semilogx(freq_mid, tau_g_ms, color="tab:green")
            ax2.set_ylabel("Group Delay (ms)")
            ax2.set_xlabel("Frequency (Hz)")
            ax2.grid(True, which="both", ls="--", alpha=0.5)
        else:
            ax2.set_xlabel("Frequency (Hz)")
            ax2.text(0.5, 0.5, "Insufficient data for group delay",
                     ha="center", va="center", transform=ax2.transAxes)
    else:
        ax1.set_xlabel("Frequency (Hz)")

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
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
