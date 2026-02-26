"""Throat impedance plotting: |Z| and angle(Z) vs frequency on dual y-axes."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plot_impedance(csv_file: str, output_file: str):
    """Plot throat impedance magnitude and phase from solver CSV output.

    Expects CSV columns: frequency, z_real, z_imag.

    Args:
        csv_file: Input CSV path.
        output_file: Output image path.
    """
    df = pd.read_csv(csv_file)
    freq = df["frequency"].values
    z_complex = df["z_real"].values + 1j * df["z_imag"].values

    z_mag = np.abs(z_complex)
    z_angle = np.degrees(np.angle(z_complex))

    fig, ax1 = plt.subplots(figsize=(10, 6))

    color_mag = "tab:blue"
    ax1.set_xlabel("Frequency (Hz)")
    ax1.set_ylabel("|Z| (Pa·s/m)", color=color_mag)
    ax1.semilogx(freq, z_mag, color=color_mag, label="|Z|")
    ax1.tick_params(axis="y", labelcolor=color_mag)

    ax2 = ax1.twinx()
    color_phase = "tab:red"
    ax2.set_ylabel("∠Z (degrees)", color=color_phase)
    ax2.semilogx(freq, z_angle, color=color_phase, linestyle="--", label="∠Z")
    ax2.tick_params(axis="y", labelcolor=color_phase)

    ax1.set_title("Throat Impedance vs Frequency")
    ax1.grid(True, which="both", ls="--", alpha=0.5)

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right")

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f"Impedance plot saved to {output_file}")


def main():
    """CLI for impedance plotting."""
    parser = argparse.ArgumentParser(description="Plot throat impedance from solver CSV.")
    parser.add_argument("csv_file", type=str, help="Input CSV with frequency, z_real, z_imag columns.")
    parser.add_argument("output_file", type=str, help="Output image path.")
    args = parser.parse_args()
    plot_impedance(args.csv_file, args.output_file)


if __name__ == "__main__":
    main()
