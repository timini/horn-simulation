"""Akabak-style directivity / radiation pattern visualization.

Reads a directivity CSV produced by the solver (columns:
``frequency``, ``theta_deg``, ``spl_db``) and generates four
diagnostic plots:

1. **Polar directivity** — SPL vs angle at selected frequencies.
2. **Directivity contour** (sonogram) — SPL heat-map over frequency and angle.
3. **Beamwidth vs frequency** — coverage angle where SPL drops by a threshold.
4. **Directivity index vs frequency** — on-axis gain relative to omnidirectional.

All plots use the shared :mod:`horn_analysis.plot_theme` styling.
"""

import argparse
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from horn_analysis import plot_theme


# -- Data loading -------------------------------------------------------------

def load_directivity(csv_file: str) -> pd.DataFrame:
    """Load and validate a directivity CSV file.

    Expected columns: ``frequency``, ``theta_deg``, ``spl_db``.

    Returns
    -------
    pd.DataFrame
        Sorted by (frequency, theta_deg).
    """
    df = pd.read_csv(csv_file)
    required = {"frequency", "theta_deg", "spl_db"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Directivity CSV missing columns: {missing}")
    return df.sort_values(["frequency", "theta_deg"]).reset_index(drop=True)


# -- Polar directivity -------------------------------------------------------

def plot_polar_directivity(
    df: pd.DataFrame,
    frequencies: Optional[List[float]] = None,
    output_file: str = "polar_directivity.png",
    *,
    n_auto: int = 6,
):
    """Polar plot of SPL vs angle at selected frequencies.

    Parameters
    ----------
    df : DataFrame
        Directivity data with columns ``frequency``, ``theta_deg``, ``spl_db``.
    frequencies : list of float, optional
        Frequencies to plot. If *None*, auto-selects ``n_auto`` log-spaced
        values from the data range.
    output_file : str
        Output image path.
    n_auto : int
        Number of frequencies to auto-select when *frequencies* is None.
    """
    plot_theme.apply_theme()

    avail_freqs = np.sort(df["frequency"].unique())
    if frequencies is None:
        indices = np.round(np.linspace(0, len(avail_freqs) - 1, n_auto)).astype(int)
        frequencies = avail_freqs[indices]
    else:
        # Snap each requested frequency to the nearest available
        frequencies = [avail_freqs[np.argmin(np.abs(avail_freqs - f))] for f in frequencies]

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(8, 8))
    ax.set_theta_zero_location("N")  # 0 degrees at top (on-axis)
    ax.set_theta_direction(-1)
    ax.set_thetamin(0)
    ax.set_thetamax(180)

    colors = plot_theme.MULTI_COLORS
    for i, freq in enumerate(frequencies):
        sub = df[df["frequency"] == freq].sort_values("theta_deg")
        theta_rad = np.radians(sub["theta_deg"].values)
        spl = sub["spl_db"].values
        label = f"{freq:.0f} Hz" if freq < 1000 else f"{freq / 1000:.1f} kHz"
        ax.plot(theta_rad, spl, color=colors[i % len(colors)], linewidth=1.3, label=label)

    ax.set_title("Polar Directivity", pad=20)
    ax.legend(loc="lower left", fontsize=7, bbox_to_anchor=(1.05, 0))

    plot_theme.save_figure(fig, output_file)


# -- Directivity contour (sonogram) ------------------------------------------

def plot_directivity_contour(
    df: pd.DataFrame,
    output_file: str = "directivity_contour.png",
    *,
    db_range: float = 40,
):
    """Frequency-angle heat-map of SPL (Akabak-style sonogram).

    Parameters
    ----------
    df : DataFrame
        Directivity data.
    output_file : str
        Output image path.
    db_range : float
        Dynamic range below peak to display (dB).
    """
    plot_theme.apply_theme()

    freqs = np.sort(df["frequency"].unique())
    angles = np.sort(df["theta_deg"].unique())

    spl_grid = np.full((len(angles), len(freqs)), np.nan)
    for j, f in enumerate(freqs):
        sub = df[df["frequency"] == f].sort_values("theta_deg")
        for i, a in enumerate(angles):
            row = sub[sub["theta_deg"] == a]
            if not row.empty:
                spl_grid[i, j] = row["spl_db"].values[0]

    vmax = np.nanmax(spl_grid)
    vmin = vmax - db_range

    fig, ax = plot_theme.create_figure(figsize=(11, 6))
    mesh = ax.pcolormesh(
        freqs, angles, spl_grid,
        cmap="inferno",
        vmin=vmin,
        vmax=vmax,
        shading="gouraud",
    )
    plot_theme.setup_freq_axis(ax, freqs.min(), freqs.max())
    ax.set_ylabel("Angle (degrees)")
    ax.set_title("Directivity Contour")
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label("SPL (dB)")
    plot_theme.setup_grid(ax)

    plot_theme.save_figure(fig, output_file)


# -- Beamwidth ----------------------------------------------------------------

def compute_beamwidth(
    df: pd.DataFrame,
    threshold_db: float = -6,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute coverage angle vs frequency.

    The beamwidth is twice the angle at which the SPL drops by
    ``threshold_db`` relative to the on-axis (0 degrees) value.

    Parameters
    ----------
    df : DataFrame
        Directivity data.
    threshold_db : float
        Drop from on-axis SPL (negative, e.g. -6).

    Returns
    -------
    (frequencies, beamwidth_deg) : tuple of 1-D arrays
    """
    freqs = np.sort(df["frequency"].unique())
    beamwidths = np.full(len(freqs), np.nan)

    for i, f in enumerate(freqs):
        sub = df[df["frequency"] == f].sort_values("theta_deg")
        theta = sub["theta_deg"].values
        spl = sub["spl_db"].values

        # On-axis SPL (theta=0 or closest)
        idx_on = np.argmin(np.abs(theta))
        spl_on = spl[idx_on]
        target = spl_on + threshold_db  # threshold_db is negative

        # Find first angle where SPL drops below target
        below = np.where(spl < target)[0]
        if len(below) > 0:
            idx_cross = below[0]
            if idx_cross > 0:
                # Linear interpolation between adjacent points
                t0, t1 = theta[idx_cross - 1], theta[idx_cross]
                s0, s1 = spl[idx_cross - 1], spl[idx_cross]
                if s1 != s0:
                    angle = t0 + (target - s0) * (t1 - t0) / (s1 - s0)
                else:
                    angle = t0
                beamwidths[i] = 2 * angle  # Full coverage angle
            else:
                beamwidths[i] = 0.0
        else:
            beamwidths[i] = 2 * theta.max()  # Never drops below threshold

    return freqs, beamwidths


def plot_beamwidth(
    df: pd.DataFrame,
    output_file: str = "beamwidth.png",
    *,
    threshold_db: float = -6,
):
    """Plot beamwidth (coverage angle) vs frequency.

    Parameters
    ----------
    df : DataFrame
        Directivity data.
    output_file : str
        Output image path.
    threshold_db : float
        dB threshold for beamwidth calculation.
    """
    plot_theme.apply_theme()

    freqs, bw = compute_beamwidth(df, threshold_db=threshold_db)

    fig, ax = plot_theme.create_figure(figsize=(10, 5))
    ax.plot(freqs, bw, color=plot_theme.COLORS["primary"], linewidth=1.4)
    plot_theme.setup_freq_axis(ax, freqs.min(), freqs.max())
    ax.set_ylabel("Beamwidth (degrees)")
    ax.set_title(f"Beamwidth ({threshold_db} dB) vs Frequency")
    plot_theme.setup_grid(ax)

    plot_theme.save_figure(fig, output_file)


# -- Directivity index -------------------------------------------------------

def compute_directivity_index(
    df: pd.DataFrame,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute axisymmetric directivity index vs frequency.

    DI = 10 * log10( p_on_axis^2 / <p^2> )

    where <p^2> is the power-weighted angular average:
        <p^2> = integral( p^2 * sin(theta) dtheta ) / integral( sin(theta) dtheta )

    Returns
    -------
    (frequencies, DI_dB) : tuple of 1-D arrays
    """
    freqs = np.sort(df["frequency"].unique())
    di = np.full(len(freqs), np.nan)

    for i, f in enumerate(freqs):
        sub = df[df["frequency"] == f].sort_values("theta_deg")
        theta_deg = sub["theta_deg"].values
        spl = sub["spl_db"].values

        # Convert SPL to linear pressure squared (relative)
        p_sq = 10 ** (spl / 10)

        theta_rad = np.radians(theta_deg)

        # On-axis value (closest to 0 degrees)
        idx_on = np.argmin(np.abs(theta_deg))
        p_sq_on = p_sq[idx_on]

        # Numerical integration using trapezoidal rule
        integrand = p_sq * np.sin(theta_rad)
        numerator = np.trapezoid(integrand, theta_rad)
        denominator = np.trapezoid(np.sin(theta_rad), theta_rad)

        if denominator > 0 and numerator > 0:
            p_sq_avg = numerator / denominator
            di[i] = 10 * np.log10(p_sq_on / p_sq_avg)

    return freqs, di


def plot_directivity_index(
    df: pd.DataFrame,
    output_file: str = "directivity_index.png",
):
    """Plot directivity index (DI) vs frequency.

    Parameters
    ----------
    df : DataFrame
        Directivity data.
    output_file : str
        Output image path.
    """
    plot_theme.apply_theme()

    freqs, di = compute_directivity_index(df)

    fig, ax = plot_theme.create_figure(figsize=(10, 5))
    ax.plot(freqs, di, color=plot_theme.COLORS["primary"], linewidth=1.4)
    plot_theme.setup_freq_axis(ax, freqs.min(), freqs.max())
    ax.set_ylabel("Directivity Index (dB)")
    ax.set_title("Directivity Index vs Frequency")
    plot_theme.setup_grid(ax)

    plot_theme.save_figure(fig, output_file)


# -- Combined report ---------------------------------------------------------

def generate_directivity_report(csv_file: str, output_dir: str):
    """Generate all four directivity plots into a directory.

    Creates:
    - ``polar_directivity.png``
    - ``directivity_contour.png``
    - ``beamwidth.png``
    - ``directivity_index.png``

    Parameters
    ----------
    csv_file : str
        Path to the directivity CSV.
    output_dir : str
        Directory for output images (created if needed).
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    df = load_directivity(csv_file)

    plot_polar_directivity(df, output_file=str(out / "polar_directivity.png"))
    plot_directivity_contour(df, output_file=str(out / "directivity_contour.png"))
    plot_beamwidth(df, output_file=str(out / "beamwidth.png"))
    plot_directivity_index(df, output_file=str(out / "directivity_index.png"))


# -- CLI ----------------------------------------------------------------------

def main():
    """Command-line interface: ``horn-directivity-plot``."""
    parser = argparse.ArgumentParser(
        description="Generate Akabak-style directivity plots from a CSV."
    )
    parser.add_argument("csv_file", help="Directivity CSV (frequency, theta_deg, spl_db)")
    parser.add_argument("--output-dir", default="directivity",
                        help="Output directory (default: directivity/)")
    args = parser.parse_args()

    generate_directivity_report(args.csv_file, args.output_dir)


if __name__ == "__main__":
    main()
