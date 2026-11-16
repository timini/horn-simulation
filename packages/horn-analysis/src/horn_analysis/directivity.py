"""Experimental infinite-baffle directivity from modal_baffled solver output.

The pattern comes from the 16 axisymmetric mouth-velocity modes that the
modal_baffled solve already exports (horn_core.modal_radiation). It assumes an
infinite rigid baffle and the far field, and has not been compared with a
measured or BEM/PML polar. It is relative to the axis, so it does not depend
on the driver: every mode is scaled by the same driver transfer.
"""
import numpy as np
import pandas as pd

from horn_core.modal_radiation import modal_baffled_directivity_db

DIRECTIVITY_LABEL = "experimental, infinite baffle"
DEFAULT_ANGLES_DEG = np.arange(5., 91., 5.)
MODE_COUNT = 16


def modal_velocity_matrix(frame):
    """Complex modal mouth velocities (frequencies x 16), or None for other models."""
    if "radiation_model" not in frame or not (frame.radiation_model == "modal_baffled").any():
        return None
    if not (frame.radiation_model == "modal_baffled").all() or "modal_mode_count" not in frame \
            or not (frame.modal_mode_count == MODE_COUNT).all():
        raise ValueError("Inconsistent modal aperture contract")
    columns = [f"modal_v_{n}_{part}" for n in range(MODE_COUNT) for part in ("real", "imag")]
    if any(name not in frame for name in columns):
        raise ValueError("Missing modal aperture velocity coefficients")
    velocity = np.column_stack([frame[f"modal_v_{n}_real"].to_numpy() + 1j * frame[f"modal_v_{n}_imag"].to_numpy()
                                for n in range(MODE_COUNT)])
    if not np.isfinite(velocity).all():
        raise ValueError("Non-finite modal aperture velocity")
    return velocity


def directivity_table(frame, angles_deg=DEFAULT_ANGLES_DEG):
    """Long table of frequency, theta_deg and spl_db relative to the axis.

    Returns None when the response was not solved with modal_baffled.
    """
    velocity = modal_velocity_matrix(frame)
    if velocity is None:
        return None
    from horn_core.duct import DEFAULT_AIR
    rho = float(frame.air_rho_kg_m3.iloc[0]) if "air_rho_kg_m3" in frame else DEFAULT_AIR.rho
    c = float(frame.air_c_m_s.iloc[0]) if "air_c_m_s" in frame else DEFAULT_AIR.c
    radius = float(np.sqrt(frame.mouth_area_m2.iloc[0] / np.pi))
    angles = np.asarray(angles_deg, dtype=float)
    f = frame.frequency.to_numpy()
    levels = modal_baffled_directivity_db(f, velocity, radius, angles, rho=rho, c=c)
    theta = np.r_[0., angles]
    levels = np.column_stack([np.zeros(len(f)), levels])
    return pd.DataFrame({"frequency": np.repeat(f, len(theta)), "theta_deg": np.tile(theta, len(f)),
                         "spl_db": levels.ravel(), "model": "modal_baffled_far_field",
                         "status": DIRECTIVITY_LABEL})


def beamwidth_figure(tables, target=None, threshold_db=-6.):
    """Beamwidth against frequency for several (label, table) pairs."""
    from horn_analysis import plot_theme
    from horn_analysis.directivity_plot import compute_beamwidth
    plot_theme.apply_theme()
    fig, ax = plot_theme.create_figure(figsize=(11, 5))
    low = high = None
    for index, (label, table) in enumerate(tables):
        f, width = compute_beamwidth(table, threshold_db=threshold_db)
        ax.plot(f, width, label=label, linewidth=1.6,
                color=plot_theme.MULTI_COLORS[index % len(plot_theme.MULTI_COLORS)])
        low = f.min() if low is None else min(low, f.min())
        high = f.max() if high is None else max(high, f.max())
    if target is not None:
        plot_theme.target_band_span(ax, target)
    if low is not None:
        plot_theme.setup_freq_axis(ax, low, high)
    ax.set_ylim(0, 185)
    ax.set_ylabel(f"Beamwidth, {threshold_db:g} dB (degrees)")
    ax.set_title(f"Beamwidth ({DIRECTIVITY_LABEL})")
    plot_theme.setup_grid(ax)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    return fig
