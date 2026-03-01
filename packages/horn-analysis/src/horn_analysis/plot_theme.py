"""Shared plot theme for professional acoustic-style plots.

Provides consistent styling across all horn-analysis plotting modules:
- Muted professional color palette
- Standard audio frequency axis ticks (20, 50, 100, 200, ..., 20k)
- Fixed dB divisions on SPL axes
- Dual-weight grid (major + minor)
- Sans-serif fonts, 150 DPI output
"""

import base64
import io

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


# -- Color palette ------------------------------------------------------------

COLORS = {
    "primary": "#1f4e79",     # steel blue
    "secondary": "#c0392b",   # muted red
    "tertiary": "#27864e",    # forest green
    "quaternary": "#7b4ea3",  # muted purple
    "quinary": "#d4831a",     # burnt orange
}

MULTI_COLORS = [
    "#1f4e79",
    "#c0392b",
    "#27864e",
    "#7b4ea3",
    "#d4831a",
    "#2a9d8f",
    "#6c5b7b",
    "#c97b3d",
]

PROFILE_STYLES = {
    "conical":     {"color": "#1f4e79", "linestyle": "-"},
    "exponential": {"color": "#27864e", "linestyle": "--"},
    "hyperbolic":  {"color": "#d4831a", "linestyle": "-."},
    "tractrix":    {"color": "#7b4ea3", "linestyle": ":"},
    "os":          {"color": "#0d9488", "linestyle": "-"},
    "lecleach":    {"color": "#be185d", "linestyle": "--"},
    "cd":          {"color": "#475569", "linestyle": "-."},
}

_DEFAULT_PROFILE_STYLE = {"color": "#6b7280", "linestyle": "-"}


def profile_style(profile: str) -> dict:
    """Return matplotlib line style dict for a horn profile name."""
    return PROFILE_STYLES.get(profile, _DEFAULT_PROFILE_STYLE)


# -- Standard audio frequency ticks -------------------------------------------

_AUDIO_TICKS = [20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000]


def _freq_label(f: float) -> str:
    """Format a frequency value: use 'k' suffix for >= 1000 Hz."""
    if f >= 1000:
        k = f / 1000
        return f"{k:g}k"
    return f"{f:g}"


def setup_freq_axis(ax, f_min: float, f_max: float):
    """Configure a log-scale frequency x-axis with standard audio ticks.

    Filters the standard audio tick set to the [f_min, f_max] data range
    and applies minor ticks between major ticks.
    """
    ax.set_xscale("log")
    ax.set_xlabel("Frequency (Hz)")

    visible = [f for f in _AUDIO_TICKS if f_min <= f <= f_max]
    if not visible:
        visible = _AUDIO_TICKS

    ax.set_xticks(visible)
    ax.set_xticklabels([_freq_label(f) for f in visible])
    ax.set_xlim(f_min, f_max)
    ax.xaxis.set_minor_locator(ticker.LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())


def setup_spl_axis(ax, spl_data, division: float = 5):
    """Configure SPL y-axis with fixed dB divisions.

    Auto-ranges from the data, snapping to the nearest ``division`` boundary.
    """
    spl_min = np.min(spl_data)
    spl_max = np.max(spl_data)
    y_lo = np.floor(spl_min / division) * division
    y_hi = np.ceil(spl_max / division) * division
    # Ensure at least one division of padding
    if y_hi - y_lo < division * 2:
        y_lo -= division
        y_hi += division
    ax.set_ylim(y_lo, y_hi)
    ax.set_yticks(np.arange(y_lo, y_hi + division, division))
    ax.set_ylabel("SPL (dB)")


# -- Grid styling --------------------------------------------------------------

def setup_grid(ax):
    """Apply dual-weight grid: heavier major lines, lighter minor lines."""
    ax.set_axisbelow(True)
    ax.grid(True, which="major", color="#cccccc", linewidth=0.6)
    ax.grid(True, which="minor", color="#e8e8e8", linewidth=0.3)


# -- Figure creation / saving --------------------------------------------------

def apply_theme():
    """Set matplotlib rcParams for a professional acoustic look."""
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 10,
        "figure.dpi": 150,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.spines.top": False,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
    })


def create_figure(figsize=(10, 6), nrows=1, ncols=1, **kwargs):
    """Create a themed figure and axes.

    Returns (fig, ax) for single-panel or (fig, axes) for multi-panel.
    """
    apply_theme()
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, **kwargs)
    return fig, axes


def save_figure(fig, path: str):
    """Apply tight_layout and save at 150 DPI."""
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_to_b64(fig) -> str:
    """Render a matplotlib figure to a base64 data-URI and close it."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    encoded = base64.b64encode(buf.read()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


# -- Annotation helpers --------------------------------------------------------

def annotate_f3(ax, f3_low: float, f3_high: float):
    """Draw vertical dashed lines at the -3 dB cutoff frequencies."""
    for f3 in (f3_low, f3_high):
        ax.axvline(f3, color="#888888", linestyle=":", linewidth=0.8, alpha=0.7)


def target_band_span(ax, target, alpha: float = 0.08):
    """Add a shaded vertical span for the target frequency band.

    ``target`` must have ``.f_low_hz`` and ``.f_high_hz`` attributes.
    """
    ax.axvspan(target.f_low_hz, target.f_high_hz, color="#6366f1", alpha=alpha, label="Target band")
