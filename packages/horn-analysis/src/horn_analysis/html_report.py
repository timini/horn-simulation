"""Self-contained HTML report generator for auto-mode pipeline results.

Produces a single HTML file with base64-embedded matplotlib plots,
ranking tables, and driver T-S parameter tables. No external
dependencies beyond matplotlib/numpy/pandas (already required).
"""

import html
from datetime import datetime, timezone
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd

from horn_core.parameters import DriverParameters
from horn_analysis.scoring import TargetSpec
from horn_analysis import plot_theme
from horn_analysis.horn_render import fig_to_b64_3d


# -- Profile badge styles (HTML-only; matplotlib styles come from plot_theme) --

_BADGE_STYLES = {
    "conical":     {"badge_bg": "#dbeafe", "badge_fg": "#1e40af"},
    "exponential": {"badge_bg": "#dcfce7", "badge_fg": "#166534"},
    "hyperbolic":  {"badge_bg": "#fef3c7", "badge_fg": "#92400e"},
    "tractrix":    {"badge_bg": "#ede9fe", "badge_fg": "#5b21b6"},
    "os":          {"badge_bg": "#ccfbf1", "badge_fg": "#115e59"},
    "lecleach":    {"badge_bg": "#fce7f3", "badge_fg": "#9d174d"},
    "cd":          {"badge_bg": "#e2e8f0", "badge_fg": "#334155"},
}

_DEFAULT_BADGE = {"badge_bg": "#f3f4f6", "badge_fg": "#374151"}


def _badge_for(profile: str) -> dict:
    return _BADGE_STYLES.get(profile, _DEFAULT_BADGE)


# -- Helpers ----------------------------------------------------------------

def _fmt(value, fmt: str = ".1f", fallback: str = "\u2014") -> str:
    """Safely format an Optional numeric value."""
    if value is None:
        return fallback
    try:
        return f"{value:{fmt}}"
    except (ValueError, TypeError):
        return fallback


# -- Plot generators --------------------------------------------------------

def _plot_coupled_spl_comparison(
    csv_pairs: List[Tuple[str, str]],
    target: TargetSpec,
) -> str:
    """Overlaid coupled SPL for top N candidates with target band."""
    fig, ax = plot_theme.create_figure(figsize=(11, 5.5))

    all_freq = []
    for csv_path, label in csv_pairs:
        df = pd.read_csv(csv_path)
        # Extract profile name from label (last word in parens)
        profile = label.rsplit("(", 1)[-1].rstrip(")") if "(" in label else ""
        style = plot_theme.profile_style(profile)
        ax.plot(df["frequency"], df["spl"], label=label,
                color=style["color"], linestyle=style["linestyle"], linewidth=1.4)
        all_freq.extend(df["frequency"].values)

    plot_theme.target_band_span(ax, target)

    if all_freq:
        plot_theme.setup_freq_axis(ax, min(all_freq), max(all_freq))
    ax.set_ylabel("SPL (dB)")
    ax.set_title("Coupled SPL \u2014 Top Candidates")
    plot_theme.setup_grid(ax)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    return plot_theme.fig_to_b64(fig)


def _plot_raw_profile_spl(
    solver_csvs: Dict[str, str],
    target: TargetSpec,
) -> str:
    """Overlaid uncoupled SPL comparing horn profiles."""
    fig, ax = plot_theme.create_figure(figsize=(11, 5.5))

    all_freq = []
    for profile, csv_path in sorted(solver_csvs.items()):
        df = pd.read_csv(csv_path)
        style = plot_theme.profile_style(profile)
        ax.plot(df["frequency"], df["spl"], label=profile.capitalize(),
                color=style["color"], linestyle=style["linestyle"], linewidth=1.4)
        all_freq.extend(df["frequency"].values)

    plot_theme.target_band_span(ax, target)

    if all_freq:
        plot_theme.setup_freq_axis(ax, min(all_freq), max(all_freq))
    ax.set_ylabel("SPL (dB)")
    ax.set_title("Raw Horn SPL by Profile (uncoupled)")
    plot_theme.setup_grid(ax)
    ax.legend(fontsize=9)
    fig.tight_layout()
    return plot_theme.fig_to_b64(fig)


def _plot_profile_impedance(solver_csvs: Dict[str, str]) -> str:
    """|Z| magnitude + phase angle overlay for each profile (dual y-axis)."""
    fig, ax1 = plot_theme.create_figure(figsize=(11, 5.5))
    ax2 = ax1.twinx()

    lines = []
    all_freq = []
    for profile, csv_path in sorted(solver_csvs.items()):
        df = pd.read_csv(csv_path)
        freq = df["frequency"].values
        z_complex = df["z_real"].values + 1j * df["z_imag"].values
        z_mag = np.abs(z_complex)
        z_angle = np.degrees(np.angle(z_complex))

        style = plot_theme.profile_style(profile)
        l1, = ax1.plot(freq, z_mag, color=style["color"],
                       linestyle=style["linestyle"], linewidth=1.3,
                       label=f"|Z| {profile}")
        l2, = ax2.plot(freq, z_angle, color=style["color"],
                       linestyle=":", linewidth=1.0, alpha=0.7,
                       label=f"\u2220Z {profile}")
        lines.extend([l1, l2])
        all_freq.extend(freq)

    if all_freq:
        plot_theme.setup_freq_axis(ax1, min(all_freq), max(all_freq))
    ax1.set_ylabel("|Z| (Pa\u00b7s/m)")
    ax2.set_ylabel("\u2220Z (degrees)")
    ax1.set_title("Throat Impedance by Profile")
    plot_theme.setup_grid(ax1)
    ax1.legend(handles=lines, fontsize=7, loc="upper right", ncol=2)
    fig.tight_layout()
    return plot_theme.fig_to_b64(fig)


def _plot_profile_phase(solver_csvs: Dict[str, str]) -> str:
    """Unwrapped phase + group delay subplots for each profile."""
    fig, (ax1, ax2) = plot_theme.create_figure(figsize=(11, 8), nrows=2, ncols=1, sharex=True)

    all_freq = []
    for profile, csv_path in sorted(solver_csvs.items()):
        df = pd.read_csv(csv_path)
        freq = df["frequency"].values
        phase_deg = df["phase_deg"].values

        phase_rad = np.radians(phase_deg)
        phase_unwrapped = np.degrees(np.unwrap(phase_rad))

        style = plot_theme.profile_style(profile)
        ax1.plot(freq, phase_unwrapped, color=style["color"],
                 linestyle=style["linestyle"], linewidth=1.3,
                 label=profile.capitalize())

        # Group delay
        if len(freq) > 1:
            dphi = np.diff(np.unwrap(phase_rad))
            df_vals = np.diff(freq)
            tau_g_ms = -dphi / (2 * np.pi * df_vals) * 1000
            freq_mid = (freq[:-1] + freq[1:]) / 2
            ax2.plot(freq_mid, tau_g_ms, color=style["color"],
                     linestyle=style["linestyle"], linewidth=1.0,
                     label=profile.capitalize())

        all_freq.extend(freq)

    if all_freq:
        f_min, f_max = min(all_freq), max(all_freq)
        plot_theme.setup_freq_axis(ax1, f_min, f_max)
        ax1.set_xlabel("")
        plot_theme.setup_freq_axis(ax2, f_min, f_max)

    ax1.set_ylabel("Phase (degrees)")
    ax1.set_title("Unwrapped Phase Response")
    plot_theme.setup_grid(ax1)
    ax1.legend(fontsize=9)

    ax2.set_ylabel("Group Delay (ms)")
    plot_theme.setup_grid(ax2)
    ax2.legend(fontsize=9)

    fig.tight_layout()
    return plot_theme.fig_to_b64(fig)


# -- Table renderers --------------------------------------------------------

def _profile_badge(profile: str) -> str:
    badge = _badge_for(profile)
    return (
        f'<span style="display:inline-block;padding:2px 8px;border-radius:4px;'
        f'font-size:0.8em;font-weight:600;'
        f'background:{badge["badge_bg"]};color:{badge["badge_fg"]}">'
        f'{html.escape(profile.capitalize())}</span>'
    )


def _render_rankings_rows(
    ranked_results: List[dict],
    drivers: Dict[str, DriverParameters],
    show_geometry: bool = False,
) -> str:
    rows = []
    for rank, r in enumerate(ranked_results, 1):
        kpi = r.get("kpi", {})
        f3l = _fmt(kpi.get("f3_low_hz"), ".0f")
        f3h = _fmt(kpi.get("f3_high_hz"), ".0f")
        drv = drivers.get(r.get("driver_id", ""))
        drv_type = html.escape(drv.driver_type or "—") if drv else "—"
        drv_size = html.escape(drv.nominal_diameter or "—") if drv else "—"
        drv_power = _fmt(getattr(drv, "power_w", None), ".0f") if drv else "—"

        geom_cols = ""
        if show_geometry:
            geom_cols = (
                f"<td>{_fmt(r.get('mouth_radius'), '.4f')}</td>"
                f"<td>{_fmt(r.get('length'), '.4f')}</td>"
            )

        rows.append(
            f"<tr>"
            f"<td>{rank}</td>"
            f"<td>{html.escape(r.get('manufacturer', ''))}</td>"
            f"<td>{html.escape(r.get('model_name', ''))}</td>"
            f"<td>{drv_type}</td>"
            f"<td>{drv_size}</td>"
            f"<td>{drv_power}</td>"
            f"<td>{_profile_badge(r.get('horn_label', ''))}</td>"
            f"{geom_cols}"
            f"<td><strong>{_fmt(r.get('composite_score'), '.3f')}</strong></td>"
            f"<td>{_fmt(r.get('bandwidth_coverage'), '.1%')}</td>"
            f"<td>{_fmt(r.get('passband_ripple_db'), '.1f')}</td>"
            f"<td>{_fmt(r.get('avg_sensitivity_db'), '.1f')}</td>"
            f"<td>{f3l} \u2014 {f3h}</td>"
            f"<td>{_fmt(kpi.get('peak_spl_db'), '.1f')}</td>"
            f"</tr>"
        )
    return "\n".join(rows)


def _render_drivers_rows(drivers: Dict[str, DriverParameters]) -> str:
    rows = []
    for drv in sorted(drivers.values(), key=lambda d: d.driver_id):
        sd_cm2 = drv.sd_m2 * 1e4
        mms_g = drv.mms_kg * 1e3
        le_mh = drv.le_h * 1e3
        ebp = drv.fs_hz / drv.qes if drv.qes and drv.qes > 0 else None
        rows.append(
            f"<tr>"
            f"<td>{html.escape(drv.manufacturer)}</td>"
            f"<td>{html.escape(drv.model_name)}</td>"
            f"<td>{_fmt(drv.fs_hz, '.0f')}</td>"
            f"<td>{_fmt(drv.re_ohm, '.1f')}</td>"
            f"<td>{_fmt(drv.bl_tm, '.1f')}</td>"
            f"<td>{_fmt(sd_cm2, '.2f')}</td>"
            f"<td>{_fmt(mms_g, '.2f')}</td>"
            f"<td>{_fmt(le_mh, '.2f')}</td>"
            f"<td>{_fmt(drv.qms, '.1f')}</td>"
            f"<td>{_fmt(drv.qes, '.2f')}</td>"
            f"<td>{_fmt(drv.qts, '.2f')}</td>"
            f"<td>{_fmt(ebp, '.0f')}</td>"
            f"</tr>"
        )
    return "\n".join(rows)


# -- HTML template ----------------------------------------------------------

_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Horn Auto-Select Report</title>
<style>
  *, *::before, *::after {{ box-sizing: border-box; }}
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
    margin: 0; padding: 0; background: #f8fafc; color: #1e293b;
  }}
  .container {{ max-width: 1100px; margin: 0 auto; padding: 24px 20px; }}
  h1 {{ font-size: 1.6em; margin: 0 0 4px; }}
  h2 {{ font-size: 1.2em; margin: 32px 0 12px; border-bottom: 2px solid #e2e8f0; padding-bottom: 6px; }}
  .subtitle {{ color: #64748b; font-size: 0.9em; margin-bottom: 20px; }}
  /* Summary cards */
  .cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 12px; margin-bottom: 24px; }}
  .card {{
    background: #fff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 14px 16px;
    text-align: center;
  }}
  .card .label {{ font-size: 0.75em; text-transform: uppercase; color: #64748b; letter-spacing: 0.05em; }}
  .card .value {{ font-size: 1.5em; font-weight: 700; color: #0f172a; margin-top: 4px; }}
  /* Tables */
  table {{ width: 100%; border-collapse: collapse; font-size: 0.85em; background: #fff; border-radius: 8px; overflow: hidden; }}
  th {{ background: #f1f5f9; text-align: left; padding: 10px 12px; font-weight: 600; white-space: nowrap; }}
  td {{ padding: 8px 12px; border-top: 1px solid #e2e8f0; }}
  tr:hover {{ background: #f8fafc; }}
  /* Plots */
  .plot {{ text-align: center; margin: 16px 0; }}
  .plot img {{ max-width: 100%; height: auto; border-radius: 6px; border: 1px solid #e2e8f0; }}
  .plot-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }}
  @media (max-width: 800px) {{ .plot-grid {{ grid-template-columns: 1fr; }} }}
  .design-summary {{ background: #fff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 16px 20px; margin-bottom: 24px; }}
  .design-summary dt {{ font-weight: 600; color: #475569; font-size: 0.85em; }}
  .design-summary dd {{ margin: 0 0 10px 0; font-size: 0.95em; }}
  .footer {{ margin-top: 40px; padding-top: 16px; border-top: 1px solid #e2e8f0; color: #94a3b8; font-size: 0.8em; text-align: center; }}
</style>
</head>
<body>
<div class="container">

<h1>Horn Auto-Select Report</h1>
<p class="subtitle">
  Target: {target_low:.0f} Hz — {target_high:.0f} Hz &nbsp;|&nbsp;
  Throat: {throat_radius:.4f} m &nbsp;|&nbsp;
  Mouth: {mouth_radius} &nbsp;|&nbsp;
  Length: {horn_length} &nbsp;|&nbsp;
  Profiles: {profiles} &nbsp;|&nbsp;
  Generated: {timestamp}
</p>

<div class="cards">
  <div class="card"><div class="label">Candidates scored</div><div class="value">{n_scored}</div></div>
  <div class="card"><div class="label">Top shown</div><div class="value">{n_top}</div></div>
  <div class="card"><div class="label">Best score</div><div class="value">{best_score}</div></div>
  <div class="card"><div class="label">Best BW coverage</div><div class="value">{best_bw}</div></div>
  <div class="card"><div class="label">Best sensitivity</div><div class="value">{best_sens}</div></div>
  <div class="card"><div class="label">Lowest ripple</div><div class="value">{best_ripple}</div></div>
</div>

{design_summary_section}

{geometry_section}

<h2>Rankings</h2>
<table>
<thead>
<tr>
  <th>#</th><th>Manufacturer</th><th>Model</th><th>Type</th><th>Size</th><th>Power (W)</th><th>Profile</th>
  {geometry_header_cols}
  <th>Score</th>
  <th>BW Cov.</th><th>Ripple (dB)</th><th>Sensitivity (dB)</th><th>f3 range (Hz)</th><th>Peak (dB)</th>
</tr>
</thead>
<tbody>
{rankings_rows}
</tbody>
</table>

<h2>Coupled SPL — Top Candidates</h2>
<div class="plot"><img src="{plot_coupled_spl}" alt="Coupled SPL comparison"></div>

<h2>Raw Horn SPL by Profile</h2>
<div class="plot"><img src="{plot_raw_spl}" alt="Raw profile SPL comparison"></div>

<h2>Impedance &amp; Phase</h2>
<div class="plot-grid">
  <div class="plot"><img src="{plot_impedance}" alt="Throat impedance"></div>
  <div class="plot"><img src="{plot_phase}" alt="Phase response"></div>
</div>

<h2>Pre-Screened Drivers</h2>
<table>
<thead>
<tr>
  <th>Manufacturer</th><th>Model</th><th>Fs (Hz)</th><th>Re (&Omega;)</th><th>BL (T&middot;m)</th>
  <th>Sd (cm&sup2;)</th><th>Mms (g)</th><th>Le (mH)</th><th>Qms</th><th>Qes</th><th>Qts</th><th>EBP</th>
</tr>
</thead>
<tbody>
{drivers_rows}
</tbody>
</table>

<div class="footer">Horn Auto-Select Report &mdash; generated {timestamp}</div>

</div>
</body>
</html>
"""


# -- Public entry point -----------------------------------------------------

def _render_design_summary(derived_geometry: dict) -> str:
    """Render the Design Summary section for fullauto mode."""
    mr = derived_geometry.get("mouth_radius_range", [])
    lr = derived_geometry.get("length_range", [])
    sr = derived_geometry.get("sim_freq_range", [])
    return (
        '<h2>Design Summary</h2>\n'
        '<div class="design-summary"><dl>'
        f'<dt>Ideal mouth radius</dt><dd>{_fmt(derived_geometry.get("ideal_mouth_radius"), ".4f")} m '
        f'(from c&#8320;/2&pi;f<sub>low</sub>)</dd>'
        f'<dt>Mouth radius range</dt><dd>{_fmt(mr[0] if mr else None, ".4f")} — '
        f'{_fmt(mr[1] if len(mr) > 1 else None, ".4f")} m (&plusmn;30%)</dd>'
        f'<dt>Length range</dt><dd>{_fmt(lr[0] if lr else None, ".4f")} — '
        f'{_fmt(lr[1] if len(lr) > 1 else None, ".4f")} m '
        f'(&lambda;/4 to &lambda;/2)</dd>'
        f'<dt>Simulation freq range</dt><dd>{_fmt(sr[0] if sr else None, ".0f")} — '
        f'{_fmt(sr[1] if len(sr) > 1 else None, ".0f")} Hz (&plusmn;0.5 octave)</dd>'
        f'<dt>Candidate count</dt><dd>{derived_geometry.get("candidate_count", "—")}</dd>'
        '</dl></div>'
    )


def generate_html_report(
    all_ranked: List[dict],
    solver_csvs: Dict[str, str],
    drivers: Dict[str, DriverParameters],
    throat_radius: float,
    target: TargetSpec,
    csv_pairs: List[Tuple[str, str]],
    top_n: int = 5,
    mouth_radius: Optional[float] = None,
    length: Optional[float] = None,
    derived_geometry: Optional[dict] = None,
) -> str:
    """Generate a self-contained HTML report string.

    Args:
        all_ranked: Ranked results (sorted by composite_score descending).
        solver_csvs: Mapping of profile name -> solver CSV path.
        drivers: Mapping of driver_id -> DriverParameters.
        throat_radius: Horn throat radius in metres.
        target: Target frequency specification.
        csv_pairs: List of (csv_path, label) for coupled SPL plots.
        top_n: Number of top candidates to include.
        mouth_radius: Horn mouth radius in metres (enables 3D geometry renders).
        length: Horn length in metres (enables 3D geometry renders).
        derived_geometry: Optional dict from geometry_designer (fullauto mode).

    Returns:
        Complete HTML document as a string.
    """
    top_results = all_ranked[:top_n]
    show_geometry = derived_geometry is not None

    # Summary card values
    best_score = _fmt(top_results[0]["composite_score"], ".3f") if top_results else "\u2014"
    best_bw = _fmt(
        max((r.get("bandwidth_coverage", 0) for r in top_results), default=None), ".1%"
    ) if top_results else "\u2014"
    best_sens = _fmt(
        max((r.get("avg_sensitivity_db", 0) for r in top_results), default=None), ".1f"
    ) if top_results else "\u2014"
    best_ripple = _fmt(
        min((r.get("passband_ripple_db", 99) for r in top_results), default=None), ".1f"
    ) if top_results else "\u2014"

    # Generate 3D horn geometry renders (when dimensions are provided)
    geometry_html = ""
    if mouth_radius is not None and length is not None:
        geom_imgs = []
        for profile in sorted(solver_csvs.keys()):
            b64 = fig_to_b64_3d(
                throat_radius=throat_radius,
                mouth_radius=mouth_radius,
                length=length,
                profile=profile,
                figsize=(12, 5),
            )
            geom_imgs.append(
                f'<div class="plot"><img src="{b64}" '
                f'alt="{html.escape(profile.capitalize())} horn geometry"></div>'
            )
        geometry_html = "\n".join(geom_imgs)

    # Generate plots
    plot_coupled_spl = _plot_coupled_spl_comparison(csv_pairs, target)
    plot_raw_spl = _plot_raw_profile_spl(solver_csvs, target)
    plot_impedance = _plot_profile_impedance(solver_csvs)
    plot_phase = _plot_profile_phase(solver_csvs)

    # Render tables
    rankings_rows = _render_rankings_rows(top_results, drivers, show_geometry=show_geometry)
    drivers_rows = _render_drivers_rows(drivers)

    # Conditional sections for fullauto
    design_summary_section = _render_design_summary(derived_geometry) if show_geometry else ""
    geometry_header_cols = '<th>Mouth R (m)</th><th>Length (m)</th>' if show_geometry else ""

    # Mouth/Length display: for fullauto show "varies", for auto show fixed value
    if show_geometry:
        mouth_display = "varies"
        length_display = "varies"
    else:
        mouth_display = f"{mouth_radius:.4f} m" if mouth_radius is not None else "—"
        length_display = f"{length:.3f} m" if length is not None else "—"

    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    # Build geometry section HTML
    if geometry_html:
        geometry_section = (
            '<h2>Horn Geometry</h2>\n'
            '<div class="plot-grid">\n'
            f'{geometry_html}\n'
            '</div>'
        )
    else:
        geometry_section = ""

    return _HTML_TEMPLATE.format_map({
        "target_low": target.f_low_hz,
        "target_high": target.f_high_hz,
        "throat_radius": throat_radius,
        "mouth_radius": mouth_display,
        "horn_length": length_display,
        "profiles": ", ".join(sorted(solver_csvs.keys())),
        "timestamp": timestamp,
        "n_scored": len(all_ranked),
        "n_top": len(top_results),
        "best_score": best_score,
        "best_bw": best_bw,
        "best_sens": best_sens,
        "best_ripple": best_ripple,
        "design_summary_section": design_summary_section,
        "geometry_header_cols": geometry_header_cols,
        "geometry_section": geometry_section,
        "rankings_rows": rankings_rows,
        "plot_coupled_spl": plot_coupled_spl,
        "plot_raw_spl": plot_raw_spl,
        "plot_impedance": plot_impedance,
        "plot_phase": plot_phase,
        "drivers_rows": drivers_rows,
    })
