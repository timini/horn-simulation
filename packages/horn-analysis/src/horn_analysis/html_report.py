"""Self-contained HTML report generator for auto-mode pipeline results.

Produces a single HTML file with base64-embedded matplotlib plots,
ranking tables, and driver T-S parameter tables. No external
dependencies beyond matplotlib/numpy/pandas (already required).
"""

from horn_analysis.kpi import format_kpi

import html
from datetime import datetime, timezone
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd

from horn_core.parameters import DriverParameters
from horn_analysis.scoring import TargetSpec, describe_scoring, size_compromise_note
from horn_analysis import plot_theme
from horn_analysis.horn_render import fig_to_b64_3d


# -- Profile badge styles (HTML-only; matplotlib styles come from plot_theme) --

_BADGE_STYLES = {
    "conical":     {"badge_bg": "#dbeafe", "badge_fg": "#1e40af"},
    "exponential": {"badge_bg": "#dcfce7", "badge_fg": "#166534"},
    "hyperbolic":  {"badge_bg": "#fef3c7", "badge_fg": "#92400e"},
    "tractrix":    {"badge_bg": "#ede9fe", "badge_fg": "#5b21b6"},
    "os":          {"badge_bg": "#ccfbf1", "badge_fg": "#115e59"},
    "scaled_tractrix": {"badge_bg": "#fce7f3", "badge_fg": "#9d174d"},
    "lecleach":    {"badge_bg": "#fce7f3", "badge_fg": "#9d174d"},
    "cd":          {"badge_bg": "#e2e8f0", "badge_fg": "#334155"},
}

_DEFAULT_BADGE = {"badge_bg": "#f3f4f6", "badge_fg": "#374151"}


def _badge_for(profile: str) -> dict:
    return _BADGE_STYLES.get(profile, _DEFAULT_BADGE)


# -- Helpers ----------------------------------------------------------------

def _fmt(value, fmt: str = ".1f", fallback: str = "n/a") -> str:
    """Safely format an Optional numeric value."""
    if value is None:
        return fallback
    try:
        return f"{value:{fmt}}"
    except (ValueError, TypeError):
        return fallback


# -- Plot generators --------------------------------------------------------

_COMPARISON_LINESTYLES = ["-", "--", "-.", ":", (0, (3, 1, 1, 1))]


def _plot_coupled_spl_comparison(
    csv_pairs: List[Tuple[str, str]],
    target: TargetSpec,
) -> str:
    """Overlaid coupled SPL for top N candidates with target band."""
    fig, ax = plot_theme.create_figure(figsize=(11, 5.5))

    all_freq = []
    for i, (csv_path, label) in enumerate(csv_pairs):
        df = pd.read_csv(csv_path)
        color = plot_theme.MULTI_COLORS[i % len(plot_theme.MULTI_COLORS)]
        ls = _COMPARISON_LINESTYLES[i % len(_COMPARISON_LINESTYLES)]
        ax.plot(df["frequency"], df["spl"], label=label,
                color=color, linestyle=ls, linewidth=1.8)
        all_freq.extend(df["frequency"].values)

    plot_theme.target_band_span(ax, target)

    if all_freq:
        plot_theme.setup_freq_axis(ax, min(all_freq), max(all_freq))
        ax.legend(fontsize=8, loc="best")
    else:
        _no_data(ax, "No feasible design in evaluated candidates")
    ax.set_ylabel("SPL (dB)")
    ax.set_title("Estimated observer output: Top Candidates")
    plot_theme.setup_grid(ax)
    fig.tight_layout()
    return plot_theme.fig_to_b64(fig)


_SHORT_TO_PROFILE = {
    "con": "conical",
    "exp": "exponential",
    "hyp": "hyperbolic",
    "tra": "tractrix",
    "os": "os",
    "lec": "lecleach",
    "sca": "scaled_tractrix",
    "cd": "cd",
}


def _pick_representative_per_profile(solver_csvs: Dict[str, str]) -> Dict[str, str]:
    """Select one representative candidate per horn profile.

    Groups candidates by profile (extracted from IDs like ``auto_hyp_0012``)
    and picks the median candidate from each group.
    """
    groups: Dict[str, list] = {}
    for cid, path in solver_csvs.items():
        parts = cid.split("_")
        short = parts[1] if len(parts) >= 2 else cid
        profile = _SHORT_TO_PROFILE.get(short, short)
        groups.setdefault(profile, []).append((cid, path))

    result: Dict[str, str] = {}
    for profile, items in sorted(groups.items()):
        items.sort()
        result[profile] = items[len(items) // 2][1]
    return result


def _horn_curves(top_results: List[dict], solver_csvs: Dict[str, str]) -> List[Tuple[str, str]]:
    """(legend, solver CSV) for each distinct horn among the ranked rows, best rank first."""
    from horn_analysis.report_tables import horn_legend, ranked_horns
    return [(horn_legend(rank, row), solver_csvs[row["horn_label"]])
            for rank, row in ranked_horns(top_results) if row["horn_label"] in solver_csvs]


def _curve_style(i):
    return {"color": plot_theme.MULTI_COLORS[i % len(plot_theme.MULTI_COLORS)],
            "linestyle": _COMPARISON_LINESTYLES[i % len(_COMPARISON_LINESTYLES)]}


def _no_data(ax, text="No ranked candidate"):
    ax.text(.5, .5, text, ha="center", va="center", transform=ax.transAxes)


def _plot_raw_profile_spl(
    curves: List[Tuple[str, str]],
    target: TargetSpec,
) -> str:
    """Horn-only level at unit throat pressure for each ranked horn."""
    fig, ax = plot_theme.create_figure(figsize=(11, 5.5))

    all_freq = []
    for i, (label, csv_path) in enumerate(curves):
        df = pd.read_csv(csv_path)
        ax.plot(df["frequency"], df["spl"], label=label, linewidth=1.4, **_curve_style(i))
        all_freq.extend(df["frequency"].values)

    plot_theme.target_band_span(ax, target)

    if all_freq:
        plot_theme.setup_freq_axis(ax, min(all_freq), max(all_freq))
        ax.legend(fontsize=9)
    else:
        _no_data(ax)
    ax.set_ylabel("SPL (dB)")
    ax.set_title("Horn only at 1 Pa throat pressure, ranked horns (uncoupled)")
    plot_theme.setup_grid(ax)
    fig.tight_layout()
    return plot_theme.fig_to_b64(fig)


def _plot_profile_impedance(curves: List[Tuple[str, str]]) -> str:
    """Acoustic throat impedance |Z| and phase for each ranked horn (dual y-axis)."""
    fig, ax1 = plot_theme.create_figure(figsize=(11, 5.5))
    ax2 = ax1.twinx()

    lines = []
    all_freq = []
    for i, (label, csv_path) in enumerate(curves):
        df = pd.read_csv(csv_path)
        freq = df["frequency"].values
        z_complex = df["z_real"].values + 1j * df["z_imag"].values
        style = _curve_style(i)
        l1, = ax1.plot(freq, np.abs(z_complex), color=style["color"],
                       linestyle=style["linestyle"], linewidth=1.3, label=label)
        ax2.plot(freq, np.degrees(np.angle(z_complex)), color=style["color"],
                 linestyle=":", linewidth=1.0, alpha=0.7)
        lines.append(l1)
        all_freq.extend(freq)

    if all_freq:
        plot_theme.setup_freq_axis(ax1, min(all_freq), max(all_freq))
        ax1.legend(handles=lines, fontsize=7, loc="upper right")
    else:
        _no_data(ax1)
    ax1.set_ylabel("|Z| (Pa\u00b7s/m), solid")
    ax2.set_ylabel("\u2220Z (degrees), dotted")
    ax1.set_title("Acoustic throat impedance of the ranked horns (not electrical)")
    plot_theme.setup_grid(ax1)
    fig.tight_layout()
    return plot_theme.fig_to_b64(fig)


def _plot_profile_phase(curves: List[Tuple[str, str]]) -> str:
    """Horn-only unwrapped phase and group delay for each ranked horn."""
    fig, (ax1, ax2) = plot_theme.create_figure(figsize=(11, 8), nrows=2, ncols=1, sharex=True)

    all_freq = []
    for i, (label, csv_path) in enumerate(curves):
        df = pd.read_csv(csv_path)
        freq = df["frequency"].values
        if "phase_deg" not in df:
            continue
        phase_rad = np.unwrap(np.radians(df["phase_deg"].values))
        style = _curve_style(i)
        ax1.plot(freq, np.degrees(phase_rad), linewidth=1.3, label=label, **style)
        if len(freq) > 1:
            tau_g_ms = -np.diff(phase_rad) / (2 * np.pi * np.diff(freq)) * 1000
            ax2.plot((freq[:-1] + freq[1:]) / 2, tau_g_ms, linewidth=1.0, **style)
        all_freq.extend(freq)

    if all_freq:
        f_min, f_max = min(all_freq), max(all_freq)
        plot_theme.setup_freq_axis(ax1, f_min, f_max)
        ax1.set_xlabel("")
        plot_theme.setup_freq_axis(ax2, f_min, f_max)
        ax1.legend(fontsize=9)
    else:
        _no_data(ax1)

    ax1.set_ylabel("Phase (degrees)")
    ax1.set_title("Horn-only phase at the mouth plane, ranked horns")
    plot_theme.setup_grid(ax1)
    ax2.set_ylabel("Horn-only group delay (ms)")
    plot_theme.setup_grid(ax2)

    fig.tight_layout()
    return plot_theme.fig_to_b64(fig)


def _plot_electrical(electrical_tables) -> str:
    from horn_analysis.electrical_output import electrical_figure
    return plot_theme.fig_to_b64(electrical_figure(electrical_tables))


# -- Table renderers --------------------------------------------------------

def _profile_badge(profile: str) -> str:
    badge = _badge_for(profile)
    return (
        f'<span style="display:inline-block;padding:2px 8px;border-radius:4px;'
        f'font-size:0.8em;font-weight:600;'
        f'background:{badge["badge_bg"]};color:{badge["badge_fg"]}">'
        f'{html.escape(profile.capitalize())}</span>'
    )


_EVIDENCE_LABELS = {
    "driver_usable_band_unknown": "Driver operating band is unknown",
    "driver_parameter_provenance_missing": "Driver parameters lack a verified source",
    "driver_interface_unverified": "Driver chamber or throat adapter has not been validated",
    "moving_mass_and_rear_load_not_separated": "Diaphragm mass and rear air load have not been characterized separately",
    "maximum_output_limits_unknown": "Power or excursion limits are missing",
    "mechanical_damping_unknown": "Mechanical damping is unknown",
    "front_chamber_volume_estimated": "Front chamber volume is a cone-geometry estimate, not a measured adapter",
    "front_chamber_volume_unknown": "Front chamber volume is unknown; the ideal area transformer is optimistic",
    "diaphragm_smaller_than_throat": "Throat is larger than the diaphragm",
    "front_chamber_not_lumped_in_band": "Target band goes above ka = 1 for the front chamber, where a uniform chamber pressure is approximate",
    "response_limited_by_estimated_front_chamber": "Meets the limits without a front chamber but not with the estimated one; supply the real adapter volume",
    "throat_above_first_asymmetric_mode": "Throat ka exceeds 1.84 at the top of the band, where an asymmetric duct mode can propagate",
}


def _render_evidence_summary(results):
    """One line per driver and distinct set of evidence gaps, in rank order."""
    lines = []
    for r in results:
        line = ("<li><strong>" + html.escape(r.get("model_name") or r.get("driver_id", "")) + "</strong>: "
                + html.escape("; ".join(_EVIDENCE_LABELS.get(g, g.replace("_", " ")) for g in r.get("evidence_gaps", []))
                              or "Independent assembly validation is pending") + ".</li>")
        if line not in lines:
            lines.append(line)
    return "<ul>" + "".join(lines) + "</ul>"


def _chamber_note(row) -> str:
    """Mark a level computed with an estimated front chamber (front_chamber_estimated)."""
    if not row.get("front_chamber_estimated"):
        return ""
    ideal = (row.get("without_front_chamber") or {}).get("passband_ripple_db")
    return ('<br><small title="front_chamber_estimated">est. chamber'
            + (f"; {_fmt(ideal, '.1f')} without" if ideal is not None else "") + '</small>')


_STATUS_LABELS = {
    "model_feasible": "model feasible",
    "insufficient_evidence": "feasible, evidence gaps",
    "infeasible": "rejected",
}


def _render_rankings_rows(
    ranked_results: List[dict],
    drivers: Dict[str, DriverParameters],
    show_geometry: bool = False,
    show_price: bool = False,
) -> str:
    from horn_analysis.report_tables import candidate_id
    rows = []
    for rank, r in enumerate(ranked_results, 1):
        kpi = r.get("kpi", {})
        # f3 relative to the target-band mean, like the ripple and mean columns.
        f3l = format_kpi(r, "f3_low_hz")
        f3h = format_kpi(r, "f3_high_hz")
        drv = drivers.get(r.get("driver_id", ""))
        drv_type = html.escape(drv.driver_type or "n/a") if drv else "n/a"
        drv_size = html.escape(drv.nominal_diameter or "n/a") if drv else "n/a"
        drv_power = _fmt(getattr(drv, "power_w", None), ".0f") if drv else "n/a"

        geom_cols = ""
        if show_geometry:
            mm = lambda key: _fmt(r[key] * 1e3 if r.get(key) is not None else None, ".1f")
            geom_cols = f"<td>{mm('throat_radius')}</td><td>{mm('mouth_radius')}</td><td>{mm('length')}</td>"

        interface = r.get("driver_interface") or {}
        chamber = interface.get("front_chamber_source")
        chamber_cell = "n/a" if chamber is None else (
            f"{_fmt(interface.get('front_chamber_volume_cm3'), '.0f')} cm&sup3; {html.escape(chamber)}"
            + (" (flag)" if r.get("estimated_chamber_reasons") else ""))
        ident = candidate_id(r) if r.get("driver_id") and r.get("horn_label") else ""
        status = r.get("eligibility_status", "")
        rows.append(
            f'<tr data-rank="{rank}" data-candidate="{html.escape(ident)}">'
            f"<td>{rank}</td>"
            f"<td><code>{html.escape(ident)}</code></td>"
            f"<td>{html.escape(r.get('manufacturer', ''))}</td>"
            f"<td>{html.escape(r.get('model_name', ''))}</td>"
            f"<td>{drv_type}</td>"
            f"<td>{drv_size}</td>"
            f"<td>{drv_power}</td>"
            f"<td>{_profile_badge(r.get('profile') or r.get('horn_label', ''))}</td>"
            f"{geom_cols}"
            f"<td><strong>{_fmt(r.get('composite_score'), '.3f')}</strong>{' (near tie)' if r.get('comparison_status') == 'near_tie' else ''}</td>"
            f"<td>{'≥' if r.get('extension_is_lower_bound') else ''}{_fmt(r.get('extension_octaves'), '+.2f')}</td>"
            f"<td>{_flat_band_cell(r)}</td>"
            f"<td>{_fmt(r.get('passband_ripple_db'), '.1f')}{_chamber_note(r)}</td>"
            f"<td>{_fmt(r.get('avg_sensitivity_db'), '.1f')}</td>"
            f"<td>{f3l} to {f3h}</td>"
            f"<td>{_fmt(kpi.get('peak_spl_db'), '.1f')}</td>"
            f"<td>{chamber_cell}</td>"
            f"<td>{html.escape(_STATUS_LABELS.get(status, status.replace('_', ' ')))}</td>"
            f"{_price_cell(r) if show_price else ''}"
            f"</tr>"
        )
    return "\n".join(rows)


def _price_cell(row) -> str:
    price = row.get("price")
    if not price:
        return "<td>no current offer</td>"
    from horn_analysis.pricing import format_price
    link = f'<a href="{html.escape(price["url"])}">{html.escape(price["source_id"])}</a>'
    return (f"<td>{html.escape(format_price(price))}<br><small>{link}, {html.escape(price['retrieved_on'])}"
            + (", in stock" if price.get("in_stock") else ", out of stock" if price.get("in_stock") is False else "")
            + "</small>" + ("<br><strong>cheapest feasible</strong>" if row.get("cheapest_feasible") else "") + "</td>")


def _render_pricing(pricing: Optional[dict]) -> str:
    """Price region, snapshot dates, sources and the cheapest feasible design."""
    if not pricing:
        return ""
    from horn_analysis.pricing import describe_pricing
    lines = describe_pricing(pricing)
    items = "".join(f"<li>{html.escape(line.strip())}</li>" for line in lines)
    return f'<div class="scoring"><h3>Prices</h3><ul>{items}</ul></div>'


def _flat_band_cell(r) -> str:
    """Widest +-3 dB flat band over the sweep, for example '≥1.45 (566-1543)'."""
    if "flat_band_3db_octaves" not in r:
        return "n/a"
    low, high = r.get("flat_band_3db_low_hz"), r.get("flat_band_3db_high_hz")
    if low is None or high is None:
        return "none"
    mark = "≥" if r.get("flat_band_3db_is_lower_bound") else ""
    return f"{mark}{r['flat_band_3db_octaves']:.2f} ({low:.0f} to {high:.0f} Hz)"


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
            f"<td>{_fmt(drv.power_w, '.0f')}</td>"
            f"<td>{_fmt(drv.sensitivity_db_1w1m, '.1f')}</td>"
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
  .table-wrap {{ overflow-x: auto; -webkit-overflow-scrolling: touch; margin: 0 -4px; padding: 0 4px; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.85em; background: #fff; border-radius: 8px; overflow: hidden; min-width: 800px; }}
  th {{ background: #f1f5f9; text-align: left; padding: 10px 12px; font-weight: 600; white-space: nowrap; }}
  td {{ padding: 8px 12px; border-top: 1px solid #e2e8f0; white-space: nowrap; }}
  tr:hover {{ background: #f8fafc; }}
  /* Plots */
  .plot {{ text-align: center; margin: 16px 0; }}
  .plot img {{ max-width: 100%; height: auto; border-radius: 6px; border: 1px solid #e2e8f0; }}
  .plot-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }}
  @media (max-width: 800px) {{ .plot-grid {{ grid-template-columns: 1fr; }} }}
  .design-summary {{ background: #fff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 16px 20px; margin-bottom: 24px; }}
  .design-summary dt {{ font-weight: 600; color: #475569; font-size: 0.85em; }}
  .design-summary dd {{ margin: 0 0 10px 0; font-size: 0.95em; }}
  .banner {{ background: #fef3c7; border: 2px solid #d97706; color: #78350f; border-radius: 8px; padding: 12px 16px; margin: 12px 0 20px; }}
  .banner strong {{ font-size: 1.1em; }}
  .footer {{ margin-top: 40px; padding-top: 16px; border-top: 1px solid #e2e8f0; color: #94a3b8; font-size: 0.8em; text-align: center; }}
</style>
</head>
<body>
<div class="container">

<h1>Horn Auto-Select Report</h1>
<div class="banner"><strong>Experimental prediction.</strong> {evidence_status}</div>
<div class="design-summary"><strong>{assessment_status}</strong>
{no_feasible_reason}
<p>Predictions are experimental. Driver suitability and the exterior model require independent validation.
{observer_summary}
Raw horn plots show mouth-plane pressure. Legacy inputs contain mouth-plane levels only.</p>
<p>Drive: {drive_voltage} V RMS. Observer: {observer_distance} m on axis from mouth plane.
Maximum target-band ripple: {ripple_limit} dB. Acoustic CAD describes air volume, not manufacturing walls.</p>
{stress_summary}
<p>Acoustic models used: {physics_summary}. Scores within 0.02 are near-ties for comparison; physical uncertainty remains unquantified.</p>
{radiation_caveat}
{evidence_summary}</div>
<p class="subtitle">
  Target: {target_low:.0f} to {target_high:.0f} Hz &nbsp;|&nbsp;
  Throat radius: {throat_display} &nbsp;|&nbsp;
  Mouth: {mouth_radius} &nbsp;|&nbsp;
  Length: {horn_length} &nbsp;|&nbsp;
  Horns: {profiles} &nbsp;|&nbsp;
  Generated: {timestamp}
</p>

<div class="cards">
  <div class="card"><div class="label">Candidates scored</div><div class="value">{n_scored}</div></div>
  <div class="card"><div class="label">Top shown</div><div class="value">{n_top}</div></div>
  <div class="card"><div class="label">Best score</div><div class="value">{best_score}</div></div>
  <div class="card"><div class="label">Best extension (oct)</div><div class="value">{best_bw}</div></div>
  <div class="card"><div class="label">Best mean output</div><div class="value">{best_sens}</div></div>
  <div class="card"><div class="label">Lowest ripple</div><div class="value">{best_ripple}</div></div>
</div>

{how_to_read}

{builder_section}

{design_summary_section}

{geometry_section}

<h2>Rankings</h2>
{scoring_section}
{pricing_section}
<div class="table-wrap">
<table>
<thead>
<tr>
  <th>#</th><th>Candidate ID (horn + driver)</th><th>Manufacturer</th><th>Model</th><th>Type</th><th>Size</th><th>Power (W)</th><th>Profile</th>
  {geometry_header_cols}
  <th>Score</th>
  <th>Extension (oct; ≥ = lower bound)</th><th>&plusmn;3 dB flat band (oct, Hz; ≥ = lower bound)</th><th>Ripple (dB)</th><th>Mean output (dB)</th><th>f3 range (Hz; ≤/≥ = sweep bound)</th><th>Peak (dB)</th><th>Front chamber</th><th>Evidence status</th>{price_header}
</tr>
</thead>
<tbody>
{rankings_rows}
</tbody>
</table>
</div>

{rejection_section}

<h2>Estimated observer output: top candidates</h2>
<div class="plot"><img src="{plot_coupled_spl}" alt="Coupled SPL comparison"></div>

{electrical_section}

<h2>Horn only, ranked horns</h2>
<div class="plot"><img src="{plot_raw_spl}" alt="Horn-only SPL of the ranked horns"></div>

<h2>Horn-only acoustic throat impedance &amp; phase</h2>
<div class="plot-grid">
  <div class="plot"><img src="{plot_impedance}" alt="Throat impedance"></div>
  <div class="plot"><img src="{plot_phase}" alt="Phase response"></div>
</div>

{directivity_section}
<h2>Pre-Screened Drivers</h2>
<div class="table-wrap">
<table>
<thead>
<tr>
  <th>Manufacturer</th><th>Model</th><th>Fs (Hz)</th><th>Re (&Omega;)</th><th>BL (T&middot;m)</th>
  <th>Sd (cm&sup2;)</th><th>Mms (g)</th><th>Le (mH)</th><th>Qms</th><th>Qes</th><th>Qts</th><th>EBP</th>
  <th>Power (W)</th><th>Sensitivity (dB, 1 W, 1 m)</th>
</tr>
</thead>
<tbody>
{drivers_rows}
</tbody>
</table>
</div>

<div class="footer">Horn Auto-Select Report, generated {timestamp}</div>

</div>
</body>
</html>
"""


# -- Public entry point -----------------------------------------------------

_BOUND_SOURCE_LABELS = {
    "user_fixed": "set by user",
    "heuristic_with_user_limit": "heuristic range with a user limit",
    "heuristic": "heuristic range",
    "driver_prescreen": "from driver prescreen",
}


def _render_scoring(scoring: dict, ranked: Optional[List[dict]] = None) -> str:
    """State the score formula, weights, exchange rates and bound origins."""
    weights = ", ".join(f"{html.escape(k)} {v:.3f}" for k, v in scoring["weights"].items())
    rates = scoring["exchange_rates"]
    rate_text = html.escape(rates["note"]) if "note" in rates else (
        f"1 dB less ripple is worth {rates['sensitivity_db_per_db_of_ripple']:.2f} dB of sensitivity; "
        f"1 octave more extension is worth {rates['sensitivity_db_per_octave_of_extension']:.2f} dB; "
        f"halving the horn volume is worth {rates['sensitivity_db_per_doubling_of_horn_volume']:.2f} dB.")
    bounds = "; ".join(
        f"{html.escape(k)}: {_BOUND_SOURCE_LABELS.get(v['source'], html.escape(v['source']))}"
        + (f" {_fmt(min(v['searched_range']), '.4f')} to {_fmt(max(v['searched_range']), '.4f')}" if v.get("searched_range") else "")
        for k, v in scoring["geometry_bounds"].items())
    constraints = "; ".join(
        f"{html.escape(k)} {html.escape(str(v['value']))} ({html.escape(v['source'].replace('_', ' '))})"
        for k, v in scoring["constraints"].items())
    extra = ""
    compromise = scoring.get("size_compromise")
    if compromise:
        extra += (f'<p><strong>Size compromise:</strong> {html.escape(size_compromise_note(compromise, ranked))}</p>')
    refinement = scoring.get("refinement")
    if refinement:
        widened = "; ".join(f'{html.escape(w["field"])} {html.escape(w["side"])} to {_fmt(w["to"][1] if w["side"] == "high" else w["to"][0], ".4f")}'
                            for w in refinement.get("widened_bounds", []))
        extra += (f'<p><strong>Refinement:</strong> {refinement["new_evaluations"]} new geometries from '
                  f'{html.escape(str(refinement["seed_candidate"]))} ({html.escape(refinement["seed_source"].replace("_", " "))})'
                  + (f'; widened heuristic bounds: {widened}' if widened else '') + '.</p>')
    if scoring.get("front_chamber_policy"):
        extra += f'<p><strong>Front chamber:</strong> {html.escape(scoring["front_chamber_policy"])}</p>'
    return (
        '<div class="design-summary">'
        f'<p><strong>Score.</strong> {html.escape(scoring["formula"])}</p>'
        f'<p><strong>Weights</strong> ({html.escape(scoring["weights_source"].replace("_", " "))}): {weights}. {rate_text}</p>'
        f'<p><strong>Constraints:</strong> {constraints}.</p>'
        f'<p><strong>Geometry bounds:</strong> {bounds}.</p>'
        + extra +
        '</div>'
    )


def _render_evidence_status(top_results: List[dict]) -> str:
    """Banner text: what the evidence status of the shown candidates is."""
    if not top_results:
        return "No design passed the model constraints. Nothing here is a recommendation."
    counts = {}
    for r in top_results:
        status = r.get("eligibility_status", "unknown")
        counts[status] = counts.get(status, 0) + 1
    parts = ", ".join(f"{n} {html.escape(_STATUS_LABELS.get(k, k.replace('_', ' ')))}" for k, n in counts.items())
    estimated = sum(1 for r in top_results if r.get("front_chamber_estimated"))
    text = (f"Evidence status of the {len(top_results)} shown: {parts}. None is physically validated; "
            "the driver, its adapter and the horn have not been measured together.")
    if estimated:
        text += (f" {estimated} use an estimated front chamber (<code>front_chamber_estimated</code>); "
                 "supply the real adapter volume to confirm them.")
    return text


def _render_how_to_read(target: TargetSpec) -> str:
    return (
        '<h2>How to read this report</h2>\n<div class="design-summary"><ol>'
        '<li>Every number is a linear model prediction (FEM horn, lumped driver). Treat it as an experimental '
        'starting point to build and measure, not a validated recommendation.</li>'
        f'<li>Rank 1 is the highest score among designs that meet the hard limits (ripple at most {target.max_ripple_db:g} dB '
        f'in {target.f_low_hz:.0f} to {target.f_high_hz:.0f} Hz, compression ratio, size and driver limits). '
        'A design whose limits are exceeded only because of an estimated front chamber is still ranked, '
        'with the exceeded limits listed in its Front chamber column (item 6). '
        'The score formula and weights are stated under Rankings. Scores within 0.02 are near ties.</li>'
        '<li>Each candidate is named <code>horn ID + driver ID</code>. Plot legends use the same <code>#rank ID</code>, '
        'so every curve maps to one ranking row.</li>'
        f'<li>Levels are on axis at {target.observation_distance_m:g} m for {target.voltage_rms:g} V RMS. '
        '|Ze| is the electrical impedance an amplifier sees, in ohms. The throat impedance plot is acoustic, horn only.</li>'
        '<li>Build from the top candidate section: the profile CSV gives the air passage wall radius against axial position '
        'in mm (z = 0 at the throat), and the STEP file is the same air volume at true size. Add your own wall thickness. '
        'File paths are relative to the folder that holds this report.</li>'
        '<li>The front chamber column states the chamber between diaphragm and throat that the prediction used: supplied, '
        'estimated from the cone geometry, or none. Estimated chambers are flagged and do not reject a design on their own.</li>'
        '<li>The rejection table counts why other designs failed and how close they came.</li>'
        '</ol></div>'
    )


def _render_builder(top_results: List[dict], builder_files: Dict[str, dict], drivers) -> str:
    """Winner dimensions in mm, its profile CSV and STEP path, and the chamber assumption."""
    if not top_results:
        return ""
    from horn_analysis.report_tables import candidate_label, chamber_assumption, has_parametric_geometry
    winner = top_results[0]
    files = builder_files.get(winner["horn_label"], {})
    mm = lambda key: _fmt(winner[key] * 1e3 if winner.get(key) is not None else None, ".1f")
    dims = (f'<dt>Throat</dt><dd>radius {mm("throat_radius")} mm</dd>'
            f'<dt>Mouth</dt><dd>radius {mm("mouth_radius")} mm</dd>'
            f'<dt>Length</dt><dd>{mm("length")} mm</dd>'
            f'<dt>Profile</dt><dd>{html.escape(str(winner.get("profile") or "imported"))}</dd>') if has_parametric_geometry(winner) else (
            f'<dt>Inlet area</dt><dd>{_fmt((winner.get("inlet_area_m2") or 0) * 1e6, ".1f")} mm&sup2;</dd>')
    step = files.get("step")
    profile_csv = files.get("profile_csv")
    rows = []
    for rank, r in enumerate(top_results, 1):
        f = builder_files.get(r["horn_label"], {})
        rows.append(f"<tr><td>{html.escape(candidate_label(rank, r))}</td>"
                    f"<td>{html.escape(chamber_assumption(r))}</td>"
                    f"<td><code>{html.escape(f.get('profile_csv') or 'n/a')}</code></td>"
                    f"<td><code>{html.escape(f.get('step') or 'n/a')}</code></td></tr>")
    return (
        '<h2>Top candidate: build data</h2>\n<div class="design-summary"><dl>'
        f'<dt>Candidate</dt><dd>{html.escape(candidate_label(1, winner))} '
        f'({html.escape(winner.get("manufacturer", ""))} {html.escape(winner.get("model_name", ""))})</dd>'
        f'{dims}'
        f'<dt>STEP file (acoustic air volume, mm; path relative to this report)</dt><dd><code>{html.escape(step or "not available")}</code></dd>'
        f'<dt>Profile CSV (z_mm, r_mm)</dt><dd><code>{html.escape(profile_csv or "not available")}</code></dd>'
        f'<dt>Front chamber assumed</dt><dd>{html.escape(chamber_assumption(winner))}</dd>'
        '</dl></div>\n<div class="table-wrap"><table><thead><tr><th>Candidate</th><th>Front chamber assumed</th>'
        '<th>Profile CSV</th><th>STEP</th></tr></thead><tbody>' + "".join(rows) + '</tbody></table></div>'
    )


def _render_rejections(summary: List[dict], n_rejected: int, no_feasible_reason: Optional[str], stage: Optional[str]) -> str:
    """Counts of rejection reasons with how far past each limit the designs were."""
    if not summary and not no_feasible_reason:
        return ""
    head = f'<h2>Rejected designs</h2>\n<div class="design-summary"><p>{n_rejected} driver-horn pairs were rejected'
    head += f" at the {html.escape(stage)} stage" if stage else ""
    head += ". A pair can fail more than one limit. Excess is how far past the limit it was; the closest miss is the smallest.</p>"
    if no_feasible_reason:
        head += f"<p><strong>Why no design is shown:</strong> {html.escape(no_feasible_reason)}</p>"
    head += "</div>"
    if not summary:
        return head
    rows = "".join(
        f'<tr data-reason="{html.escape(e["reason"])}"><td><code>{html.escape(e["reason"])}</code></td><td>{html.escape(e["description"])}</td>'
        f'<td>{e["count"]}</td><td>{_fmt(e["closest_excess"], ".2f")}</td><td>{_fmt(e["largest_excess"], ".2f")}</td>'
        f'<td>{html.escape(e["unit"]) or "n/a"}</td></tr>' for e in summary)
    return (head + '<div class="table-wrap"><table><thead><tr><th>Reason</th><th>Meaning</th><th>Count</th>'
            '<th>Closest miss</th><th>Largest excess</th><th>Unit</th></tr></thead><tbody>' + rows + '</tbody></table></div>')


def _sweep_margins(derived_geometry: dict) -> str:
    """Octaves swept below and above the target band, from the recorded values."""
    import math
    sr = derived_geometry.get("sim_freq_range") or []
    low, high = derived_geometry.get("target_f_low"), derived_geometry.get("target_f_high")
    if len(sr) != 2 or not low or not high or sr[0] <= 0:
        return ""
    return f" ({math.log2(low / sr[0]):.1f} octave below, {math.log2(sr[1] / high):.1f} above)"


def _render_design_summary(derived_geometry: dict) -> str:
    """Render the Design Summary section with optimization parameters."""
    tr = derived_geometry.get("throat_radius_range", [])
    throat_desc = " to ".join(_fmt(v, ".4f") for v in tr) + " m" if tr else "Unknown"
    mr = derived_geometry.get("mouth_radius_range", [])
    lr = derived_geometry.get("length_range", [])
    sr = derived_geometry.get("sim_freq_range", [])

    # Determine fixed vs optimized for mouth radius and length
    mr_fixed = len(mr) == 2 and abs(mr[1] - mr[0]) < 1e-6
    lr_fixed = len(lr) == 2 and abs(lr[1] - lr[0]) < 1e-6

    if mr_fixed:
        mouth_desc = f'{_fmt(mr[0], ".3f")} m <span style="color:#059669;font-weight:600">(fixed)</span>'
    else:
        mouth_desc = (
            f'{_fmt(mr[0] if mr else None, ".4f")} to '
            f'{_fmt(mr[1] if len(mr) > 1 else None, ".4f")} m '
            f'<span style="color:#2563eb;font-weight:600">(optimized)</span>'
        )

    if lr_fixed:
        length_desc = f'{_fmt(lr[0], ".3f")} m <span style="color:#059669;font-weight:600">(fixed)</span>'
    else:
        length_desc = (
            f'{_fmt(lr[0] if lr else None, ".4f")} to '
            f'{_fmt(lr[1] if len(lr) > 1 else None, ".4f")} m '
            f'<span style="color:#2563eb;font-weight:600">(optimized, &lambda;/4 to &lambda;/2)</span>'
        )

    return (
        '<h2>Design Summary</h2>\n'
        '<div class="design-summary"><dl>'
        f'<dt>Target frequency band</dt><dd>{_fmt(derived_geometry.get("target_f_low"), ".0f")} to '
        f'{_fmt(derived_geometry.get("target_f_high"), ".0f")} Hz</dd>'
        f'<dt>Throat radius range</dt><dd>{throat_desc}</dd>'
        f'<dt>Mouth radius</dt><dd>{mouth_desc}</dd>'
        f'<dt>Length</dt><dd>{length_desc}</dd>'
        f'<dt>Ideal mouth radius</dt><dd>{_fmt(derived_geometry.get("ideal_mouth_radius"), ".4f")} m '
        f'(from c&#8320;/2&pi;f<sub>low</sub>)</dd>'
        f'<dt>Simulation freq range</dt><dd>{_fmt(sr[0] if sr else None, ".0f")} to '
        f'{_fmt(sr[1] if len(sr) > 1 else None, ".0f")} Hz{_sweep_margins(derived_geometry)}</dd>'
        f'<dt>Profiles</dt><dd>7 (conical, exponential, hyperbolic, tractrix, OS, axially scaled tractrix, CD)</dd>'
        f'<dt>Geometry candidates</dt><dd>{derived_geometry.get("candidate_count", "n/a")}</dd>'
        '</dl></div>'
    )


def _embed_png(path: str) -> str:
    import base64
    with open(path, "rb") as handle:
        return "data:image/png;base64," + base64.b64encode(handle.read()).decode("ascii")


def _render_electrical_section(electrical_png: Optional[str]) -> str:
    if not electrical_png:
        return ""
    import base64
    with open(electrical_png, "rb") as handle:
        encoded = base64.b64encode(handle.read()).decode("ascii")
    return ('<h2>Electrical impedance, excursion and group delay</h2>\n'
            f'<div class="plot"><img src="data:image/png;base64,{encoded}" '
            'alt="Electrical impedance, excursion and group delay of the top candidates"></div>\n')


FLANGED_PISTON_CAVEAT = (
    "Known radiation-model error: the flanged_piston mouth boundary misses the independent MMM "
    "infinite-baffle reference by up to 1.05 dB and 5.7 degrees in throat impedance (limits 0.5 dB / 5 degrees), "
    "which can change driver coupling and ranking. modal_baffled passes that reference and is the default "
    "for generated axisymmetric horns.")


def _render_directivity(directivity_tables, target, top_results) -> str:
    if directivity_tables:
        from horn_analysis.directivity import beamwidth_figure, DIRECTIVITY_LABEL
        image = plot_theme.fig_to_b64(beamwidth_figure(directivity_tables, target))
        return ('<h2>Directivity (' + html.escape(DIRECTIVITY_LABEL) + ')</h2>\n'
                '<p>Far-field pattern from the solved axisymmetric mouth velocity, relative to the axis, '
                'for an infinite rigid baffle. A finite baffle or cabinet changes it by several dB, and it has '
                'not been compared with a measured polar. Use it to compare profiles, not as a coverage claim.</p>\n'
                f'<div class="plot"><img src="{image}" alt="Beamwidth against frequency"></div>\n')
    if top_results:
        return ('<h2>Directivity</h2>\n<p>Not available: off-axis pressure needs the modal_baffled '
                'radiation model.</p>\n')
    return ""



def _render_stress_summary(target, top_results) -> str:
    """State the test power and each shown candidate's excursion and power at it."""
    power = getattr(target, "test_power_w", None)
    if power is None:
        return ""
    lines = []
    for rank, result in enumerate(top_results, 1):
        if result.get("test_power_band_max_excursion_mm") is None:
            continue
        xmax = result.get("test_power_xmax_margin_mm")
        rating = result.get("test_power_rating_margin_w")
        lines.append(
            f"<li>#{rank} {html.escape(str(result.get('horn_label')))} + {html.escape(str(result.get('driver_id')))}: "
            f"{result['test_voltage_rms']:.2f} V RMS; in-band excursion {result['test_power_band_max_excursion_mm']:.2f} mm"
            + (f", Xmax margin {xmax:.2f} mm" if xmax is not None else ", Xmax unknown")
            + f"; in-band input power up to {result['test_power_band_max_input_power_w']:.1f} W"
            + (f", rating margin {rating:.1f} W" if rating is not None else ", rating unknown")
            + "</li>")
    return (f"<p>Stress check: {power:g} W into each driver's nominal impedance (Re when none is recorded). "
            "In-band excursion above Xmax or in-band input power above the driver rating at this level rejects a design. "
            "Linear model, no crossover high-pass; not a thermal or large-signal prediction.</p>"
            + (f"<ul>{''.join(lines)}</ul>" if lines else ""))

def generate_html_report(
    all_ranked: List[dict],
    solver_csvs: Dict[str, str],
    drivers: Dict[str, DriverParameters],
    throat_radius: float | None,
    target: TargetSpec,
    csv_pairs: List[Tuple[str, str]],
    top_n: int = 5,
    mouth_radius: Optional[float] = None,
    length: Optional[float] = None,
    derived_geometry: Optional[dict] = None,
    total_candidates: Optional[int] = None,
    total_scored: Optional[int] = None,
    lem_results: Optional[dict] = None,
    no_feasible_reason: Optional[str] = None,
    scoring: Optional[dict] = None,
    electrical_png: Optional[str] = None,
    directivity_tables: Optional[List[Tuple[str, "pd.DataFrame"]]] = None,
    electrical_tables: Optional[list] = None,
    rejected: Optional[List[dict]] = None,
    rejection_stage: Optional[str] = None,
    builder_files: Optional[Dict[str, dict]] = None,
    pricing: Optional[dict] = None,
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
        total_candidates: Total geometry candidates simulated.
        total_scored: Total driver-horn combinations scored (before top-N).
        electrical_png: Optional path to the electrical plot (|Ze|, excursion,
            group delay) of the top candidates; embedded when given.
        directivity_tables: (horn label, relative far-field table) pairs from
            horn_analysis.directivity, shown as an experimental section.
        electrical_tables: (legend, electrical_table frame, Xmax m) per ranked
            candidate, for the |Ze|, excursion and group delay plot.
        rejected: Rejected rows, summarised by reason with margins.
        rejection_stage: Where those rows were rejected (for example "FEM").
        builder_files: horn_label -> {"profile_csv": path, "step": path}.
        pricing: Run pricing summary (horn_analysis.pricing.annotate_ranking);
            adds a price column and a prices note. Price is not scored.

    Returns:
        Complete HTML document as a string.
    """
    top_results = all_ranked[:top_n]
    show_geometry = derived_geometry is not None

    # Summary card values
    best_score = _fmt(top_results[0]["composite_score"], ".3f") if top_results else "n/a"
    best_bw = _fmt(
        max((r.get("extension_octaves") for r in top_results if r.get("extension_octaves") is not None), default=None), "+.2f"
    ) if top_results else "n/a"
    best_sens = _fmt(
        max((r.get("avg_sensitivity_db", 0) for r in top_results), default=None), ".1f"
    ) if top_results else "n/a"
    best_ripple = _fmt(
        min((r.get("passband_ripple_db", 99) for r in top_results), default=None), ".1f"
    ) if top_results else "n/a"

    # Generate 3D horn geometry renders from top-ranked candidates
    geometry_html = ""
    if mouth_radius is not None and length is not None and not show_geometry:
        # Fixed geometry: render each profile
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
    elif top_results:
        # Variable geometry (auto/fullauto): render unique geometries from top results
        seen = set()
        geom_imgs = []
        for r in top_results:
            profile = r.get("profile", "") or r.get("horn_label", "")
            horn_label = r.get("horn_label", "")
            mr = r.get("mouth_radius")
            ln = r.get("length")
            tr = r.get("throat_radius", throat_radius)
            if profile and mr and ln and mr > tr:
                key = (profile, round(tr, 6), round(mr, 6), round(ln, 6))
                if key not in seen:
                    seen.add(key)
                    title = (
                        f"{profile.capitalize()} ({horn_label})"
                        f": throat {tr*1000:.1f}, mouth {mr*1000:.1f}, L {ln*1000:.0f} mm"
                    )
                    b64 = fig_to_b64_3d(
                        throat_radius=tr,
                        mouth_radius=mr,
                        length=ln,
                        profile=profile,
                        title=title,
                        figsize=(12, 5),
                    )
                    geom_imgs.append(
                        f'<div class="plot"><img src="{b64}" '
                        f'alt="{html.escape(title)}"></div>'
                    )
        geometry_html = "\n".join(geom_imgs)

    # Horn-only plots use the ranked horns, labelled with their best rank.
    horn_curves = _horn_curves(top_results, solver_csvs)

    # Generate plots
    plot_coupled_spl = _plot_coupled_spl_comparison(csv_pairs, target)
    plot_raw_spl = _plot_raw_profile_spl(horn_curves, target)
    plot_impedance = _plot_profile_impedance(horn_curves)
    plot_phase = _plot_profile_phase(horn_curves)
    electrical_section = ""
    if electrical_tables:
        # Embed the PNG written next to the report when there is one, so the
        # page shows the same figure as auto_electrical.png.
        image = _embed_png(electrical_png) if electrical_png else _plot_electrical(electrical_tables)
        electrical_section = (
            '<h2>Electrical impedance, excursion and group delay of the top candidates</h2>\n'
            '<p class="subtitle">|Ze| in ohms as the amplifier sees it; peak excursion against Xmax (dotted); '
            'system group delay on axis with the straight-line time of flight removed. '
            'Per-frequency values are in the coupled_NN CSV files.</p>\n'
            f'<div class="plot"><img src="{image}" alt="Electrical impedance, excursion and group delay"></div>'
        )
    elif electrical_png:
        electrical_section = _render_electrical_section(electrical_png)
    from horn_analysis.report_tables import rejection_summary
    rejected = rejected or []
    rejection_section = _render_rejections(rejection_summary(rejected, target, drivers), len(rejected),
                                           no_feasible_reason if not top_results else None, rejection_stage)

    # Render tables
    rankings_rows = _render_rankings_rows(top_results, drivers, show_geometry=show_geometry, show_price=bool(pricing))
    drivers_rows = _render_drivers_rows(drivers)

    # Conditional sections for fullauto
    design_summary_section = _render_design_summary(derived_geometry) if show_geometry else ""

    # LEM prescreen stats
    if lem_results:
        lem_total = lem_results.get("total_evaluated", "n/a")
        lem_pairs = lem_results.get("total_pairs", "n/a")
        lem_passed = len(lem_results.get("filtered_candidate_ids", []))
        design_summary_section += (
            '<h2>LEM Prescreening</h2>\n'
            '<div class="design-summary"><dl>'
            f'<dt>Candidates evaluated (LEM/Webster)</dt><dd>{lem_total}</dd>'
            f'<dt>Driver-horn pairs scored</dt><dd>{lem_pairs}</dd>'
            f'<dt>Candidates passed to FEM</dt><dd>{lem_passed}</dd>'
            '</dl></div>'
        )
    geometry_header_cols = '<th>Throat R (mm)</th><th>Mouth R (mm)</th><th>Length (mm)</th>' if show_geometry else ""

    radii = (derived_geometry or {}).get("throat_radius_range") or [r["throat_radius"] for r in all_ranked if r.get("throat_radius") is not None]
    if radii:
        throat_display = (f"{min(radii):.4f} to {max(radii):.4f} m (search range)"
                          if min(radii) != max(radii) else f"{radii[0]:.4f} m")
    else:
        throat_display = "not available" if show_geometry or throat_radius is None else f"{throat_radius:.4f} m"

    # Equal search endpoints describe a fixed dimension even in unified auto mode.
    if show_geometry:
        def dimension_display(key, fixed_value, precision):
            bounds = derived_geometry.get(key)
            if bounds:
                lo, hi = min(bounds), max(bounds)
                if lo == hi:
                    return f"{lo:.{precision}f} m"
                return f"{lo:.{precision}f} to {hi:.{precision}f} m (search range)"
            return f"{fixed_value:.{precision}f} m" if fixed_value is not None else "not available"

        mouth_display = dimension_display("mouth_radius_range", mouth_radius, 4)
        length_display = dimension_display("length_range", length, 3)
    else:
        mouth_display = f"{mouth_radius:.4f} m" if mouth_radius is not None else "n/a"
        length_display = f"{length:.3f} m" if length is not None else "n/a"

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

    n_scored = total_scored if total_scored is not None else len(all_ranked)

    return _HTML_TEMPLATE.format_map({
        "no_feasible_reason": "<p>"+html.escape(no_feasible_reason)+"</p>" if no_feasible_reason and not top_results else "",
        "evidence_status": _render_evidence_status(top_results),
        "how_to_read": _render_how_to_read(target),
        "builder_section": _render_builder(top_results, builder_files or {}, drivers),
        "rejection_section": rejection_section,
        "electrical_section": electrical_section,
        "assessment_status": "Experimental candidates, not validated recommendations" if top_results else "No feasible design in the evaluated set",
        "observer_summary": ("The observer integrates the nonuniform axisymmetric aperture velocity in an infinite baffle."
                             if top_results and all(r.get("output_metric")=="modal_baffled_on_axis" for r in top_results)
                             else "The observer estimate assumes a uniformly moving circular aperture in an infinite baffle."),
        "physics_summary": html.escape(", ".join(sorted({str(r.get("loss_model", "lossless"))+" / "+str(r.get("radiation_model", "legacy unknown")) for r in top_results})) or "No candidates"),
        "drive_voltage": target.voltage_rms,
        "stress_summary": _render_stress_summary(target, top_results),
        "observer_distance": target.observation_distance_m,
        "ripple_limit": target.max_ripple_db,
        "evidence_summary": _render_evidence_summary(top_results),
        "target_low": target.f_low_hz,
        "target_high": target.f_high_hz,
        "throat_display": throat_display,
        "mouth_radius": mouth_display,
        "horn_length": length_display,
        "profiles": ", ".join(sorted(solver_csvs.keys())),
        "timestamp": timestamp,
        "n_scored": n_scored,
        "n_top": len(top_results),
        "best_score": best_score,
        "best_bw": best_bw,
        "best_sens": best_sens,
        "best_ripple": best_ripple,
        "design_summary_section": design_summary_section,
        "geometry_header_cols": geometry_header_cols,
        "geometry_section": geometry_section,
        "rankings_rows": rankings_rows,
        "price_header": "<th>Price (as listed, shipping excluded)</th>" if pricing else "",
        "pricing_section": _render_pricing(pricing),
        "scoring_section": _render_scoring(scoring or describe_scoring(target, derived_geometry), all_ranked),
        "plot_coupled_spl": plot_coupled_spl,
        "plot_raw_spl": plot_raw_spl,
        "plot_impedance": plot_impedance,
        "plot_phase": plot_phase,
        "drivers_rows": drivers_rows,
        "directivity_section": _render_directivity(directivity_tables, target, top_results),
        "radiation_caveat": "<p>" + html.escape(FLANGED_PISTON_CAVEAT) + "</p>"
                            if any(r.get("radiation_model") == "flanged_piston" for r in top_results) else "",
    })
