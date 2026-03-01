"""Self-contained HTML report for single-mode pipeline results.

Reads existing PNG files from disk and base64-encodes them into a single
HTML file with KPI cards and professional styling matching the auto-mode report.
"""

import argparse
import base64
import html
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


# -- Profile badge styles (copied from html_report.py; modules kept independent) --

_BADGE_STYLES = {
    "conical":     {"badge_bg": "#dbeafe", "badge_fg": "#1e40af"},
    "exponential": {"badge_bg": "#dcfce7", "badge_fg": "#166534"},
    "hyperbolic":  {"badge_bg": "#fef3c7", "badge_fg": "#92400e"},
}

_DEFAULT_BADGE = {"badge_bg": "#f3f4f6", "badge_fg": "#374151"}


def _badge_for(profile: str) -> dict:
    return _BADGE_STYLES.get(profile, _DEFAULT_BADGE)


def _profile_badge(profile: str) -> str:
    badge = _badge_for(profile)
    return (
        f'<span style="display:inline-block;padding:2px 8px;border-radius:4px;'
        f'font-size:0.8em;font-weight:600;'
        f'background:{badge["badge_bg"]};color:{badge["badge_fg"]}">'
        f'{html.escape(profile.capitalize())}</span>'
    )


def _fmt(value, fmt: str = ".1f", fallback: str = "\u2014") -> str:
    """Safely format an Optional numeric value."""
    if value is None:
        return fallback
    try:
        return f"{value:{fmt}}"
    except (ValueError, TypeError):
        return fallback


def _img_b64(path: str) -> str:
    """Read a PNG file and return a data-URI string."""
    data = Path(path).read_bytes()
    encoded = base64.b64encode(data).decode("ascii")
    return f"data:image/png;base64,{encoded}"


# -- HTML template -------------------------------------------------------------

_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Horn Simulation Report</title>
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

<h1>Horn Simulation Report</h1>
<p class="subtitle">
  {profile_badge} &nbsp;|&nbsp;
  Throat: {throat_mm:.1f} mm &nbsp;|&nbsp;
  Mouth: {mouth_mm:.1f} mm &nbsp;|&nbsp;
  Length: {length_mm:.1f} mm &nbsp;|&nbsp;
  {freq_range}
  Generated: {timestamp}
</p>

<div class="cards">
  <div class="card"><div class="label">Peak SPL (dB)</div><div class="value">{peak_spl}</div></div>
  <div class="card"><div class="label">Peak Freq (Hz)</div><div class="value">{peak_freq}</div></div>
  <div class="card"><div class="label">f3 Low (Hz)</div><div class="value">{f3_low}</div></div>
  <div class="card"><div class="label">f3 High (Hz)</div><div class="value">{f3_high}</div></div>
  <div class="card"><div class="label">Bandwidth (Hz)</div><div class="value">{bandwidth_hz}</div></div>
  <div class="card"><div class="label">Bandwidth (oct)</div><div class="value">{bandwidth_oct}</div></div>
  <div class="card"><div class="label">Ripple (dB)</div><div class="value">{ripple}</div></div>
  <div class="card"><div class="label">Avg Sensitivity (dB)</div><div class="value">{avg_sens}</div></div>
</div>

<h2>Horn Geometry</h2>
<div class="plot-grid">
  <div class="plot"><img src="{horn_3d_src}" alt="3D horn render"></div>
  <div class="design-summary">
    <dl>
      <dt>Profile</dt><dd>{profile_badge_dd}</dd>
      <dt>Throat radius</dt><dd>{throat_radius_m:.4f} m</dd>
      <dt>Mouth radius</dt><dd>{mouth_radius_m:.4f} m</dd>
      <dt>Horn length</dt><dd>{length_m:.3f} m</dd>
      <dt>Throat diameter</dt><dd>{throat_mm:.1f} mm</dd>
      <dt>Mouth diameter</dt><dd>{mouth_mm:.1f} mm</dd>
    </dl>
  </div>
</div>

<h2>Frequency Response</h2>
<div class="plot"><img src="{spl_src}" alt="Frequency response"></div>

<h2>Acoustic Dashboard</h2>
<div class="plot"><img src="{dashboard_src}" alt="Dashboard"></div>

<h2>Impedance &amp; Phase</h2>
<div class="plot-grid">
  <div class="plot"><img src="{impedance_src}" alt="Impedance"></div>
  <div class="plot"><img src="{phase_src}" alt="Phase response"></div>
</div>

{directivity_section}

<div class="footer">Horn Simulation Report &mdash; generated {timestamp}</div>

</div>
</body>
</html>
"""


def generate_single_report(
    throat_radius: float,
    mouth_radius: float,
    length: float,
    profile: str,
    kpis: dict,
    spl_png: str,
    impedance_png: str,
    phase_png: str,
    dashboard_png: str,
    horn_3d_png: str,
    final_csv: Optional[str] = None,
    directivity_polar_png: Optional[str] = None,
    directivity_contour_png: Optional[str] = None,
    beamwidth_png: Optional[str] = None,
    directivity_index_png: Optional[str] = None,
) -> str:
    """Generate a self-contained HTML report string for a single-mode run.

    All PNG paths are read from disk and base64-encoded into the HTML.
    """
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    # Freq range from CSV if available
    freq_range = ""
    if final_csv:
        try:
            import pandas as pd
            df = pd.read_csv(final_csv)
            fmin = df["frequency"].min()
            fmax = df["frequency"].max()
            freq_range = f"{fmin:.0f} \u2014 {fmax:.0f} Hz &nbsp;|&nbsp; "
        except Exception:
            pass

    # Directivity section
    directivity_section = ""
    has_dir = all(p is not None for p in [
        directivity_polar_png, directivity_contour_png,
        beamwidth_png, directivity_index_png,
    ])
    if has_dir:
        directivity_section = (
            '<h2>Directivity Analysis</h2>\n'
            '<div class="plot-grid">\n'
            f'  <div class="plot"><img src="{_img_b64(directivity_polar_png)}" alt="Polar directivity"></div>\n'
            f'  <div class="plot"><img src="{_img_b64(directivity_contour_png)}" alt="Directivity contour"></div>\n'
            '</div>\n'
            '<div class="plot-grid">\n'
            f'  <div class="plot"><img src="{_img_b64(beamwidth_png)}" alt="Beamwidth"></div>\n'
            f'  <div class="plot"><img src="{_img_b64(directivity_index_png)}" alt="Directivity index"></div>\n'
            '</div>'
        )

    throat_mm = throat_radius * 2 * 1000
    mouth_mm = mouth_radius * 2 * 1000
    length_mm = length * 1000

    return _HTML_TEMPLATE.format_map({
        "profile_badge": _profile_badge(profile),
        "profile_badge_dd": _profile_badge(profile),
        "throat_radius_m": throat_radius,
        "mouth_radius_m": mouth_radius,
        "length_m": length,
        "throat_mm": throat_mm,
        "mouth_mm": mouth_mm,
        "length_mm": length_mm,
        "freq_range": freq_range,
        "timestamp": timestamp,
        # KPI cards
        "peak_spl": _fmt(kpis.get("peak_spl_db"), ".1f"),
        "peak_freq": _fmt(kpis.get("peak_freq_hz"), ".0f"),
        "f3_low": _fmt(kpis.get("f3_low_hz"), ".0f"),
        "f3_high": _fmt(kpis.get("f3_high_hz"), ".0f"),
        "bandwidth_hz": _fmt(kpis.get("bandwidth_hz"), ".0f"),
        "bandwidth_oct": _fmt(kpis.get("bandwidth_octaves"), ".2f"),
        "ripple": _fmt(kpis.get("passband_ripple_db"), ".1f"),
        "avg_sens": _fmt(kpis.get("avg_sensitivity_db"), ".1f"),
        # Embedded images
        "horn_3d_src": _img_b64(horn_3d_png),
        "spl_src": _img_b64(spl_png),
        "impedance_src": _img_b64(impedance_png),
        "phase_src": _img_b64(phase_png),
        "dashboard_src": _img_b64(dashboard_png),
        # Conditional section
        "directivity_section": directivity_section,
    })


def main():
    parser = argparse.ArgumentParser(description="Generate single-mode HTML report")
    parser.add_argument("--kpis", required=True, help="Path to kpis.json")
    parser.add_argument("--final-csv", default=None, help="Path to final_results.csv")
    parser.add_argument("--throat-radius", type=float, required=True)
    parser.add_argument("--mouth-radius", type=float, required=True)
    parser.add_argument("--length", type=float, required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--spl-png", required=True)
    parser.add_argument("--impedance-png", required=True)
    parser.add_argument("--phase-png", required=True)
    parser.add_argument("--dashboard-png", required=True)
    parser.add_argument("--horn-3d-png", required=True)
    parser.add_argument("--polar-png", default=None)
    parser.add_argument("--contour-png", default=None)
    parser.add_argument("--beamwidth-png", default=None)
    parser.add_argument("--di-png", default=None)
    parser.add_argument("--output", required=True, help="Output HTML path")
    args = parser.parse_args()

    kpis = json.loads(Path(args.kpis).read_text())

    report_html = generate_single_report(
        throat_radius=args.throat_radius,
        mouth_radius=args.mouth_radius,
        length=args.length,
        profile=args.profile,
        kpis=kpis,
        spl_png=args.spl_png,
        impedance_png=args.impedance_png,
        phase_png=args.phase_png,
        dashboard_png=args.dashboard_png,
        horn_3d_png=args.horn_3d_png,
        final_csv=args.final_csv,
        directivity_polar_png=args.polar_png,
        directivity_contour_png=args.contour_png,
        beamwidth_png=args.beamwidth_png,
        directivity_index_png=args.di_png,
    )

    out = Path(args.output)
    out.write_text(report_html)
    print(f"Report written to {out}")


if __name__ == "__main__":
    main()
