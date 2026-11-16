#!/usr/bin/env python3
"""Compare horn_core.webster TMM input impedance with openwind (Inria).

openwind is GPL-3 and is never imported by the project packages. Run this
script in the separate environment from validation/cross_software/requirements.txt
with packages/horn-core/src on PYTHONPATH (``just validate-cross-software``).

Both codes are one-dimensional plane-wave models. openwind solves the
same physics with its own 1D finite elements on the true (conical-piece)
bore, so the comparison isolates our staircase TMM, our loss formulas and
our radiation load. Air properties are taken from openwind at 20 C and
passed to our code, so the gas model is not part of the difference.
"""
import argparse
import json
import math
import time
import warnings
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks

from horn_core.duct import AirProperties
from horn_core.profiles import PROFILE_NAMES, get_radius_func, natural_tractrix_length
from horn_core.webster import compute_horn_transfer_tmm, piston_radiation_impedance

ROOT = Path(__file__).resolve().parents[2]
THROAT, MOUTH, LENGTH = 0.025, 0.150, 0.300
TEMPERATURE_C = 20.0
F_MIN, F_MAX, F_STEP = 50.0, 5000.0, 2.0
PRODUCTION_SEGMENTS = 200   # lem_prescreen default
FINE_SEGMENTS = 2000
BORE_POINTS = 801           # cone pieces given to openwind for curved flares

# Frozen before the first comparison was run. Do not loosen to obtain a pass.
# "same_physics": identical propagation, loss model family and radiation load.
# "radiation_model": our exact Struve piston load against openwind's
#   causal fit to the flanged pipe (Silva et al. 2009), inside ka <= 3.
TOLERANCES = {
    "same_physics": {"db_median": 0.1, "db_p95": 0.5, "db_max": 2.0,
                     "phase_median_deg": 1.0, "phase_p95_deg": 5.0, "phase_max_deg": 20.0,
                     "peak_freq_rel": 0.005},
    "radiation_model": {"db_median": 0.25, "db_p95": 1.0, "db_max": 2.0,
                        "phase_median_deg": 2.0, "phase_p95_deg": 8.0, "phase_max_deg": 20.0,
                        "peak_freq_rel": 0.01, "ka_max": 3.0},
}

# (case id, tolerance class, our radiation, openwind radiation, our losses, openwind losses)
CASES = [
    ("matched_lossless", "same_physics", "plane_wave", "total_transmission", "lossless", False),
    ("closed_bessel", "same_physics", "closed", "closed", "kirchhoff", "bessel"),
    ("piston_lossless", "same_physics", "flanged_piston", "piston_data", "lossless", False),
    ("piston_bessel", "same_physics", "flanged_piston", "piston_data", "kirchhoff", "bessel"),
    ("piston_thin_layer", "same_physics", "flanged_piston", "piston_data", "boundary_layer", "wl"),
    ("flanged_fit_bessel", "radiation_model", "flanged_piston", "flanged_non_causal", "kirchhoff", "bessel"),
]


def horn_length(profile):
    if profile == "tractrix":
        return natural_tractrix_length(THROAT, MOUTH)
    return LENGTH


def openwind_bore(profile, length, force_points=False):
    """Main bore in openwind's list format, in metres and radii."""
    if profile == "conical":
        return [[0.0, length, THROAT, MOUTH, "linear"]]
    if profile == "exponential" and not force_points:
        return [[0.0, length, THROAT, MOUTH, "exponential"]]
    radius = get_radius_func(profile, THROAT, MOUTH, length)
    # Tractrix wall is vertical at the mouth, so cluster points there.
    s = np.linspace(0.0, 1.0, BORE_POINTS)
    z = length * (1 - (1 - s) ** 2) if profile == "tractrix" else length * s
    return [[float(zi), float(radius(zi))] for zi in z]


def openwind_air():
    from openwind.continuous.physics import Physics
    phys = Physics(TEMPERATURE_C)
    return AirProperties(c=float(phys.c(0)), rho=float(phys.rho(0)), gamma=float(phys.gamma(0)),
                         viscosity=float(phys.mu(0)), conductivity=float(phys.kappa(0)),
                         heat_capacity=float(phys.Cp(0)))


def run_openwind(freqs, bore, radiation, losses, air, method="FEM", spherical_waves=False):
    from openwind import ImpedanceComputation
    if radiation == "piston_data":
        k = 2 * np.pi * freqs / air.c
        data = (freqs, np.array([piston_radiation_impedance(ki, MOUTH) for ki in k]))
        radiation = ("from_data", (data, TEMPERATURE_C, MOUTH))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = ImpedanceComputation(freqs, bore, temperature=TEMPERATURE_C, losses=losses,
                                      radiation_category=radiation, compute_method=method,
                                      spherical_waves=spherical_waves)
    return np.asarray(result.impedance) / result.Zc


def run_tmm(freqs, profile, length, radiation, loss_model, air, segments):
    out = compute_horn_transfer_tmm(freqs, get_radius_func(profile, THROAT, MOUTH, length), length,
                                    THROAT, MOUTH, n_segments=segments, radiation_model=radiation,
                                    loss_model=loss_model, air=air)
    return (out["z_real"] + 1j * out["z_imag"]) / (air.rho * air.c)


def peaks(freqs, z):
    """Refined frequencies of |Z| maxima (parabolic fit in dB)."""
    level = 20 * np.log10(np.abs(z))
    idx, _ = find_peaks(level, prominence=1.0)
    out = []
    for i in idx:
        if 0 < i < len(level) - 1:
            a, b, c = level[i - 1], level[i], level[i + 1]
            den = a - 2 * b + c
            shift = 0.5 * (a - c) / den if den != 0 else 0.0
            out.append(float(freqs[i] + shift * (freqs[1] - freqs[0])))
    return out


def statistics(freqs, z, ref, tol, mask=None):
    mask = np.ones(len(freqs), bool) if mask is None else mask
    db = np.abs(20 * np.log10(np.abs(z[mask]) / np.abs(ref[mask])))
    ph = np.abs(np.degrees(np.angle(z[mask] / ref[mask])))
    ref_peaks, our_peaks = peaks(freqs[mask], ref[mask]), peaks(freqs[mask], z[mask])
    peak_rows = []
    for fp in ref_peaks:
        if our_peaks:
            nearest = min(our_peaks, key=lambda q: abs(q - fp))
            peak_rows.append({"reference_hz": fp, "ours_hz": nearest, "rel_error": (nearest - fp) / fp})
    stats = {
        "db_median": float(np.median(db)), "db_p95": float(np.percentile(db, 95)), "db_max": float(db.max()),
        "phase_median_deg": float(np.median(ph)), "phase_p95_deg": float(np.percentile(ph, 95)),
        "phase_max_deg": float(ph.max()),
        "peak_count_reference": len(ref_peaks), "peak_count_ours": len(our_peaks),
        "peak_freq_rel_max": float(max((abs(r["rel_error"]) for r in peak_rows), default=0.0)),
        "peaks": peak_rows,
        "f_range_hz": [float(freqs[mask][0]), float(freqs[mask][-1])],
    }
    checks = {key: stats[key] <= limit for key, limit in tol.items() if key in stats}
    checks["peak_freq_rel"] = stats["peak_freq_rel_max"] <= tol["peak_freq_rel"]
    checks["peak_count"] = len(ref_peaks) == len(our_peaks)
    stats["checks"] = checks
    stats["pass"] = all(checks.values())
    return stats


def plot(path, freqs, curves, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    for label, z, style in curves:
        a1.semilogx(freqs, 20 * np.log10(np.abs(z)), style, label=label, lw=1)
        a2.semilogx(freqs, np.degrees(np.angle(z)), style, lw=1)
    a1.set_ylabel("|Z/Zc| (dB)")
    a2.set_ylabel("phase (deg)")
    a2.set_xlabel("frequency (Hz)")
    a1.set_title(title)
    a1.legend(fontsize=8)
    a1.grid(True, which="both", alpha=.3)
    a2.grid(True, which="both", alpha=.3)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path, default=ROOT / "data/validation/cross_software_openwind.json")
    parser.add_argument("--plot-dir", type=Path, default=ROOT / "docs/validation/cross_software")
    parser.add_argument("--profiles", default=",".join(PROFILE_NAMES))
    parser.add_argument("--save-curves", type=Path, help="Also write the openwind curves (npz) for the FEM comparison")
    args = parser.parse_args()
    import openwind
    args.plot_dir.mkdir(parents=True, exist_ok=True)
    air = openwind_air()
    freqs = np.arange(F_MIN, F_MAX + F_STEP / 2, F_STEP)
    ka = 2 * np.pi * freqs / air.c * MOUTH
    results, curves_out = [], {}
    for profile in args.profiles.split(","):
        length = horn_length(profile)
        bore = openwind_bore(profile, length)
        entry = {"profile": profile, "throat_radius_m": THROAT, "mouth_radius_m": MOUTH, "length_m": length,
                 "openwind_bore": "exact segment" if profile in ("conical", "exponential")
                 else f"{BORE_POINTS} points, conical pieces", "cases": {}}
        for case, tol_class, our_rad, ow_rad, our_loss, ow_loss in CASES:
            tol = TOLERANCES[tol_class]
            t0 = time.time()
            ref = run_openwind(freqs, bore, ow_rad, ow_loss, air)
            t_ow = time.time() - t0
            z = run_tmm(freqs, profile, length, our_rad, our_loss, air, PRODUCTION_SEGMENTS)
            z_fine = run_tmm(freqs, profile, length, our_rad, our_loss, air, FINE_SEGMENTS)
            mask = ka <= tol["ka_max"] if "ka_max" in tol else None
            tol_checked = {k: v for k, v in tol.items() if k != "ka_max"}
            row = {"tolerance_class": tol_class, "ours": {"radiation": our_rad, "loss_model": our_loss},
                   "openwind": {"radiation": ow_rad, "losses": ow_loss, "method": "FEM"},
                   "openwind_seconds": round(t_ow, 2),
                   f"tmm_{PRODUCTION_SEGMENTS}": statistics(freqs, z, ref, tol_checked, mask),
                   f"tmm_{FINE_SEGMENTS}": statistics(freqs, z_fine, ref, tol_checked, mask)}
            if mask is not None:
                row["full_band_report_only"] = statistics(freqs, z_fine, ref, tol_checked)
                row["full_band_report_only"].pop("checks"); row["full_band_report_only"].pop("pass")
            if case == "matched_lossless":
                # openwind's own TMM on cone pieces checks the reference itself
                # (its TMM accepts cones only, so the exponential is sampled).
                tmm_bore = openwind_bore(profile, length, force_points=profile == "exponential")
                ref_tmm = run_openwind(freqs, tmm_bore, ow_rad, ow_loss, air, method="TMM")
                row["openwind_fem_vs_openwind_tmm"] = statistics(freqs, ref_tmm, ref, tol_checked)
            entry["cases"][case] = row
            curves_out[f"{profile}/{case}/openwind"] = ref
            curves_out[f"{profile}/{case}/tmm{PRODUCTION_SEGMENTS}"] = z
            if case in ("piston_bessel", "closed_bessel"):
                plot(args.plot_dir / f"openwind_{profile}_{case}.png", freqs,
                     [("openwind FEM", ref, "k-"), (f"TMM {PRODUCTION_SEGMENTS}", z, "r--"),
                      (f"TMM {FINE_SEGMENTS}", z_fine, "b:")],
                     f"{profile}: {case} (r_t {THROAT*1e3:.0f} mm, r_m {MOUTH*1e3:.0f} mm, L {length*1e3:.1f} mm)")
            s = row[f"tmm_{PRODUCTION_SEGMENTS}"]
            print(f"{profile:16s} {case:20s} {'PASS' if s['pass'] else 'FAIL'} "
                  f"dB med {s['db_median']:.3f} p95 {s['db_p95']:.3f} max {s['db_max']:.3f} "
                  f"ph p95 {s['phase_p95_deg']:.2f} peak {s['peak_freq_rel_max']:.4f} "
                  f"fine pass {row[f'tmm_{FINE_SEGMENTS}']['pass']}", flush=True)
        results.append(entry)
    summary = {
        "tool": "openwind", "tool_version": openwind.__version__, "tool_license": "GPL-3.0",
        "comparison": "plane-wave input impedance Z/Zc at the throat, 1D vs 1D",
        "frequencies_hz": {"min": F_MIN, "max": F_MAX, "step": F_STEP},
        "air_from_openwind_at_c": TEMPERATURE_C,
        "air": vars(air), "tolerances": TOLERANCES,
        "results": results,
        "all_pass_production_tmm": all(c[f"tmm_{PRODUCTION_SEGMENTS}"]["pass"]
                                       for r in results for c in r["cases"].values()),
        "all_pass_fine_tmm": all(c[f"tmm_{FINE_SEGMENTS}"]["pass"] for r in results for c in r["cases"].values()),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")
    if args.save_curves:
        np.savez_compressed(args.save_curves, freqs=freqs, **{k.replace("/", "__"): v for k, v in curves_out.items()})
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
