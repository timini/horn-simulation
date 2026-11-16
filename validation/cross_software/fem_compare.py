#!/usr/bin/env python3
"""Compare horn-solver FEM throat impedance with openwind and our TMM.

Each argument is a run directory from scripts/run_pipeline.py (single mode,
--radiation_model flanged_piston --loss_model lossless). The FEM, openwind and
the TMM then share the same baffled-piston mouth load, so differences come
from the 3D interior field: curved wavefronts and, at high ka, higher modes.

Lossless Z/Zc depends on frequency only through k, so openwind (which fixes
its own air from temperature) is evaluated at f * c_openwind / c_fem.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

import openwind_compare as ow
from horn_core.duct import AirProperties
from horn_core.profiles import get_radius_func
from horn_core.webster import compute_horn_transfer_tmm

ROOT = Path(__file__).resolve().parents[2]

# Frozen before the first FEM run finished. FEM is 3D, the others are 1D.
TOLERANCES = {
    "low_ka": {"ka_max": 1.5, "db_median": 0.2, "db_p95": 0.5, "db_max": 1.0, "phase_p95_deg": 3.0},
    "high_ka": "report only: plane-wave 1D models omit higher-order modes and wavefront curvature",
}


def load_fem(run_dir):
    spec = json.loads((run_dir / "outputs/resolved_specification.json").read_text())["parameters"]
    path = run_dir / "outputs/final_results.csv"
    return path, spec, list(csv.DictReader(open(path)))


def stats(z, ref, mask):
    db = np.abs(20 * np.log10(np.abs(z[mask]) / np.abs(ref[mask])))
    ph = np.abs(np.degrees(np.angle(z[mask] / ref[mask])))
    return {"points": int(mask.sum()), "db_median": float(np.median(db)), "db_p95": float(np.percentile(db, 95)),
            "db_max": float(db.max()), "phase_median_deg": float(np.median(ph)),
            "phase_p95_deg": float(np.percentile(ph, 95)), "phase_max_deg": float(ph.max())}


def compare(run_dir, plot_dir):
    path, p, rows = load_fem(run_dir)
    profile, r_t, r_m, L = p["profile"], float(p["throat_radius"]), float(p["mouth_radius"]), float(p["length"])
    f = np.array([float(r["frequency"]) for r in rows])
    rho, c = float(rows[0].get("air_rho_kg_m3", 1.225)), float(rows[0].get("air_c_m_s", 343.0))
    if rows[0].get("radiation_model") != "flanged_piston" or rows[0].get("loss_model") != "lossless":
        raise ValueError("FEM comparison needs a lossless flanged_piston run")
    z_fem = np.array([float(r["z_real"]) + 1j * float(r["z_imag"]) for r in rows]) / (rho * c)
    air = AirProperties(c=c, rho=rho)
    tmm = compute_horn_transfer_tmm(f, get_radius_func(profile, r_t, r_m, L), L, r_t, r_m, n_segments=2000,
                                    radiation_model="flanged_piston", air=air)
    z_tmm = (tmm["z_real"] + 1j * tmm["z_imag"]) / (rho * c)
    ow_air = ow.openwind_air()
    old = (ow.THROAT, ow.MOUTH)
    ow.THROAT, ow.MOUTH = r_t, r_m
    try:
        z_ow = ow.run_openwind(f * ow_air.c / c, ow.openwind_bore(profile, L), "piston_data", False, ow_air)
        # openwind's spherical-wave option uses curved wavefront lengths in
        # flaring pipes; it shows how much of the FEM difference is the
        # plane-wavefront assumption shared by both 1D models.
        z_ow_sph = ow.run_openwind(f * ow_air.c / c, ow.openwind_bore(profile, L), "piston_data", False, ow_air,
                                   spherical_waves=True)
    finally:
        ow.THROAT, ow.MOUTH = old
    ka = 2 * np.pi * f / c * r_m
    tol = TOLERANCES["low_ka"]
    low, high = ka <= tol["ka_max"], ka > tol["ka_max"]
    out = {"run": str(run_dir.relative_to(ROOT)) if run_dir.is_relative_to(ROOT) else str(run_dir),
           "fem_csv": str(path.relative_to(run_dir)), "profile": profile,
           "geometry": {"throat_radius_m": r_t, "mouth_radius_m": r_m, "length_m": L},
           "mesh_cells_max": max(int(float(r.get("mesh_cells", 0))) for r in rows),
           "frequency_hz": [float(f[0]), float(f[-1])], "points": len(f), "air": {"rho": rho, "c": c}}
    for name, ref in (("fem_vs_openwind", z_ow), ("fem_vs_tmm", z_tmm), ("tmm_vs_openwind", None),
                      ("fem_vs_openwind_spherical_waves", z_ow_sph)):
        a, b = (z_tmm, z_ow) if ref is None else (z_fem, ref)
        lo = stats(a, b, low) if low.any() else None
        entry = {"low_ka": lo, "high_ka_report_only": stats(a, b, high) if high.any() else None}
        if lo is not None:
            entry["low_ka_checks"] = {k: lo[k] <= v for k, v in tol.items() if k in lo}
            entry["low_ka_pass"] = all(entry["low_ka_checks"].values())
        out[name] = entry
    # dB difference against mouth ka, for the doc.
    edges = [0, 0.5, 1, 1.5, 2, 3, 5, 8, 15]
    for key, ref in (("fem_vs_openwind_by_ka", z_ow), ("fem_vs_openwind_spherical_by_ka", z_ow_sph)):
        out[key] = [{"ka": [a, b], **stats(z_fem, ref, (ka > a) & (ka <= b))}
                    for a, b in zip(edges[:-1], edges[1:]) if ((ka > a) & (ka <= b)).any()]
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    for label, z, style in (("horn-solver FEM (3D)", z_fem, "k-"), ("openwind (1D FEM)", z_ow, "r--"),
                            ("openwind, spherical waves", z_ow_sph, "g-."),
                            ("TMM 2000", z_tmm, "b:")):
        a1.semilogx(f, 20 * np.log10(np.abs(z)), style, label=label, lw=1)
        a2.semilogx(f, np.degrees(np.angle(z)), style, lw=1)
    f15 = 1.5 * c / (2 * np.pi * r_m)
    for ax in (a1, a2):
        ax.axvline(f15, color="0.6", lw=.8)
        ax.grid(True, which="both", alpha=.3)
    a1.set_title(f"{profile}: throat impedance, baffled piston mouth (grey line: ka = 1.5)")
    a1.set_ylabel("|Z/Zc| (dB)")
    a2.set_ylabel("phase (deg)")
    a2.set_xlabel("frequency (Hz)")
    a1.legend(fontsize=8)
    fig.tight_layout()
    plot_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_dir / f"fem_{profile}.png", dpi=110)
    plt.close(fig)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("runs", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, default=ROOT / "data/validation/cross_software_fem.json")
    parser.add_argument("--plot-dir", type=Path, default=ROOT / "docs/validation/cross_software")
    args = parser.parse_args()
    results = [compare(r.resolve(), args.plot_dir) for r in args.runs]
    for r in results:
        s = r["fem_vs_openwind"]
        print(r["profile"], "low ka", json.dumps(s["low_ka"]), "pass", s.get("low_ka_pass"))
        print(r["profile"], "high ka", json.dumps(s["high_ka_report_only"]))
    args.output.write_text(json.dumps({"tolerances": TOLERANCES, "results": results}, indent=1) + "\n")


if __name__ == "__main__":
    main()
