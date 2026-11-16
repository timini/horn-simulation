#!/usr/bin/env python3
"""Compare horn-solver FEM throat impedance with openwind and our TMM.

Each argument is a run directory from scripts/run_pipeline.py (single mode,
--radiation_model flanged_piston --loss_model lossless). All three codes then
use the exact baffled-piston impedance at the mouth, but not in the same way:
the FEM applies it pointwise as a local Robin condition on the flat mouth disc,
the 1D codes as a lumped load on the mean mouth velocity. The FEM also drives
a flat throat disc at uniform pressure and has a curved interior field.

--refinement COARSE FINE compares two FEM runs of one horn at their common
frequencies, to test whether the FEM mesh explains the difference.

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

# Written before this script was first run; committed (a7de55e, 16:08 UTC)
# after the conical FEM run had finished (16:01 UTC) and before the
# exponential run finished. FEM is 3D, the others are 1D.
TOLERANCES = {
    "low_ka": {"ka_max": 1.5, "db_median": 0.2, "db_p95": 0.5, "db_max": 1.0, "phase_p95_deg": 3.0},
    "high_ka": "report only: plane-wave 1D models omit higher-order modes and wavefront curvature",
    # Review fix pass, committed before the refined run started: the mesh is
    # ruled out as the cause if refining changes |Z| by at most 0.1 dB and
    # 1 degree for ka <= 1.5 (the BEM refinement limit), far below the
    # 0.4 to 1.6 dB FEM against 1D differences.
    "mesh_refinement": {"ka_max": 1.5, "db_max": 0.1, "phase_max_deg": 1.0},
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


def provenance(run_dir):
    """Source revision, container image ids and mesh settings of a run."""
    path = run_dir / "manifest.json"
    if not path.exists():
        return None
    m = json.loads(path.read_text())
    spec = json.loads((run_dir / "outputs/resolved_specification.json").read_text())["parameters"]
    return {"source_revision": m.get("source_revision"), "started_at": m.get("started_at"),
            "finished_at": m.get("finished_at"), "arguments": m.get("arguments"),
            "container_image_ids": {k: v.get("id") for k, v in (m.get("containers") or {}).items()},
            "mesh_size_m": spec.get("mesh_size"), "element_degree": spec.get("element_degree"),
            "max_mesh_nodes": spec.get("max_mesh_nodes")}


def fem_impedance(run_dir):
    _, _, rows = load_fem(run_dir)
    rho, c = float(rows[0].get("air_rho_kg_m3", 1.225)), float(rows[0].get("air_c_m_s", 343.0))
    f = np.array([float(r["frequency"]) for r in rows])
    z = np.array([float(r["z_real"]) + 1j * float(r["z_imag"]) for r in rows]) / (rho * c)
    cells = np.array([int(float(r.get("mesh_cells", 0))) for r in rows])
    return f, z, cells, c


def refinement(coarse_dir, fine_dir):
    """FEM against FEM on a finer mesh, at the frequencies both runs share."""
    fc, zc, cc, c = fem_impedance(coarse_dir)
    ff, zf, cf, _ = fem_impedance(fine_dir)
    pc = json.loads((coarse_dir / "outputs/resolved_specification.json").read_text())["parameters"]
    pf = json.loads((fine_dir / "outputs/resolved_specification.json").read_text())["parameters"]
    for key in ("profile", "throat_radius", "mouth_radius", "length", "radiation_model", "loss_model"):
        if str(pc.get(key)) != str(pf.get(key)):
            raise ValueError(f"refinement runs differ in {key}")
    ic, jf = [], []
    for i, f in enumerate(fc):
        j = int(np.argmin(np.abs(ff - f)))
        if abs(ff[j] - f) <= 1e-6 * f:
            ic.append(i)
            jf.append(j)
    ic, jf = np.array(ic), np.array(jf)
    f = fc[ic]
    ka = 2 * np.pi * f / c * float(pc["mouth_radius"])
    tol = TOLERANCES["mesh_refinement"]
    low = ka <= tol["ka_max"]
    lo = stats(zf[jf], zc[ic], low)
    edges = [0, 0.5, 1, 1.5, 2, 3]
    return {"profile": pc["profile"], "coarse_run": _rel(coarse_dir), "fine_run": _rel(fine_dir),
            "coarse": provenance(coarse_dir), "fine": provenance(fine_dir),
            "coarse_mesh_cells_in_shared_band": [int(cc[ic].min()), int(cc[ic].max())],
            "fine_mesh_cells_in_shared_band": [int(cf[jf].min()), int(cf[jf].max())],
            "shared_points": len(f), "frequency_hz": [float(f[0]), float(f[-1])],
            "ka_max_shared": float(ka.max()),
            "fine_vs_coarse_low_ka": lo,
            "checks": {"db_max": lo["db_max"] <= tol["db_max"], "phase_max_deg": lo["phase_max_deg"] <= tol["phase_max_deg"]},
            "fine_vs_coarse_by_ka": [{"ka": [a, b], **stats(zf[jf], zc[ic], (ka > a) & (ka <= b))}
                                     for a, b in zip(edges[:-1], edges[1:]) if ((ka > a) & (ka <= b)).any()]}


def _rel(run_dir):
    return str(run_dir.relative_to(ROOT)) if run_dir.is_relative_to(ROOT) else str(run_dir)


def compare(run_dir, plot_dir, plot=True):
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
    out = {"run": _rel(run_dir), "provenance": provenance(run_dir),
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
    if not plot:
        return out
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
    parser.add_argument("--refinement", nargs=2, type=Path, action="append", default=[],
                        metavar=("COARSE", "FINE"), help="Two runs of one horn, the second on a finer mesh")
    parser.add_argument("--output", type=Path, default=ROOT / "data/validation/cross_software_fem.json")
    parser.add_argument("--plot-dir", type=Path, default=ROOT / "docs/validation/cross_software")
    args = parser.parse_args()
    results = [compare(r.resolve(), args.plot_dir) for r in args.runs]
    for r in results:
        s = r["fem_vs_openwind"]
        print(r["profile"], "low ka", json.dumps(s["low_ka"]), "pass", s.get("low_ka_pass"))
        print(r["profile"], "high ka", json.dumps(s["high_ka_report_only"]))
    refined = []
    for coarse, fine in args.refinement:
        r = refinement(coarse.resolve(), fine.resolve())
        r["fine_vs_1d"] = compare(fine.resolve(), args.plot_dir, plot=False)
        refined.append(r)
        print(r["profile"], "refinement", json.dumps(r["fine_vs_coarse_low_ka"]), json.dumps(r["checks"]))
    args.output.write_text(json.dumps({"tolerances": TOLERANCES, "results": results, "mesh_refinement": refined},
                                      indent=1) + "\n")


if __name__ == "__main__":
    main()
