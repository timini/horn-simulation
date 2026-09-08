#!/usr/bin/env python3
"""Reproduce the full pipe reference audit, including production FEM convergence.

Run in the solver container from the repository root. The output directory must
be new: predictions and the protocol are frozen before reading reference values.
This is an unfitted nominal-geometry impedance audit, not assembly certification.
"""
import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import multiprocessing
from pathlib import Path
import subprocess
import time

import numpy as np
import pandas as pd

from validate_duct_physics import AIR, CASES, metrics
from horn_core.webster import compute_horn_transfer_tmm


UNIQUE_CASES = ("Brass_O", "Wood_O", "Brass_C", "Cone_O", "Cone_C")
MODEL_CASE = {**{c: c for c in CASES}, "3D_O": "Wood_O",
              "Wood_C": "Brass_C", "3D_C": "Brass_C"}
ALL_CASES = {**CASES, "Cylinder_U": (.007, .007, 0., "finite_flange"),
             "Cone_U": (.005, .0113, 0., "finite_flange")}
NUMERICAL_CASE = {
    "Cylinder_closed": "Brass_C", "Cone_closed": "Cone_C",
    "Cylinder_finite_flanged_width2mm": "Brass_O",
    "Cylinder_finite_flanged_width7mm": "Wood_O",
    "Cylinder_unflanged": "Cylinder_U", "Cone_unflanged": "Cone_U",
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def prediction(frequency, case, segments=400):
    throat, mouth, width, radiation = ALL_CASES[case]
    response = compute_horn_transfer_tmm(
        frequency, lambda z: throat + (mouth-throat)*z/.18, .18, throat, mouth,
        n_segments=segments, air=AIR, loss_model="boundary_layer",
        radiation_model=radiation, flange_width=width)
    return (response["z_real"] + 1j*response["z_imag"])/(AIR.rho*AIR.c)


def compare(predicted, reference):
    predicted, reference = np.asarray(predicted), np.asarray(reference)
    if predicted.shape != reference.shape or not predicted.size:
        raise ValueError("Comparison needs nonempty, matching shapes")
    if not np.isfinite(predicted).all() or not np.isfinite(reference).all():
        raise ValueError("Nonfinite impedance")
    result = metrics(predicted, reference)
    phase = np.abs(np.angle(predicted*np.conj(reference), deg=True))
    result.update(median_phase_error_deg=float(np.median(phase)),
                  p95_phase_error_deg=float(np.percentile(phase, 95)))
    return result


def interpolate(frame, frequency, real, imag):
    f = frame.frequency.to_numpy()
    values = frame[[real, imag]].to_numpy()
    if (len(f) < 2 or not np.isfinite(f).all() or not np.isfinite(values).all()
            or np.any(np.diff(f) <= 0) or f[0] > frequency[0] or f[-1] < frequency[-1]):
        raise ValueError("Reference must be finite, strictly ordered and cover the entire fixed band")
    return np.interp(frequency, f, values[:, 0]) + 1j*np.interp(frequency, f, values[:, 1])


def fem_impedance(frame):
    return (frame.z_real.to_numpy()+1j*frame.z_imag.to_numpy())/(AIR.rho*AIR.c)


def fem_health(frame):
    supplied = frame.input_acoustic_power_w.to_numpy()
    outgoing = (frame.mouth_acoustic_power_w + frame.viscous_wall_power_w
                + frame.thermal_wall_power_w).to_numpy()
    balance = np.abs(supplied-outgoing)/np.maximum(np.abs(supplied), 1e-30)
    finite = bool(np.isfinite(frame.select_dtypes(include="number")).all().all())
    result = {"finite": finite, "frequencies": len(frame),
              "mesh_cells": int(frame.mesh_cells.iloc[0]),
              "max_relative_residual": float(frame.relative_residual.max()),
              "min_converged_reason": int(frame.converged_reason.min()),
              "max_relative_power_balance_error": float(balance.max()),
              "minimum_real_normalized_impedance": float(fem_impedance(frame).real.min()),
              "nonnegative_power": bool((supplied >= 0).all() and (outgoing >= 0).all()
                  and (frame.mouth_acoustic_power_w >= -1e-15).all()
                  and (frame.viscous_wall_power_w >= -1e-15).all()
                  and (frame.thermal_wall_power_w >= -1e-15).all())}
    result["passed"] = bool(finite and result["max_relative_residual"] <= 1e-8
        and result["min_converged_reason"] > 0 and balance.max() <= 1e-7
        and result["minimum_real_normalized_impedance"] >= 0 and result["nonnegative_power"])
    return result


def run_fem(out, case, h, count):
    import gmsh
    from horn_solver.solver import run_simulation_from_step
    throat, mouth, width, radiation = CASES[case]
    step = out/f"{case}.step"
    if not step.exists():
        gmsh.initialize()
        try:
            gmsh.model.add(case)
            if throat == mouth:
                gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, .18, throat)
            else:
                gmsh.model.occ.addCone(0, 0, 0, 0, 0, .18, throat, mouth)
            gmsh.model.occ.synchronize()
            gmsh.write(str(step))
        finally:
            gmsh.finalize()
    csv = out/f"{case}-h{h:g}-n{count}.csv"
    print(f"AUDIT: production FEM {case}, mesh {h:g} m, {count} frequencies", flush=True)
    run_simulation_from_step(str(step), (110, 3900), count, {"length": .18}, str(csv),
        3900, mesh_size=h, element_degree=2, loss_model="boundary_layer",
        minimum_wall_scale=throat, radiation_model=radiation, flange_width=width, air=AIR)
    return pd.read_csv(csv)


def run_fem_case(task):
    out, case = task
    return case, tuple(run_fem(out, case, h, count)
                       for h, count in ((.004, 481), (.006, 121), (.003, 121)))


def aggregate(rows, key):
    result = {}
    for value in sorted({r[key] for r in rows}):
        group = [r for r in rows if r[key] == value]
        result[value] = {"curves": len(group)}
        for model in ("tmm_dense", "fem_481"):
            result[value][model] = {
                "magnitude_passes": sum(r[model]["magnitude_agreement"] for r in group),
                "median_p95_db": float(np.median([r[model]["p95_db"] for r in group])),
                "median_p95_phase_error_deg": float(np.median([r[model]["p95_phase_error_deg"] for r in group])),
                "median_p95_scaled_complex_error": float(np.median([r[model]["p95_scaled_complex_error"] for r in group])),
            }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1,
                        help="Independent solver processes (1–5); no shared meshes or state")
    args = parser.parse_args()
    out, root = args.output_dir, args.reference_dir
    if out.exists():
        parser.error("Use a new output directory to preserve previous evidence")
    if not 1 <= args.workers <= 5:
        parser.error("workers must be between 1 and 5")
    manifest = json.loads((root/"manifest.json").read_text())
    if sha(args.archive) != manifest["reference"]["sha256"]:
        raise ValueError("Source archive checksum mismatch")
    out.mkdir(parents=True)
    started = time.perf_counter()
    sources = sorted(set(Path("packages").glob("*/src/**/*.py")))
    sources += [Path(__file__), Path(__file__).with_name("validate_duct_physics.py")]
    protocol = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source": manifest["reference"], "manifest_sha256": sha(root/"manifest.json"),
        "air": asdict(AIR), "cases": ALL_CASES, "length_m": .18,
        "frequency_band_hz": [110, 3900], "tmm_samples": 20001,
        "fem_samples": 481, "convergence_samples": 121,
        "mesh_sizes_m": [.006, .004, .003], "element_degree": 2,
        "solver_workers": args.workers,
        "tmm_segments": [200, 400, 800], "fitted_parameters": [],
        "magnitude_limits_db": {"median": 2, "p95": 4},
        "mesh_convergence_limits": {"max_impedance_change_db": .5,
                                     "p95_scaled_complex_change": .05},
        "phase_gate": None, "independent_simulation_scope": "All 40 curves; differing loss and wavefront models are diagnostic comparisons, not identical-model certification.",
        "held_out_scope": "Historical split retained; these data were already examined in the previous audit and are not newly unseen validation data.",
        "model_aliases": MODEL_CASE,
        "model_alias_reason": "Identical nominal air-domain geometry and rigid-wall boundary conditions; material compliance is not modeled.",
        "source_code_sha256": {str(p): sha(p) for p in sources},
    }
    try:
        protocol["git_revision"] = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        protocol["git_revision"] = None
    write_json(out/"protocol.json", protocol)
    frequency = np.geomspace(110, 3900, 20001)
    predicted, sampled = {}, {}
    segmentation = {}
    for case in ALL_CASES:
        predicted[case] = prediction(frequency, case)
        z = predicted[case]
        pd.DataFrame({"frequency": frequency, "z_real_normalized": z.real,
                      "z_imag_normalized": z.imag}).to_csv(out/f"{case}-prediction.csv", index=False)
        if case in CASES:
            for n in (121, 241, 481, 961):
                f = np.geomspace(110, 3900, n)
                sampled[case, n] = prediction(f, case)
                s = sampled[case, n]
                pd.DataFrame({"frequency": f, "z_real_normalized": s.real,
                              "z_imag_normalized": s.imag}).to_csv(out/f"{case}-tmm-n{n}.csv", index=False)
        if case in UNIQUE_CASES:
            segmentation[case] = {"200_vs_400": compare(prediction(frequency, case, 200), z),
                                  "400_vs_800": compare(z, prediction(frequency, case, 800))}
    fem, convergence, health = {}, {}, {}
    # Spawn, rather than fork, so each dolfinx/MPI/gmsh runtime is independent.
    with multiprocessing.get_context("spawn").Pool(args.workers) as pool:
        responses = pool.map(run_fem_case, [(out, case) for case in UNIQUE_CASES])
    for case, (middle, coarse, fine) in responses:
        fem[case] = middle
        for name, frame in (("middle", middle), ("coarse", coarse), ("fine", fine)):
            health[f"{case}-{name}"] = fem_health(frame)
        subset = middle.iloc[::4]
        np.testing.assert_allclose(subset.frequency, fine.frequency, rtol=1e-12)
        change = compare(fem_impedance(subset), fem_impedance(fine))
        convergence[case] = {"coarse_vs_fine": compare(fem_impedance(coarse), fem_impedance(fine)),
            "middle_vs_fine": change,
            "passed": change["max_db"] <= .5 and change["p95_scaled_complex_error"] <= .05,
            "fem_vs_tmm_481": compare(fem_impedance(middle), prediction(middle.frequency.to_numpy(), case))}
    # Freeze every prediction before opening any of the reference curves.
    write_json(out/"prediction_hashes.json", {p.name: sha(p) for p in sorted(out.glob("*.csv"))})
    rows, numerical = [], []
    for curve in manifest["curves"]:
        if curve.get("duplicate_of"):
            raise ValueError("Unexpected duplicate: explicitly account for it before changing coverage")
        path = root/curve["csv"]
        if sha(path) != curve["csv_sha256"]:
            raise ValueError(f"Changed reference: {path}")
        frame = pd.read_csv(path)
        row = {k: curve[k] for k in ("source_member", "csv_sha256", "configuration", "operator")}
        if curve["reference_kind"] == "measurement":
            case = curve["configuration"]
            measured = interpolate(frame, frequency, "z_real_normalized", "z_imag_normalized")
            ff = fem[MODEL_CASE[case]]
            row["tmm_dense"] = compare(predicted[case], measured)
            row["fem_481"] = compare(fem_impedance(ff), interpolate(frame, ff.frequency.to_numpy(),
                "z_real_normalized", "z_imag_normalized"))
            row["historical_role"] = "development" if case == "Brass_O" else "held_out"
            row["sampling"] = {}
            for n in (121, 241, 481, 961):
                f = np.geomspace(110, 3900, n)
                row["sampling"][str(n)] = compare(sampled[case, n],
                    interpolate(frame, f, "z_real_normalized", "z_imag_normalized"))
            rows.append(row)
        elif curve["reference_kind"] == "independent_simulation":
            case = NUMERICAL_CASE[curve["configuration"]]
            z = interpolate(frame, frequency, "z_real_pa_s_per_m3", "z_imag_pa_s_per_m3")
            z *= np.pi*ALL_CASES[case][0]**2/(AIR.rho*AIR.c)
            row["tmm_dense"] = compare(predicted[case], z)
            numerical.append(row)
        else:
            raise ValueError(f"Unexpected reference kind: {curve['reference_kind']}")
    if len(rows) != 299 or len(numerical) != 40:
        raise ValueError("Unexpected reference coverage; audit requires 299 measured and 40 numerical curves")
    resolution = {str(n): {
        "magnitude_passes": sum(r["sampling"][str(n)]["magnitude_agreement"] for r in rows),
        "classification_changes_vs_dense": sum(r["sampling"][str(n)]["magnitude_agreement"] != r["tmm_dense"]["magnitude_agreement"] for r in rows),
        "largest_p95_error_change_db": max(abs(r["sampling"][str(n)]["p95_db"]-r["tmm_dense"]["p95_db"]) for r in rows),
    } for n in (121, 241, 481, 961)}
    result = {"protocol_sha256": sha(out/"protocol.json"),
        "prediction_hashes_sha256": sha(out/"prediction_hashes.json"),
        "source_archive_sha256": sha(args.archive), "measurements": rows,
        "independent_simulations": numerical, "by_case": aggregate(rows, "configuration"),
        "by_operator": aggregate(rows, "operator"), "by_historical_role": aggregate(rows, "historical_role"),
        "fem_health": health, "mesh_convergence": convergence,
        "tmm_segmentation": segmentation, "frequency_sampling": resolution,
        "totals": {"measured_curves": len(rows), "independent_simulation_curves": len(numerical),
            "tmm_magnitude_passes": sum(r["tmm_dense"]["magnitude_agreement"] for r in rows),
            "fem_sampled_magnitude_passes": sum(r["fem_481"]["magnitude_agreement"] for r in rows),
            "independent_simulation_magnitude_matches": sum(r["tmm_dense"]["magnitude_agreement"] for r in numerical),
            "healthy_fem_runs": sum(h["passed"] for h in health.values()),
            "mesh_converged_cases": sum(c["passed"] for c in convergence.values())},
        "physical_horn_driver_assembly_validated": False,
        "elapsed_seconds": time.perf_counter()-started}
    write_json(out/"audit.json", result)
    print("AUDIT COMPLETE: " + json.dumps(result["totals"]), flush=True)


if __name__ == "__main__":
    main()
