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
from tempfile import TemporaryDirectory

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


def verify_reference_source(manifest, archive, catalog_path=None):
    catalog_path = catalog_path or Path(__file__).resolve().parents[1]/"data/validation/references.json"
    catalog = json.loads(Path(catalog_path).read_text())
    expected = next(r for r in catalog["references"] if r["id"] == "ernoult-pipe-impedance-v2")
    reference = manifest["reference"]
    fields = ("id", "sha256", "source_url", "download_url", "archive_name", "format", "conditions")
    if any(reference.get(k) != expected[k] for k in fields):
        raise ValueError("Reference identity disagrees with pinned catalog")
    if sha(archive) != expected["sha256"]:
        raise ValueError("Source archive checksum mismatch")


def verify_imported_curves(manifest, archive, root):
    from horn_analysis.reference_data import import_pipe_archive
    # Rebuild from the already authenticated archive: manifest-supplied hashes
    # alone cannot authenticate imported values or prove that no curve is missing.
    with TemporaryDirectory(prefix="horn-reference-check-") as temporary:
        expected = import_pipe_archive(Path(archive).read_bytes(), manifest["reference"], temporary)
        supplied = [{k: v for k, v in c.items() if k != "duplicate_of"}
                    for c in manifest["curves"]]
        if supplied != expected["curves"]:
            raise ValueError("Imported curve manifest disagrees with source archive")
        for curve in expected["curves"]:
            if sha(Path(root)/curve["csv"]) != curve["csv_sha256"]:
                raise ValueError("Imported curve values disagree with source archive")


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


def checked_fem_health(frame, label):
    result = fem_health(frame)
    if not result["passed"]:
        raise RuntimeError(f"Unreliable FEM run {label}; excluded from validation: {result}")
    return result


def checked_fem_grid(frame, count, band, label):
    frequency = frame.frequency.to_numpy()
    if (len(frequency) != count or not np.isfinite(frequency).all()
            or not np.allclose(frequency, np.geomspace(*band, count), rtol=1e-12, atol=0)):
        raise RuntimeError(f"Unexpected FEM frequency grid {label}; excluded from validation")


def checked_mesh_convergence(change, limits, label):
    values = [change["max_db"], change["p95_scaled_complex_error"]]
    if (not np.isfinite(values).all() or values[0] > limits["max_impedance_change_db"]
            or values[1] > limits["p95_scaled_complex_change"]):
        raise RuntimeError(f"Unconverged FEM mesh {label}; excluded from validation: {change}")
    return True


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
    verify_reference_source(manifest, args.archive)
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
        "phase_gate": None, "independent_simulation_scope": "All 40 curves evaluated only at their supplied frequencies within the fixed band. Differing loss and wavefront models are diagnostic comparisons, not identical-model certification. Sparse references do not support full-band conclusions.",
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
    # Reference frequencies are metadata. Do not invent response values between
    # sparse numerical samples (the two 3D references have only 3/11 in-band points).
    numerical_predictions = {}
    numerical_cache = {}
    for curve in manifest["curves"]:
        if curve["reference_kind"] != "independent_simulation":
            continue
        path = root/curve["csv"]
        if sha(path) != curve["csv_sha256"]:
            raise ValueError("Changed numerical reference")
        f = pd.read_csv(path, usecols=["frequency"]).frequency.to_numpy()
        f = f[(f >= 110) & (f <= 3900)]
        if len(f) < 2:
            raise ValueError("Numerical reference needs at least two in-band samples")
        case = NUMERICAL_CASE[curve["configuration"]]
        key = (case, f.tobytes())
        if key not in numerical_cache:
            numerical_cache[key] = prediction(f, case)
        z = numerical_cache[key]
        numerical_predictions[curve["csv"]] = (f, z)
        pd.DataFrame({"frequency": f, "z_real_normalized": z.real,
                      "z_imag_normalized": z.imag}).to_csv(
            out/f"numerical-{curve['csv_sha256'][:16]}-prediction.csv", index=False)
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
    pool = multiprocessing.get_context("spawn").Pool(args.workers)
    try:
        responses = pool.map(run_fem_case, [(out, case) for case in UNIQUE_CASES])
    except BaseException:
        pool.terminate()
        raise
    else:
        # Pool.__exit__ terminates even successful workers, which triggers
        # PETSc's SIGTERM/MPI_Abort handler after valid results are returned.
        pool.close()
    finally:
        pool.join()
    for case, (middle, coarse, fine) in responses:
        fem[case] = middle
        for name, frame in (("middle", middle), ("coarse", coarse), ("fine", fine)):
            count = protocol["fem_samples"] if name == "middle" else protocol["convergence_samples"]
            checked_fem_grid(frame, count, protocol["frequency_band_hz"], f"{case}-{name}")
            health[f"{case}-{name}"] = checked_fem_health(frame, f"{case}-{name}")
        subset = middle.iloc[::4]
        np.testing.assert_allclose(subset.frequency, fine.frequency, rtol=1e-12)
        change = compare(fem_impedance(subset), fem_impedance(fine))
        convergence[case] = {"coarse_vs_fine": compare(fem_impedance(coarse), fem_impedance(fine)),
            "middle_vs_fine": change,
            "passed": checked_mesh_convergence(change, protocol["mesh_convergence_limits"], case),
            "fem_vs_tmm_481": compare(fem_impedance(middle), prediction(middle.frequency.to_numpy(), case))}
    # Freeze every prediction before opening any of the reference curves.
    write_json(out/"prediction_hashes.json", {p.name: sha(p) for p in sorted(out.glob("*.csv"))})
    # Only now parse reference values to authenticate their import. Frequency
    # columns used above are metadata; no response values informed predictions.
    verify_imported_curves(manifest, args.archive, root)
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
            f, prediction_at_samples = numerical_predictions[curve["csv"]]
            selected = frame[frame.frequency.between(110, 3900)]
            np.testing.assert_array_equal(f, selected.frequency.to_numpy())
            z = selected.z_real_pa_s_per_m3.to_numpy()+1j*selected.z_imag_pa_s_per_m3.to_numpy()
            z *= np.pi*ALL_CASES[case][0]**2/(AIR.rho*AIR.c)
            row["sample_count"] = len(f)
            row["sample_frequency_range_hz"] = [float(f[0]), float(f[-1])]
            row["maximum_frequency_gap_hz"] = float(np.diff(f).max())
            row["comparison_scope"] = "Provided reference frequencies only; no full-band certification"
            row["tmm_at_reference_frequencies"] = compare(prediction_at_samples, z)
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
            "independent_simulation_sampled_magnitude_matches": sum(r["tmm_at_reference_frequencies"]["magnitude_agreement"] for r in numerical),
            "healthy_fem_runs": sum(h["passed"] for h in health.values()),
            "mesh_converged_cases": sum(c["passed"] for c in convergence.values())},
        "physical_horn_driver_assembly_validated": False,
        "elapsed_seconds": time.perf_counter()-started}
    write_json(out/"audit.json", result)
    print("AUDIT COMPLETE: " + json.dumps(result["totals"]), flush=True)


if __name__ == "__main__":
    main()
