"""Checks that the scientific audit cannot hide invalid or misleading evidence."""
import importlib.util
from pathlib import Path
import sys
import json
import subprocess
import hashlib

import numpy as np
import pandas as pd
import pytest

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location("reference_audit", SCRIPTS/"validate_reference_audit.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
import render_reference_audit as renderer


def test_phase_wrap_does_not_turn_two_degrees_into_358():
    result = audit.compare(np.exp(1j*np.deg2rad([179.])), np.exp(1j*np.deg2rad([-179.])))
    assert result["p95_phase_error_deg"] == pytest.approx(2)


def test_magnitude_pass_does_not_hide_wrong_phase():
    result = audit.compare([1j, 1j], [-1j, -1j])
    assert result["magnitude_agreement"]
    assert result["p95_phase_error_deg"] == pytest.approx(180)
    assert result["p95_scaled_complex_error"] == pytest.approx(1)


@pytest.mark.parametrize("frequency,real", [([120, 3900], [1, 1]),
    ([110, 3800], [1, 1]), ([110, 110, 3900], [1, 1, 1]),
    ([3900, 110], [1, 1]), ([110, 3900], [1, np.nan])])
def test_invalid_reference_cannot_be_extrapolated_or_silently_reordered(frequency, real):
    frame = pd.DataFrame({"frequency": frequency, "real": real, "imag": np.zeros(len(real))})
    with pytest.raises(ValueError, match="entire fixed band"):
        audit.interpolate(frame, np.array([110, 3900]), "real", "imag")


def healthy_frame():
    return pd.DataFrame({"frequency": [110., 3900.], "z_real": [1., 1.], "z_imag": [0., 0.],
        "input_acoustic_power_w": [1., 1.], "mouth_acoustic_power_w": [.5, .5],
        "viscous_wall_power_w": [.3, .3], "thermal_wall_power_w": [.2, .2],
        "relative_residual": [1e-12, 1e-12], "converged_reason": [4, 4], "mesh_cells": [100, 100]})


@pytest.mark.parametrize("column,value", [("z_real", -1), ("relative_residual", 1e-5),
    ("converged_reason", -3), ("input_acoustic_power_w", 2),
    ("viscous_wall_power_w", -.3), ("z_imag", np.inf)])
def test_solver_health_rejects_unreliable_results(column, value):
    frame = healthy_frame()
    assert audit.fem_health(frame)["passed"]
    frame.loc[0, column] = value
    with np.errstate(invalid="ignore"):
        assert not audit.fem_health(frame)["passed"]
        with pytest.raises(RuntimeError, match="excluded from validation"):
            audit.checked_fem_health(frame, "middle/coarse/fine")


def test_all_simulation_configurations_have_explicit_normalization_geometry():
    for case in audit.NUMERICAL_CASE.values():
        assert audit.ALL_CASES[case][0] > 0
    assert audit.ALL_CASES["Cone_U"][3] == "finite_flange"
    assert audit.ALL_CASES["Cone_U"][2] == 0


def test_changed_archive_is_rejected_before_creating_evidence(tmp_path):
    reference = tmp_path/"reference"
    reference.mkdir()
    catalog = json.loads((SCRIPTS.parent/"data/validation/references.json").read_text())
    pinned = next(r for r in catalog["references"] if r["id"] == "ernoult-pipe-impedance-v2")
    (reference/"manifest.json").write_text(json.dumps({"reference": pinned}))
    archive = tmp_path/"reference.zip"
    archive.write_bytes(b"changed archive")
    output = tmp_path/"output"
    result = subprocess.run([sys.executable, str(SCRIPTS/"validate_reference_audit.py"),
        "--reference-dir", str(reference), "--archive", str(archive), "--output-dir", str(output)],
        capture_output=True, text=True)
    assert result.returncode != 0
    assert "Source archive checksum mismatch" in result.stderr
    assert not output.exists()


def test_existing_evidence_cannot_be_overwritten(tmp_path):
    marker = tmp_path/"protocol.json"
    marker.write_bytes(b"preserve existing evidence")
    result = subprocess.run([sys.executable, str(SCRIPTS/"validate_reference_audit.py"),
        "--reference-dir", "unused", "--archive", "unused", "--output-dir", str(tmp_path)],
        capture_output=True, text=True)
    assert result.returncode != 0
    assert "new output directory" in result.stderr
    assert marker.read_bytes() == b"preserve existing evidence"


def test_replaced_prediction_and_hash_manifest_are_rejected_together(tmp_path):
    csv = tmp_path/"prediction.csv"
    manifest = tmp_path/"prediction_hashes.json"
    csv.write_bytes(b"original")
    manifest.write_text(json.dumps({csv.name: renderer.file_sha(csv)}))
    result = {"prediction_hashes_sha256": renderer.file_sha(manifest)}
    renderer.verify_prediction_hashes(tmp_path, result)
    csv.write_bytes(b"different run")
    manifest.write_text(json.dumps({csv.name: renderer.file_sha(csv)}))
    with pytest.raises(ValueError, match="hash manifest"):
        renderer.verify_prediction_hashes(tmp_path, result)


def test_archived_plot_uses_frozen_air_and_rejects_mismatched_csv():
    frame = pd.DataFrame({"z_real": [8.], "z_imag": [16.],
                          "air_rho_kg_m3": [2.], "air_c_m_s": [4.]})
    np.testing.assert_array_equal(renderer.archived_fem_impedance(frame, {"rho": 2., "c": 4.}), [1+2j])
    with pytest.raises(ValueError, match="frozen air"):
        renderer.archived_fem_impedance(frame, {"rho": 1., "c": 4.})


def test_import_coverage_is_verified_from_files(tmp_path):
    directory = tmp_path/"diy"
    directory.mkdir()
    archive = tmp_path/"diy.zip"
    archive.write_bytes(b"archive contents")
    csv = directory/"response.csv"
    csv.write_bytes(b"frequency,spl\n100,90\n")
    dataset = {"reference": {"id": "diy", "archive_name": archive.name,
        "sha256": renderer.file_sha(archive), "format": "zip_frd_zma"},
        "curves": [{"csv": csv.name, "csv_sha256": renderer.file_sha(csv),
                    "source_sha256": hashlib.sha256(b"source curve").hexdigest()}]}
    (directory/"manifest.json").write_text(json.dumps(dataset))
    path = tmp_path/"inventory.json"
    path.write_text(json.dumps({"datasets": [dataset], "unique_curve_files": 1}))
    result = renderer.verify_import_inventory(path)
    assert (result["archives"], result["curves"], result["diy_curves"]) == (1, 1, 1)
    csv.write_bytes(b"changed")
    with pytest.raises(ValueError, match="Changed imported curve"):
        renderer.verify_import_inventory(path)


@pytest.mark.parametrize("change", [{"max_db": .6, "p95_scaled_complex_error": .01},
    {"max_db": .1, "p95_scaled_complex_error": .06},
    {"max_db": np.nan, "p95_scaled_complex_error": .01}])
def test_unconverged_mesh_cannot_contribute_to_validation(change):
    limits = {"max_impedance_change_db": .5, "p95_scaled_complex_change": .05}
    assert audit.checked_mesh_convergence({"max_db": .1, "p95_scaled_complex_error": .01}, limits, "good")
    with pytest.raises(RuntimeError, match="excluded from validation"):
        audit.checked_mesh_convergence(change, limits, "bad")


def test_old_search_configuration_cannot_be_presented_as_current_benchmark():
    with pytest.raises(ValueError, match="Unexpected search benchmark configuration"):
        renderer.validate_search_configuration({}, {"geometry_count": 12, "driver_count": 2})


def test_search_case_count_must_match_frozen_configuration():
    protocol = {"benchmark": "finite_grid_manufacturer_motors_v2",
        "bands": [[800, 1600, "development"], [1000, 1200, "held_out"],
                  [1200, 2000, "held_out"], [600, 1000, "held_out"]],
        "geometry_count": 20, "driver_count": 3, "frequencies": 100, "mesh_size": .008,
        "shortlist_budget": 10, "score_regret_limit": .02, "top_ten_recall_limit": .9,
        "minimum_feasible_pairs_per_case": 10, "physical_validation_passed": False,
        "sim_band": [600/2**.5, 2000*2**.5]}
    result = {**protocol, "candidates": [{}]*20, "drivers": [{}]*3, "cases": [{}]*3}
    with pytest.raises(ValueError, match="case/count mismatch"):
        renderer.validate_search_configuration(result, protocol)


def resonance_protocol():
    return {"geometry": {"length_m": .535, "throat_diameter_m": .018, "mouth_diameter_m": .08},
        "measured_sixth_resonance_hz": 1712.,
        "assumed_air": {"c": 343., "rho": 1.225, "gamma": 1.4, "viscosity": 1.81e-5, "conductivity": .0257, "heat_capacity": 1006.},
        "loss_model": "boundary_layer", "radiation_model": "finite_flange",
        "assumed_flange_width_m": 0., "frequency_step_hz": .1,
        "source_url": "https://doi.org/10.5050/KSNVE.2014.24.7.537",
        "fitted_parameters": [], "physical_assembly_validated": False}


def test_different_horn_cannot_use_the_published_resonance_label(tmp_path):
    protocol = resonance_protocol()
    protocol["geometry"]["length_m"] = .136
    (tmp_path/"protocol.json").write_text(json.dumps(protocol))
    (tmp_path/"comparison.json").write_text(json.dumps(protocol))
    with pytest.raises(ValueError, match="Unexpected published resonance protocol"):
        renderer.load_verified_resonance(tmp_path/"comparison.json")


def test_changed_auxiliary_resonance_prediction_is_rejected(tmp_path):
    protocol = resonance_protocol()
    result = {**protocol, "predictions": {"400": {"prediction_sha256": "0"*64}, "800": {}}}
    (tmp_path/"protocol.json").write_text(json.dumps(protocol))
    (tmp_path/"comparison.json").write_text(json.dumps(result))
    (tmp_path/"prediction-400.csv").write_bytes(b"different prediction")
    with pytest.raises(ValueError, match="Changed resonance prediction"):
        renderer.load_verified_resonance(tmp_path/"comparison.json")


def test_modified_aggregate_cannot_override_independently_pinned_digest(tmp_path):
    path = tmp_path/"audit.json"
    path.write_text(json.dumps({"measured_passes": 237}))
    pinned = renderer.file_sha(path)
    assert renderer.load_pinned_json(path, pinned)["measured_passes"] == 237
    path.write_text(json.dumps({"measured_passes": 299}))
    with pytest.raises(ValueError, match="independent pinned digest"):
        renderer.load_pinned_json(path, pinned)


def test_exhaustive_search_cannot_duplicate_a_pair_to_hide_an_omission():
    rows = [{"horn_label": c, "driver_id": d} for c in ("a", "b") for d in ("x", "y")]
    case = {"exhaustive": rows, "screening": {"filtered_candidate_ids": ["a"]}}
    renderer.validate_search_pairs(case, {"a", "b"}, {"x", "y"}, 1)
    rows[-1] = rows[0]
    with pytest.raises(ValueError, match="Cartesian product"):
        renderer.validate_search_pairs(case, {"a", "b"}, {"x", "y"}, 1)


@pytest.mark.parametrize("shortlist", [["a", "b"], ["a", "a"], ["unknown"], []])
def test_search_rejects_overbudget_duplicate_unknown_or_empty_shortlist(shortlist):
    case = {"exhaustive": [{"horn_label": c, "driver_id": "x"} for c in ("a", "b")],
            "screening": {"filtered_candidate_ids": shortlist}}
    with pytest.raises(ValueError, match="shortlist"):
        renderer.validate_search_pairs(case, {"a", "b"}, {"x"}, 1)


@pytest.mark.parametrize("key,value", [("assumed_air", {"c": 350.}),
    ("loss_model", "lossless"), ("radiation_model", "flanged_piston"),
    ("assumed_flange_width_m", .01)])
def test_resonance_requires_frozen_physical_model(tmp_path, key, value):
    protocol = resonance_protocol()
    protocol[key] = value
    (tmp_path/"protocol.json").write_text(json.dumps(protocol))
    (tmp_path/"comparison.json").write_text(json.dumps(protocol))
    with pytest.raises(ValueError, match="Unexpected published resonance protocol"):
        renderer.load_verified_resonance(tmp_path/"comparison.json")


def test_replaced_archive_and_manifest_cannot_replace_catalog_identity(tmp_path):
    catalog = json.loads((Path(__file__).resolve().parents[2]/"data/validation/references.json").read_text())
    reference = next(r for r in catalog["references"] if r["id"] == "ernoult-pipe-impedance-v2")
    archive = tmp_path/"replacement.zip"
    archive.write_bytes(b"replacement")
    reference["sha256"] = audit.sha(archive)
    with pytest.raises(ValueError, match="pinned catalog"):
        audit.verify_reference_source({"reference": reference}, archive)
