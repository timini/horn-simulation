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
    (reference/"manifest.json").write_text(json.dumps({"reference": {"sha256": "0"*64}}))
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
