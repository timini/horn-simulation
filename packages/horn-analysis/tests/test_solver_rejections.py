"""A per-band solver rejection becomes a candidate rejection, not a failed run."""
import json
import sys

import pytest

from horn_analysis import merge
from horn_analysis.rank_pipeline import rejection_rows


def band_record(path, reason, detail):
    path.write_text(json.dumps({"status": "rejected", "rejection_reasons": [reason],
                                "rejection_detail": detail, "max_freq_hz": 8000.}))
    return path


def test_rejected_band_rejects_the_whole_candidate(tmp_path, monkeypatch):
    first = band_record(tmp_path / "results_os_0001_7.rejected.json", "mesh_node_limit_exceeded", "too many")
    second = band_record(tmp_path / "results_os_0001_6.rejected.json", "throat_area_guard_failed", "faceted")
    output = tmp_path / "os_0001_rejected.json"
    monkeypatch.setattr(sys, "argv", ["merge", "--rejections", str(first), str(second),
                                      "--candidate-id", "os_0001", "--output", str(output)])
    merge.main()
    record = json.loads(output.read_text())
    assert record["candidate_id"] == "os_0001"
    assert record["rejection_reasons"] == ["mesh_node_limit_exceeded", "throat_area_guard_failed"]
    assert len(record["bands"]) == 2

    rows = rejection_rows(record, {"profile": "os", "throat_radius": "0.0125", "mouth_radius": "0.1",
                                   "length": "0.2"}, {"b", "a"})
    assert [r["driver_id"] for r in rows] == ["a", "b"]
    assert all(r["model_feasible"] is False and r["horn_label"] == "os_0001" for r in rows)
    assert rows[0]["throat_radius"] == 0.0125
    assert "mesh_node_limit_exceeded" in rows[0]["rejection_reasons"]


def test_rejection_files_must_name_a_reason(tmp_path):
    bad = tmp_path / "results_x_0.rejected.json"
    bad.write_text(json.dumps({"status": "rejected", "rejection_reasons": []}))
    with pytest.raises(ValueError, match="reason"):
        merge.merge_rejections([bad], candidate_id="x", output=tmp_path / "out.json")


def _band(path, low, high, mouth_p=1.):
    import numpy as np
    import pandas as pd
    pd.DataFrame({"frequency": np.geomspace(low, high, 3), "spl": 94., "z_real": 420., "z_imag": 0.,
                  "schema_version": 2, "inlet_area_m2": .01, "mouth_area_m2": .02,
                  "mouth_u_real": .0001, "mouth_u_imag": 0., "mouth_p_real": mouth_p, "mouth_p_imag": 0.,
                  "radiation_model": "plane_wave", "phasor_convention": "exp(+iwt)_rms",
                  "bc_mode": "dirichlet"}).to_csv(path, index=False)
    return path


def _merge_argv(tmp_path, files, *extra):
    return ["merge", "--files", *map(str, files), "--num-bands", "2", "--min-freq", "100",
            "--max-freq", "1000", "--num-intervals", "6", "--min-points-per-octave", "0",
            "--output", str(tmp_path / "c_results.csv"), *extra]


def test_bands_that_do_not_join_reject_the_candidate_with_the_mismatch_size(tmp_path, monkeypatch):
    edge = 1000 ** .5 * 10
    # Mouth pressure steps from 1.00 to 1.15 Pa at the join, about the 15.6%
    # jump seen between an unresolved 10 elements per wavelength band and its
    # finer neighbour in a real auto run.
    files = [_band(tmp_path / "results_c_0.csv", 100, edge),
             _band(tmp_path / "results_c_1.csv", edge, 1000, mouth_p=1.15)]
    rejection = tmp_path / "c_rejected.json"
    monkeypatch.setattr(sys, "argv", _merge_argv(tmp_path, files, "--candidate-id", "c",
                                                 "--rejection-file", str(rejection)))
    merge.main()
    assert not (tmp_path / "c_results.csv").exists()
    record = json.loads(rejection.read_text())
    assert record["status"] == "rejected"
    assert record["rejection_reasons"] == ["band_boundary_mismatch"]
    worst = record["band_boundary_mismatch"]
    assert worst["quantity"] == "mouth_p"
    assert worst["frequency_hz"] == pytest.approx(edge)
    assert worst["jump"] == pytest.approx(.15 / 1.15)
    assert worst["excess_ratio"] > 1
    assert "mouth_p" in record["rejection_detail"] and "Hz" in record["rejection_detail"]
    rows = rejection_rows(record, {"profile": "os", "length": "0.2"}, {"d"})
    assert rows[0]["rejection_reasons"] == ["band_boundary_mismatch"]


def test_mismatch_still_fails_without_a_rejection_file(tmp_path, monkeypatch):
    edge = 1000 ** .5 * 10
    files = [_band(tmp_path / "results_c_0.csv", 100, edge),
             _band(tmp_path / "results_c_1.csv", edge, 1000, mouth_p=1.15)]
    monkeypatch.setattr(sys, "argv", _merge_argv(tmp_path, files))
    with pytest.raises(merge.BandBoundaryMismatch, match="complex mouth_p mismatch"):
        merge.main()


def test_jump_within_tolerance_still_merges(tmp_path, monkeypatch):
    edge = 1000 ** .5 * 10
    files = [_band(tmp_path / "results_c_0.csv", 100, edge),
             _band(tmp_path / "results_c_1.csv", edge, 1000, mouth_p=1.04)]
    rejection = tmp_path / "c_rejected.json"
    monkeypatch.setattr(sys, "argv", _merge_argv(tmp_path, files, "--candidate-id", "c",
                                                 "--rejection-file", str(rejection)))
    merge.main()
    assert (tmp_path / "c_results.csv").exists() and not rejection.exists()


def test_other_merge_faults_still_stop_the_run(tmp_path, monkeypatch):
    files = [_band(tmp_path / "results_c_0.csv", 100, 1000 ** .5 * 10)]
    monkeypatch.setattr(sys, "argv", _merge_argv(tmp_path, files, "--candidate-id", "c",
                                                 "--rejection-file", str(tmp_path / "r.json")))
    with pytest.raises(ValueError, match="Incomplete simulation"):
        merge.main()
