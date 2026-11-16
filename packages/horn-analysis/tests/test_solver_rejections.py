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
