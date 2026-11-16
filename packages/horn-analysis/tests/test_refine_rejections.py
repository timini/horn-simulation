"""Refinement must drop a candidate whose mesh is rejected, not stop the run."""
import json
import sys
import types
from pathlib import Path

import pytest

meshing = pytest.importorskip("horn_solver.meshing")
pytest.importorskip("horn_geometry.generator")

ROOT = Path(__file__).resolve().parents[3]


@pytest.mark.parametrize("error,reason", [
    (meshing.MeshTooLargeError("too many nodes"), "mesh_node_limit_exceeded"),
    (meshing.ThroatAreaGuardError("faceted throat"), "throat_area_guard_failed"),
])
def test_oversized_refinement_mesh_is_a_candidate_rejection(tmp_path, monkeypatch, error, reason):
    calls = []

    def run_simulation_from_step(*args, **kwargs):
        calls.append(args)
        raise error

    fake = types.ModuleType("horn_solver.solver")
    fake.run_simulation_from_step = run_simulation_from_step
    monkeypatch.setitem(sys.modules, "horn_solver.solver", fake)

    seed = dict(horn_label="conical_0001", driver_id="18sound-6nmb420", profile="conical",
                throat_radius=0.025, mouth_radius=0.1, length=0.15, composite_score=0.5,
                model_feasible=True)
    (tmp_path / "ranked.json").write_text(json.dumps(
        {"results": [seed], "rejected": [], "total_scored": 1, "total_candidates": 1}))
    (tmp_path / "design.json").write_text(json.dumps(
        {"target_f_low": 800., "target_f_high": 1600., "sim_freq_range": [566., 2263.],
         "mouth_radius_range": [0.08, 0.12], "length_range": [0.12, 0.18],
         "throat_radius_range": [0.025, 0.025]}))
    (tmp_path / "prescreen.json").write_text(json.dumps(
        {"drivers": ["18sound-6nmb420"], "throat_radii_m": [0.025]}))
    (tmp_path / "conical_0001_results.csv").write_text("frequency\n1\n")
    out = tmp_path / "refinement"
    monkeypatch.setattr(sys, "argv", [
        "refine", "--ranked-json", str(tmp_path / "ranked.json"),
        "--design-json", str(tmp_path / "design.json"),
        "--prescreen-json", str(tmp_path / "prescreen.json"),
        "--drivers-db", str(ROOT / "data/drivers-curated"),
        "--solver-csvs", str(tmp_path / "conical_0001_results.csv"),
        "--output-dir", str(out), "--budget", "1", "--num-bands", "1",
        "--num-frequencies", "4", "--min-points-per-octave", "0"])
    from horn_analysis import refine
    refine.main()
    assert calls, "refinement should have attempted a FEM evaluation"
    data = json.loads((out / "ranked_results.json").read_text())
    assert [r["horn_label"] for r in data["results"]] == ["conical_0001"]
    rejected = [r for r in data["rejected"] if reason in r["rejection_reasons"]]
    assert rejected and all(r["model_feasible"] is False for r in rejected)
    assert data["refinement"]["fem_evaluations"] == 0
