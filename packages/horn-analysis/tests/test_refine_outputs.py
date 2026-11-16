"""Refinement must satisfy its Nextflow outputs even when the solver rejected every candidate."""
import itertools
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]


def _required_refine_outputs():
    text = (ROOT / "main.nf").read_text()
    block = re.search(r"process refine_ranked_candidates \{(.*?)\n    script:", text, re.S).group(1)
    outputs = block.split("output:", 1)[1]
    return [m.group(1) for m in re.finditer(r'path "([^"]+)"(.*)', outputs) if "optional: true" not in m.group(2)]


def _expand_braces(pattern):
    parts = re.split(r"\{([^}]*)\}", pattern)
    choices = [[p] if i % 2 == 0 else p.split(",") for i, p in enumerate(parts)]
    return ["".join(c) for c in itertools.product(*choices)]


def test_every_required_refine_output_exists_when_all_candidates_were_rejected(tmp_path, monkeypatch):
    rejection = {"candidate_id": "conical_0001", "reason": "mesh_node_limit_exceeded"}
    rejected_json = tmp_path / "conical_0001_rejected.json"
    rejected_json.write_text(json.dumps(rejection))
    (tmp_path / "ranked.json").write_text(json.dumps(
        {"results": [], "rejected": [], "status": "no_feasible_design", "total_scored": 0, "total_candidates": 1}))
    (tmp_path / "design.json").write_text(json.dumps({"target_f_low": 800., "target_f_high": 1600.}))
    (tmp_path / "prescreen.json").write_text(json.dumps({"drivers": []}))
    out_root = tmp_path / "work"
    out_root.mkdir()
    monkeypatch.chdir(out_root)
    monkeypatch.setattr(sys, "argv", [
        "refine", "--ranked-json", str(tmp_path / "ranked.json"),
        "--design-json", str(tmp_path / "design.json"),
        "--prescreen-json", str(tmp_path / "prescreen.json"),
        "--drivers-db", str(ROOT / "data/drivers-curated"),
        "--solver-csvs", str(rejected_json),
        "--output-dir", "refinement", "--budget", "6"])
    from horn_analysis import refine
    refine.main()

    required = _required_refine_outputs()
    assert "refinement/ranked_results.json" in required
    for pattern in required:
        matches = [p for alt in _expand_braces(pattern) for p in out_root.glob(alt)]
        assert matches, f"required output {pattern} has no file, so Nextflow would stop the run"
    audit = json.loads((out_root / "refinement/search_audit.json").read_text())
    assert audit["reason"] == "no_feasible_seed"
