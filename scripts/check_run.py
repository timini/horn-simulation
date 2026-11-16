#!/usr/bin/env python3
"""Check a launcher run directory: status, outputs and recorded provenance.

Usage: check_run.py RUN_DIR --expect PATH [--expect PATH ...] [--params-file FILE]
                    [--min-scored N]
Exits non-zero with one line per failed check. Used by the CI end-to-end job.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


RANKED = "auto/refinement/ranked_results.json"


def check(run_dir, expected, params_file=None, min_scored=None):
    problems = []
    manifest = json.loads((run_dir/"manifest.json").read_text())
    if manifest.get("status") != "completed" or manifest.get("exit_code") != 0:
        problems.append(f"run status {manifest.get('status')} exit {manifest.get('exit_code')}")
    for relative in expected:
        path = run_dir/"outputs"/relative
        if not path.is_file() or path.stat().st_size == 0:
            problems.append(f"missing or empty output: outputs/{relative}")
        elif f"outputs/{relative}" not in manifest.get("output_sha256", {}):
            problems.append(f"output not sealed in manifest: outputs/{relative}")
    sources = manifest.get("source_sha256", {})
    for name in ("main.nf", "nextflow.config", "scripts/run_pipeline.py", "uv.lock"):
        if sources.get(name) != sha256(ROOT/name):
            problems.append(f"source hash missing or wrong for {name}")
    if not (run_dir/"source.tar.gz").is_file():
        problems.append("source archive missing")
    for name, image in manifest.get("image_sources", {}).items():
        if image.get("stale"):
            problems.append(f"image {name} not built from this checkout")
    if set(manifest.get("image_sources", {})) != {"horn-solver", "horn-analysis", "horn-geometry"}:
        problems.append("image source labels not recorded")
    if params_file:
        recorded = manifest.get("input_sha256", {}).get("-params-file", {})
        if list(recorded.values()) != [sha256(Path(params_file).resolve())]:
            problems.append("params file hash not recorded")
    if min_scored is not None:
        ranked = run_dir/"outputs"/RANKED
        scored = json.loads(ranked.read_text()).get("total_scored", 0) if ranked.is_file() else 0
        if scored < min_scored:
            problems.append(f"{RANKED} scored {scored} candidates, expected at least {min_scored}")
    return problems


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--expect", action="append", default=[])
    parser.add_argument("--params-file")
    parser.add_argument("--min-scored", type=int,
                        help=f"require at least N scored candidates in outputs/{RANKED}")
    options = parser.parse_args()
    problems = check(options.run_dir, options.expect, options.params_file, options.min_scored)
    for problem in problems:
        print(f"FAIL: {problem}", file=sys.stderr)
    if not problems:
        print(f"OK: {options.run_dir} completed with {len(options.expect)} expected outputs and provenance")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
