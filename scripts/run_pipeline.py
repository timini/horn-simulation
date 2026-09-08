#!/usr/bin/env python3
"""Launch an isolated Nextflow run with source, data and container provenance.

Use --run-dir PATH -resume to resume exactly that run. Arguments after the
wrapper options are forwarded to Nextflow as argument tokens, never a shell.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
import re
import shutil
from pathlib import Path
import subprocess
import sys
import tarfile
import uuid

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_hashes(run_dir=None):
    files = subprocess.check_output(["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"], cwd=ROOT).decode().split("\0")
    return {name: digest(ROOT/name) for name in sorted(set(files))
            if name and not name.startswith((".nf-test/", ".nextflow/", "results/", "work/"))
            and (ROOT/name).is_file()
            and (run_dir is None or not (ROOT/name).resolve().is_relative_to(Path(run_dir).resolve()))}


def inspect_images():
    images = {}
    for name in ("horn-solver", "horn-analysis", "horn-geometry"):
        data = json.loads(subprocess.check_output(["docker", "image", "inspect", name+":latest"]))[0]
        images[name] = {"id": data["Id"], "repo_digests": data.get("RepoDigests", [])}
    return images


def java_environment():
    """Select an available Java 17–22 for Nextflow without global changes."""
    environment = os.environ.copy()
    candidates = []
    for key in ("NXF_JAVA_HOME", "JAVA_HOME"):
        if environment.get(key): candidates.append(Path(environment[key])/"bin/java")
    if environment.get("JAVA_CMD"): candidates.append(Path(environment["JAVA_CMD"]))
    if shutil.which("java"): candidates.append(Path(shutil.which("java")))
    if sys.platform == "darwin":
        for version in ("21", "17"):
            for prefix in ("/opt/homebrew", "/usr/local"):
                candidates.append(Path(prefix)/f"opt/openjdk@{version}/libexec/openjdk.jdk/Contents/Home/bin/java")
            found = subprocess.run(["/usr/libexec/java_home", "-F", "-v", version], capture_output=True, text=True)
            if found.returncode == 0: candidates.append(Path(found.stdout.strip())/"bin/java")
    for executable in candidates:
        if not executable.is_file(): continue
        result = subprocess.run([str(executable), "-version"], capture_output=True, text=True)
        match = re.search(r'version "(\d+)', result.stderr + result.stdout)
        if result.returncode == 0 and match and 17 <= int(match[1]) <= 22:
            environment["NXF_JAVA_HOME"] = str(executable.resolve().parent.parent)
            return environment
    return environment  # Nextflow reports its normal diagnostic if unavailable.


def nextflow_identity(environment):
    """Identify both the launcher and the engine it actually selects."""
    executable = shutil.which("nextflow", path=environment.get("PATH"))
    if not executable:
        raise ValueError("Nextflow executable unavailable")
    executable = Path(executable).resolve()
    output = subprocess.check_output([str(executable), "-version"], env=environment,
                                     text=True, stderr=subprocess.STDOUT)
    match = re.search(r"version\s+([\w.+-]+)\s+build\s+(\d+)", output)
    if not match:
        raise ValueError("Cannot identify the Nextflow engine version/build")
    return {"executable": str(executable), "launcher_sha256": digest(executable),
            "version": match[1], "build": match[2]}


def main():
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument("--run-dir", type=Path)
    options, forwarded = parser.parse_known_args()
    if "--outdir" in forwarded or any(x.startswith("--outdir=") for x in forwarded):
        parser.error("Use --run-dir with this launcher; it owns the output directory")
    resume = "-resume" in forwarded
    if resume and options.run_dir is None:
        parser.error("Resume requires --run-dir pointing to an existing run")
    run_dir = (options.run_dir or ROOT/"results"/(datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid.uuid4().hex[:8])).resolve()
    manifest_path = run_dir/"manifest.json"
    if any(x in forwarded for x in ("-work-dir", "-w", "-log", "-c", "-C")):
        parser.error("The launcher owns work, log and container configuration paths")
    images, hashes = inspect_images(), source_hashes(run_dir)
    run_environment = java_environment()
    engine = nextflow_identity(run_environment)
    # Freeze the engine selected by a mutable Nextflow launcher for this run.
    run_environment["NXF_VER"] = engine["version"]
    previous = None
    resume_tokens = []
    if resume:
        previous = json.loads(manifest_path.read_text())
        if previous.get("nextflow_engine") != engine:
            parser.error("Nextflow engine or launcher changed or was not recorded; start a new run")
        if previous["source_sha256"] != hashes or previous["containers"] != images:
            parser.error("Source/data or container images changed; start a new run")
        old_args = previous["arguments"]
        new_args = [x for x in forwarded if x != "-resume"]
        if new_args and new_args != old_args:
            parser.error("Resume arguments differ from the original specification")
        session_id = previous.get("nextflow_session_id")
        if not session_id:
            parser.error("Original Nextflow session ID unavailable; start a new run")
        forwarded = old_args
        resume_tokens = ["-resume", session_id]
    elif run_dir.exists():
        parser.error("Output directory already exists; use -resume or a new directory")
    input_hashes = {}
    for key, default in (("--drivers_db", "data/drivers"), ("--step_file", None), ("--horn_3d_png", None)):
        value = default
        for i, token in enumerate(forwarded):
            if token == key and i+1 < len(forwarded): value = forwarded[i+1]
            elif token.startswith(key+"="): value = token.split("=",1)[1]
        if value:
            path = (ROOT/value).resolve()
            paths = sorted(path.rglob("*")) if path.is_dir() else [path]
            input_hashes[key] = {str(f): digest(f) for f in paths if f.is_file()}
    if resume and previous.get("input_sha256") != input_hashes:
        parser.error("External input files changed; start a new run")
    run_dir.mkdir(parents=True, exist_ok=True)
    # Pin actual image IDs for this execution, even if mutable tags later move.
    config = run_dir/"containers.config"
    mappings = {
        "horn-geometry": "generate_geometry|generate_candidate_geometry",
        "horn-solver": "run_simulation|run_candidate_simulation|run_simulation_directivity|refine_ranked_candidates",
        "horn-analysis": "merge_results|extract_kpis|generate_plots|generate_impedance_plot|generate_phase_plot|generate_dashboard|couple_with_driver|render_horn_3d|merge_directivity_results|generate_directivity_plots|generate_single_report.*|prescreen_drivers|report_no_drivers|derive_auto_geometry|lem_prescreen|merge_candidate_results|score_and_rank|generate_auto_report",
    }
    config.write_text("process {\n"+"\n".join(f"  withName: /{pattern}/ {{ container = '{images[name]['id']}' }}" for name,pattern in mappings.items())+"\n}\n")
    command = [engine["executable"], "-log", str(run_dir/"nextflow.log"), "run", str(ROOT/"main.nf"), "-profile", "docker", "-c", str(config), "-work-dir", str(run_dir/"work"), *forwarded, *resume_tokens, "--outdir", str(run_dir/"outputs")]
    manifest = {
        "schema_version": 1, "started_at": datetime.now(timezone.utc).isoformat(),
        "source_revision": subprocess.check_output(["git", "rev-parse", "HEAD"],cwd=ROOT,text=True).strip(),
        "source_sha256": hashes, "input_sha256": input_hashes, "containers": images, "arguments": forwarded,
        "command": command, "status": "running",
        "nextflow_engine": engine,
        "nextflow_java_home": run_environment.get("NXF_JAVA_HOME"), "validation_status": "experimental",
        "previous_attempt": previous.get("started_at") if previous else None,
    }
    manifest_path.write_text(json.dumps(manifest,indent=2))
    (run_dir/"source.patch").write_bytes(subprocess.check_output(["git","diff","--binary","HEAD"],cwd=ROOT))
    if not resume:
        with tarfile.open(run_dir/"source.tar.gz", "w:gz") as archive:
            for name in hashes:
                archive.add(ROOT/name, arcname=name, recursive=False)
    print(f"Run directory: {run_dir}", flush=True)
    code = 1
    try:
        code = subprocess.call(command,cwd=ROOT,env=run_environment)
    finally:
        log_path = run_dir/"nextflow.log"
        sessions = re.findall(r"Session UUID: ([a-f0-9-]+)", log_path.read_text()) if log_path.exists() else []
        manifest["nextflow_session_id"] = sessions[-1] if sessions else (previous or {}).get("nextflow_session_id")
        manifest["exit_code"] = code
        manifest["status"] = "completed" if code == 0 else "failed"
        manifest["finished_at"] = datetime.now(timezone.utc).isoformat()
        manifest_path.write_text(json.dumps(manifest,indent=2))
    return code


if __name__ == "__main__":
    sys.exit(main())
