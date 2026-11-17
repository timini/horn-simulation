docker_packages := "horn-solver horn-geometry horn-analysis"
local_packages := "horn-core horn-geometry horn-analysis"

# Display help
default:
    @just --list

# Build all Docker images, labelled with the hash of their sources: just build [tag]
build tag="latest":
    #!/usr/bin/env bash
    set -euo pipefail
    for pkg in {{docker_packages}}; do
        docker build -t "$pkg:{{tag}}" --target production \
            --label "horn.source_sha256=$(python3 scripts/image_source.py "$pkg")" \
            -f "./packages/$pkg/Dockerfile" .
    done

# Run all package tests in Docker (build then test)
test:
    #!/usr/bin/env bash
    set -euo pipefail
    for pkg in {{docker_packages}}; do
        docker build -t "$pkg:test" --target test -f "./packages/$pkg/Dockerfile" .
    done
    for pkg in {{docker_packages}}; do
        echo "Running tests for $pkg..."
        docker run --rm "$pkg:test" pytest "/app/packages/$pkg/tests"
    done

# Build and test a single package in Docker: just test-package horn-solver
test-package pkg:
    docker build -t "{{pkg}}:test" --target test -f "./packages/{{pkg}}/Dockerfile" .
    docker run --rm "{{pkg}}:test" pytest "/app/packages/{{pkg}}/tests" -v

# Run the tests that need no Docker (horn-drivers, horn-core, horn-analysis, horn-geometry, solver meshing)
test-local:
    #!/usr/bin/env bash
    set -euo pipefail
    echo "Running horn-drivers tests..."
    .venv/bin/python -m pytest packages/horn-drivers/tests/ -v
    echo "Running horn-core tests..."
    .venv/bin/python -m pytest packages/horn-core/tests/ -v
    echo "Running horn-analysis tests..."
    .venv/bin/python -m pytest packages/horn-analysis/tests/ -v
    echo "Running horn-geometry tests..."
    .venv/bin/python -m pytest packages/horn-geometry/tests/ -v
    echo "Running horn-solver meshing tests (gmsh only, no dolfinx)..."
    .venv/bin/python -m pytest packages/horn-solver/tests/test_meshing.py -v

# Validate the shipped driver databases (CI requires zero errors)
validate-drivers:
    #!/usr/bin/env bash
    set -euo pipefail
    .venv/bin/horn-validate-drivers --db data/drivers | tail -3
    .venv/bin/horn-validate-drivers --db data/drivers-curated | tail -3

# Regenerate docs/DRIVER_COVERAGE.md from the driver and price files (no network)
driver-coverage:
    .venv/bin/python scripts/driver_coverage.py --out docs/DRIVER_COVERAGE.md

# Run the pipeline through the launcher in single mode: just run --throat_radius 0.025 --mouth_radius 0.1 --length 0.2
run *ARGS:
    python3 scripts/run_pipeline.py {{ARGS}}

# Run auto mode (band search and driver ranking): just run-auto --target_f_low 800 --target_f_high 1600
run-auto *ARGS:
    python3 scripts/run_pipeline.py --mode auto {{ARGS}}

# Same as run-auto; kept for older commands
run-fullauto *ARGS:
    python3 scripts/run_pipeline.py --mode auto {{ARGS}}

# Run Nextflow tests
test-nextflow *ARGS:
    scripts/nf_test.sh test tests/ {{ARGS}}

# Clean up Docker images and Nextflow files
clean:
    #!/usr/bin/env bash
    set -euo pipefail
    echo "Cleaning up Docker images..."
    for pkg in {{docker_packages}}; do
        docker rmi -f "$pkg:latest" "$pkg:test" || true
    done
    echo "Cleaning up Nextflow files..."
    rm -rf .nextflow work .nextflow.log* trace.txt timeline.html report.html dag.dot dag.png

# Locate the most recently completed run (override with --status running/failed/any).
latest-run *args:
    python3 scripts/runs.py latest {{args}}

# Read-only inventory of completed, failed, running and historical folders.
runs *args:
    python3 scripts/runs.py list {{args}}

# Cross-software checks (openwind, bempp-cl, Hornresp sheets); see docs/CROSS_SOFTWARE_VALIDATION.md
validate-cross-software:
    #!/usr/bin/env bash
    set -euo pipefail
    if [ ! -x .venv-cross-software/bin/python ]; then
        uv venv .venv-cross-software --python 3.12
        VIRTUAL_ENV=.venv-cross-software uv pip install -r validation/cross_software/requirements.txt
    fi
    export PYTHONPATH=packages/horn-core/src NUMBA_NUM_THREADS="${NUMBA_NUM_THREADS:-4}"
    .venv-cross-software/bin/python -W ignore validation/cross_software/openwind_compare.py
    .venv-cross-software/bin/python -W ignore validation/cross_software/bempp_compare.py
    .venv/bin/python validation/hornresp/hornresp_records.py

# Compare archived horn-solver FEM runs with openwind: just validate-cross-software-fem results/<conical-run> results/<exponential-run>
validate-cross-software-fem *RUNS:
    PYTHONPATH=packages/horn-core/src .venv-cross-software/bin/python -W ignore validation/cross_software/fem_compare.py {{RUNS}}
