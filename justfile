docker_packages := "horn-solver horn-geometry horn-analysis"
local_packages := "horn-core horn-geometry horn-analysis"

# Display help
default:
    @just --list

# Build all Docker images
build:
    #!/usr/bin/env bash
    set -euo pipefail
    for pkg in {{docker_packages}}; do
        docker build --platform linux/amd64 -t "$pkg:latest" --target production -f "./packages/$pkg/Dockerfile" .
    done

# Run all package tests in Docker (build then test)
test:
    #!/usr/bin/env bash
    set -euo pipefail
    for pkg in {{docker_packages}}; do
        docker build --platform linux/amd64 -t "$pkg:test" --target test -f "./packages/$pkg/Dockerfile" .
    done
    for pkg in {{docker_packages}}; do
        echo "Running tests for $pkg..."
        docker run --rm "$pkg:test" pytest "/app/packages/$pkg/tests"
    done

# Build and test a single package in Docker: just test-package horn-solver
test-package pkg:
    docker build --platform linux/amd64 -t "{{pkg}}:test" --target test -f "./packages/{{pkg}}/Dockerfile" .
    docker run --platform linux/amd64 --rm "{{pkg}}:test" pytest "/app/packages/{{pkg}}/tests" -v

# Run tests locally (no Docker — works for horn-core, horn-geometry, horn-analysis)
test-local:
    #!/usr/bin/env bash
    set -euo pipefail
    echo "Running horn-core tests..."
    python -m pytest packages/horn-core/tests/ -v
    echo "Running horn-analysis tests..."
    python -m pytest packages/horn-analysis/tests/ -v
    echo "Running horn-geometry tests..."
    cd packages/horn-geometry && python -m pytest tests/ -v

# Run the Nextflow pipeline (single mode, default params)
run *ARGS:
    nextflow run main.nf -profile docker {{ARGS}}

# Run auto mode (unified optimizer): just run-auto --target_f_low 250 --target_f_high 6500
# Fixed geometry: --mouth_radius 0.15 --length 0.3 (only varies profile)
# Free geometry: omit mouth_radius/length to derive from frequency band
run-auto *ARGS:
    nextflow run main.nf -profile docker --mode auto {{ARGS}}

# Alias for auto mode with all geometry derived (backward compat)
run-fullauto *ARGS:
    nextflow run main.nf -profile docker --mode auto {{ARGS}}

# Run Nextflow tests
test-nextflow:
    nf-test test tests/

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
