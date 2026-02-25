packages := "horn-solver horn-geometry horn-analysis"

# Display help
default:
    @just --list

# Build all Docker images
build:
    #!/usr/bin/env bash
    set -euo pipefail
    for pkg in {{packages}}; do
        docker build -t "$pkg:latest" --target production -f "./packages/$pkg/Dockerfile" .
    done

# Run all package tests (build then test)
test:
    #!/usr/bin/env bash
    set -euo pipefail
    for pkg in {{packages}}; do
        docker build -t "$pkg:test" --target test -f "./packages/$pkg/Dockerfile" .
    done
    for pkg in {{packages}}; do
        echo "Running tests for $pkg..."
        docker run --rm "$pkg:test" pytest "/app/packages/$pkg/tests"
    done

# Build and test a single package: just test-package horn-solver
test-package pkg:
    docker build -t "{{pkg}}:test" --target test -f "./packages/{{pkg}}/Dockerfile" .
    docker run --rm "{{pkg}}:test" pytest "/app/packages/{{pkg}}/tests" -v

# Run the Nextflow pipeline
run:
    nextflow run main.nf -profile docker

# Run Nextflow tests
test-nextflow:
    nf-test test tests/main.nf.test

# Clean up Docker images and Nextflow files
clean:
    #!/usr/bin/env bash
    set -euo pipefail
    echo "Cleaning up Docker images..."
    for pkg in {{packages}}; do
        docker rmi -f "$pkg:latest" "$pkg:test" || true
    done
    echo "Cleaning up Nextflow files..."
    rm -rf .nextflow work .nextflow.log* trace.txt timeline.html report.html dag.dot dag.png
