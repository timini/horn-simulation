# Horn Simulation Pipeline

Nextflow-orchestrated FEM acoustic horn simulation with Docker-containerised packages.

## Quick reference

```bash
just                  # list all commands
just build            # build all Docker images (horn-solver, horn-geometry, horn-analysis)
just test-local       # run tests locally (horn-core, horn-geometry, horn-analysis)
just test             # run tests in Docker
just run              # single mode pipeline (default params)
just run-auto         # auto mode (all 7 profiles, driver ranking)
just run-fullauto     # fullauto mode (derives geometry from frequency band)
just clean            # remove Docker images + Nextflow work dirs
```

## Running the pipeline

Always use `-profile docker` (the justfile does this automatically). The solver requires `dolfinx` which only exists in the Docker container.

```bash
# Single mode
just run --profile conical --throat_radius 0.025 --mouth_radius 0.15 --length 0.3

# Auto mode (mid-horn example)
just run-auto --target_f_low 250 --target_f_high 6500 --throat_radius 0.025 --mouth_radius 0.2 --length 0.3

# Fullauto mode (geometry derived from frequency band)
just run-fullauto --target_f_low 500 --target_f_high 4000

# Resume a failed/interrupted run
just run -resume
```

## Project structure

```
main.nf                          # Nextflow pipeline (single/auto/fullauto workflows)
nextflow.config                  # Docker container mappings per process
justfile                         # Build, test, run commands
packages/
  horn-core/                     # Shared data structures, enums, candidate generation (pure Python, no Docker)
  horn-geometry/                 # STEP file generation via gmsh (Docker: horn-geometry)
  horn-solver/                   # FEM solver via dolfinx (Docker: horn-solver)
  horn-analysis/                 # Plots, reports, scoring, rendering (Docker: horn-analysis)
```

## Horn profiles

7 flare profiles: `conical`, `exponential`, `hyperbolic`, `tractrix`, `os`, `lecleach`, `cd`

## Testing

```bash
# Local (fast, no Docker needed — covers horn-core, horn-geometry, horn-analysis)
just test-local

# Single package in Docker
just test-package horn-solver

# All packages in Docker
just test
```

horn-geometry tests need `gmsh` (installed locally). horn-solver tests need `dolfinx` (Docker only).

## Docker images

Build from repo root (context is `.`):
```bash
docker build -t horn-geometry:latest --target production -f ./packages/horn-geometry/Dockerfile .
```

Note: `gmsh` wheels don't have native `aarch64` Linux builds. On Apple Silicon, build with `--platform linux/amd64`.

## Java for Nextflow

Nextflow requires Java 8-22. If the system Java is too new, source `~/.zshrc` which sets the correct `JAVA_CMD`.
