# Horn Simulation Pipeline

Nextflow-orchestrated FEM acoustic horn simulation with Docker-containerised packages.

Read [AGENTS.md](AGENTS.md) first. It defines the product goal, proportionate validation and the **one-review-round** delivery policy. Address substantive feedback, verify the fixes and finish the delivery; do not start repeated review cycles.

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

7 flare profiles: `conical`, `exponential`, `hyperbolic` (hypex T=0, catenoidal), `tractrix`, `os`, `scaled_tractrix`, `cd`. `lecleach` is a deprecated alias of `scaled_tractrix`; it was never a Le Cleac'h horn. Profile maths lives only in `horn_core.profiles`.

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

All three images build natively on linux/arm64 and linux/amd64: the gmsh-using images take gmsh from the dolfinx base image rather than a PyPI wheel. `just build` labels images with their source hash; the launcher warns when an image is out of date (see docs/RUN_MANAGEMENT.md).

## Java for Nextflow

Nextflow requires Java 8-22. If the system Java is too new, source `~/.zshrc` which sets the correct `JAVA_CMD`.
