# Horn Simulation Pipeline

Nextflow-orchestrated FEM acoustic horn simulation with Docker-containerised packages.

Read [AGENTS.md](AGENTS.md) first. It defines the product goal, proportionate validation and the **one-review-round** delivery policy. Address substantive feedback, verify the fixes and finish the delivery; do not start repeated review cycles.

Documentation starts at [docs/index.md](docs/index.md). Commands and every parameter are in [docs/USAGE.md](docs/USAGE.md); current scope and open work in [docs/STATUS.md](docs/STATUS.md).

## Quick reference

```bash
just                       # list all commands
just build                 # build the Docker images (horn-solver, horn-geometry, horn-analysis)
just test-local            # tests without Docker: horn-drivers, horn-core, horn-analysis, horn-geometry, solver meshing
just test                  # package tests in Docker
just test-package horn-solver
just test-nextflow         # nf-test workflow tests (pinned nf-test via scripts/nf_test.sh)
just validate-drivers      # driver databases; CI requires zero errors
just run                   # single mode (one horn)
just run-auto              # auto mode (band search, all 7 profiles, driver ranking)
just latest-run            # newest completed run folder
just clean                 # remove Docker images and Nextflow files in the repo root
```

`just run-fullauto` and `--mode fullauto` are aliases of auto mode, kept for old commands. There are two modes: single and auto.

## Environment

- Python: `uv sync --all-packages --extra test`. A plain `uv sync` installs no workspace packages. `just test-local` and `just validate-drivers` use `.venv`.
- Nextflow 24.10.x is required (`manifest.nextflowVersion`). It needs Java 17 to 22. The launcher picks a suitable Java itself. To run Nextflow or nf-test directly, set `JAVA_HOME` and `NXF_JAVA_HOME` to a Java 21 install (with Homebrew: `/opt/homebrew/opt/openjdk@21`).
- The solver needs dolfinx, which exists only in the `horn-solver` image. horn-geometry tests and the solver meshing tests need gmsh, which `.venv` has.

## Running the pipeline

Always go through the launcher (`just run`, `just run-auto` or `python3 scripts/run_pipeline.py`). It creates `results/<run-id>/`, records provenance and owns `--outdir`.

```bash
# Single mode
just run --profile conical --throat_radius 0.025 --mouth_radius 0.15 --length 0.3

# Auto mode
just run-auto --target_f_low 800 --target_f_high 1600 --drivers_db data/drivers-curated

# Resume a run
python3 scripts/run_pipeline.py --run-dir results/<run-id> -resume
```

## Project structure

```
main.nf                          # Nextflow pipeline (single and auto workflows)
nextflow.config                  # Docker container mappings per process, Nextflow version pin
justfile                         # build, test and run commands
scripts/                         # launcher, run lookup, validation scripts; tests in scripts/tests
packages/
  horn-core/                     # profiles, frequency grid, Webster TMM, geometry designer, driver interface (no Docker)
  horn-drivers/                  # driver loader, validator, scraper, horn-add-driver (no Docker)
  horn-geometry/                 # STEP generation via gmsh (Docker: horn-geometry)
  horn-solver/                   # FEM solver via dolfinx (Docker: horn-solver)
  horn-analysis/                 # prescreen, coupling, scoring, refinement, plots, reports (Docker: horn-analysis)
data/                            # driver databases and validation records
validation/                      # cross-software comparisons, Hornresp kit
docs/                            # documentation; docs/history holds dated records
```

## Horn profiles

7 flare profiles: `conical`, `exponential`, `hyperbolic` (hypex T=0, catenoidal), `tractrix`, `os`, `scaled_tractrix`, `cd`. `lecleach` is a deprecated alias of `scaled_tractrix`; it was never a Le Cleac'h horn. Profile maths lives only in `horn_core.profiles`.

## Docker images

Build from repo root (context is `.`):

```bash
docker build -t horn-geometry:latest --target production -f ./packages/horn-geometry/Dockerfile .
```

Base images are pinned by multi-arch digest, so builds are native on Apple Silicon (arm64) and amd64. The gmsh-using images take gmsh from the dolfinx base image. `just build` labels images with their source hash; the launcher warns when an image is out of date. Python dependencies come from `uv.lock`; run `uv lock` after changing them.

## Documentation

- `scripts/tests/test_docs.py` fails when a Markdown relative link is broken, a `main.nf` parameter is missing from `docs/USAGE.md` or has a different default there, or a `just` recipe is not listed. Update USAGE.md when adding a parameter or recipe.
- Keep dated plans and superseded reports in `docs/history/`. Do not edit their conclusions.
