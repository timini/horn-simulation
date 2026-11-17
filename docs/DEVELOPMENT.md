# Development

## Layout

```text
main.nf, nextflow.config   Nextflow pipeline (single and auto workflows) and container mapping
justfile                   build, test and run commands
scripts/                   launcher, run lookup, validation scripts and their tests
packages/
  horn-core/               profiles, frequency grid, Webster TMM, geometry designer, driver interface,
                           modal radiation maths (pure Python, no Docker)
  horn-drivers/            driver loader, validator, scraper and horn-add-driver
  horn-geometry/           STEP generation with gmsh/OpenCASCADE (image horn-geometry)
  horn-solver/             dolfinx FEM solver and gmsh meshing (image horn-solver)
  horn-analysis/           prescreen, coupling, scoring, refinement, plots and reports (image horn-analysis)
tests/                     nf-test workflow tests and their data
data/                      driver databases and validation records
validation/                cross-software comparisons and the Hornresp kit
examples/                  preserved example outputs
site/                      GitHub Pages site
docs/                      documentation; docs/history holds dated records
```

The packages form a uv workspace. Profile maths lives only in `horn_core.profiles`.

## Local environment

```sh
uv sync --all-packages --extra test
just test-local
```

A plain `uv sync` installs no workspace packages. The local environment covers everything except the dolfinx solver, which only exists in the `horn-solver` image. horn-geometry and the solver meshing tests need gmsh, which the environment installs from PyPI.

## Docker images

```sh
just build            # horn-solver, horn-geometry and horn-analysis, tag latest
just build dev        # tag dev; run with --image-tag dev
just test             # build the test stages and run every Docker test suite
just test-package horn-solver
```

Images build from the repository root with `--target production` (or `test`). Base images are pinned by multi-arch digest, so builds are native on both linux/amd64 and linux/arm64, including Apple Silicon. The solver and geometry images are based on `dolfinx/dolfinx:v0.8.0` and take gmsh from it; horn-analysis uses `python:3.10-slim`.

Third-party Python packages are installed from `uv.lock` with hashes. The dolfinx-based images keep the base image's numpy (1.26.4), mpi4py and gmsh, because dolfinx, petsc4py and the system gmsh are compiled against them; horn-analysis gets numpy 2.2.6. Analysis code therefore runs under both, since the solver image also carries horn-analysis for refinement. `bempp-cl` is pinned to 0.3.1 in the solver Dockerfile and only its standalone operator tests use it. After changing dependencies, run `uv lock` and commit `uv.lock`; CI checks it.

`just build` labels each image `horn.source_sha256` with a hash of its source files (`scripts/image_source.py`). The launcher compares that label with the checkout and warns when an image is stale.

On x86 the solver entrypoint sets `OPENBLAS_CORETYPE=Nehalem` unless you set it yourself. This keeps one OpenBLAS kernel across host CPUs. It was introduced as a workaround for intermittent CI geometry failures whose upstream cause is not confirmed. Independent throat-area and mesh-to-CAD area checks reject corrupt geometry before a response is reported.

## Nextflow

The pipeline needs Nextflow 24.10.x (`manifest.nextflowVersion` in `nextflow.config`) and Java 17 to 22. The launcher and `scripts/nf_test.sh` find a suitable Java on their own. To run Nextflow or the nf-test binary directly, point `JAVA_HOME` and `NXF_JAVA_HOME` at a Java 21 install, for example `export JAVA_HOME=/opt/homebrew/opt/openjdk@21 NXF_JAVA_HOME=$JAVA_HOME` with Homebrew.

```sh
just test-nextflow
```

This runs nf-test 0.9.5 through `scripts/nf_test.sh`, which downloads it to `.tools/` and checks its sha256 on first use. Set `NF_TEST` to use another copy. The tests are in `tests/*.nf.test` and need the Docker images built.

## Tests

| Command | Covers |
| --- | --- |
| `just test-local` | horn-drivers, horn-core, horn-analysis, horn-geometry and the solver meshing tests, from `.venv` |
| `.venv/bin/python -m pytest scripts/tests` | launcher, run lookup, validation scripts, and the documentation checks |
| `just test` | horn-solver, horn-geometry and horn-analysis inside their test images |
| `just test-nextflow` | the workflow tests |
| `just validate-drivers` | the shipped driver databases |

Solver tests marked `validation` are slow FEM checks; CI skips them on pull requests and runs them on `main`. See [VALIDATION](VALIDATION.md) for the CI jobs.

`scripts/tests/test_docs.py` checks that every relative link in the Markdown files resolves, that every `params.*` in `main.nf` has a row in [USAGE](USAGE.md) with the same default, and that every `just` recipe is listed there. Add the row when you add a parameter.

## Adding a pipeline parameter

1. Declare it in `main.nf` with a default and a short comment.
2. Pass it to the process that needs it, and validate it in the `workflow` block if bad values would give a wrong answer rather than an error.
3. Add a row to [USAGE](USAGE.md).
4. Add an nf-test case if it changes the workflow shape.

## Conventions

- RMS phasors with exp(+iwt), SI units, levels in dB re 20 uPa. See [ACOUSTIC_CONTRACT](ACOUSTIC_CONTRACT.md).
- A failing comparison is a result. Do not loosen a frozen tolerance to make it pass; record the failure.
- Keep validation records under `data/validation/` with the script, source revision and hashes that produced them.
- Read [AGENTS.md](../AGENTS.md) before delivering a change: proportionate validation and one review round.
