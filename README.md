# Horn Loudspeaker Simulation Pipeline

[![CI](https://github.com/timini/horn-simulation/actions/workflows/ci.yml/badge.svg)](https://github.com/timini/horn-simulation/actions/workflows/ci.yml)
![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)

Give it a frequency band and it designs horns for it. The pipeline generates horn geometries, solves the acoustic field inside each with the finite element method, couples real loudspeaker drivers through their Thiele-Small parameters, and ranks every driver and horn pair. It writes a self-contained HTML report, the horn's STEP geometry and the predicted response, impedance and excursion.

It covers conventional single-driver, rigid, axisymmetric, front-loaded horns. Multiple-entry horns belong in the sister project [MEH Design Studio](https://github.com/timini/meh-design-studio).

**See it first:** the [example site](https://timini.github.io/horn-simulation/) and the [examples folder](examples/README.md) have finished reports you can open without installing anything.

## Status

Experimental. The software works end to end, but its predictions have not been checked against a measured horn and driver assembly. Reports label every result as an experimental prediction, and a run may return no feasible design.

What has been checked:

- The interior solver agrees with an independent finite element code, analytical cavity modes and an independent mode-matching code for the default mouth model.
- The input impedance of one measured exponential horn is matched for ka 1 to 5.
- The fast screen agrees with openwind, and keeps the best designs found by exhaustive FEM on a test grid.

What has not:

- Absolute sound level for a real driver on a real horn, and the order in which drivers are ranked.
- Driver adapters and phase plugs. No shipped driver has a verified interface.
- Directivity, finite baffles and free-standing horns.
- Some comparisons fail and are kept as failures, for example the Hajicek and Neely conical horn standard.

Details: [STATUS](docs/STATUS.md) and [VALIDATION](docs/VALIDATION.md).

## Quick start

### Prerequisites

- [Docker](https://docs.docker.com/get-docker/), running. Give it at least 8 GB of memory.
- [just](https://github.com/casey/just).
- Python 3.10 or newer, as `python3`.
- Java 17 to 22 for Nextflow; 21 is recommended. On macOS: `brew install openjdk@21`. The launcher and `just test-nextflow` find it without changing your shell. To call `nextflow` yourself, point `JAVA_HOME` at it: `export JAVA_HOME=/opt/homebrew/opt/openjdk@21 NXF_JAVA_HOME=$JAVA_HOME`.
- [Nextflow](https://www.nextflow.io/) 24.10.x. The pipeline refuses other versions. To install 24.10.5:

  ```sh
  curl -fL https://github.com/nextflow-io/nextflow/releases/download/v24.10.5/nextflow -o nextflow
  chmod +x nextflow
  mkdir -p ~/.local/bin && mv nextflow ~/.local/bin/
  ```

  `~/.local/bin` is not on the macOS PATH by default. Add it (`export PATH="$HOME/.local/bin:$PATH"` in `~/.zshrc`), or move the file to a folder that is already on your PATH.

- Optional, for the local tests and driver tools: [uv](https://docs.astral.sh/uv/).
- Optional: nf-test. `just test-nextflow` downloads the pinned version itself.

### Build and run

```sh
git clone https://github.com/timini/horn-simulation.git
cd horn-simulation
just build
```

Simulate one horn (about a minute):

```sh
just run --profile exponential --throat_radius 0.025 --mouth_radius 0.1 --length 0.2 \
  --min_freq 500 --max_freq 2000 --num_bands 2
```

Design a horn for a band (a few minutes):

```sh
just run-auto --target_f_low 800 --target_f_high 1600 --drivers_db data/drivers-curated \
  --lem_top_n 3 --refinement_budget 2 --num_bands 2 --num_intervals 101 --mesh_size 0.01
```

This is the [worked example](docs/WORKED_EXAMPLE.md). A small search budget keeps it quick.

Each run gets its own folder under `results/`. Find the latest with:

```sh
just latest-run
```

Open `outputs/single_report.html` (single mode) or `outputs/auto/report/auto_report.html` (auto mode) inside that folder.

### Test

```sh
uv sync --all-packages --extra test
just test-local          # package tests without Docker
.venv/bin/python -m pytest scripts/tests   # launcher, scripts and documentation checks
just test                # package tests in Docker
just test-nextflow       # workflow tests
```

## Modes

**Single mode** (`just run`) simulates one horn you describe, or a STEP file you import, and optionally couples one driver to it.

**Auto mode** (`just run-auto`) takes a target band and optional limits, searches seven profiles and a range of sizes, and ranks driver and horn pairs. `just run-fullauto` is the same command under its old name.

All commands and parameters: [USAGE](docs/USAGE.md).

## Features

- Seven horn profiles from one shared implementation: [PHYSICS, profiles](docs/PHYSICS.md#profiles)
- 3D finite element solver with mesh rules tied to wavelength: [PHYSICS](docs/PHYSICS.md)
- Modal mouth radiation in an infinite baffle, plus local radiation models: [MODAL_APERTURE](docs/MODAL_APERTURE.md)
- Optional thermoviscous wall losses: [LOSS_PHYSICS_VALIDATION](docs/LOSS_PHYSICS_VALIDATION.md)
- Band-driven search with a fast screen, FEM shortlist and bounded refinement: [SEARCH_AND_SCORING](docs/SEARCH_AND_SCORING.md)
- Hard limits and a documented weighted score: [SEARCH_AND_SCORING, score](docs/SEARCH_AND_SCORING.md#score)
- Size caps, fixed dimensions, minimum output and profile restriction: [USAGE, auto search](docs/USAGE.md#auto-search)
- Driver coupling with front and rear chambers: [DRIVERS, chambers](docs/DRIVERS.md#front-and-rear-chambers)
- Driver databases with provenance, a validator and a datasheet entry tool: [DRIVERS](docs/DRIVERS.md)
- Electrical impedance, excursion, power and group delay per candidate: [OUTPUTS](docs/OUTPUTS.md#coupled-csv-columns)
- Experimental directivity and beamwidth: [OUTPUTS, directivity](docs/OUTPUTS.md#directivity-csv)
- True-size STEP files and builder profile CSVs: [OUTPUTS, STEP files](docs/OUTPUTS.md#step-files)
- STEP import with unit handling: [USAGE, imported STEP](docs/USAGE.md#imported-step-geometry)
- Self-contained HTML reports with rejection reasons and margins: [OUTPUTS](docs/OUTPUTS.md)
- Isolated runs with provenance, checked resume and run lookup: [USAGE, launcher](docs/USAGE.md#launcher)
- Cross-software checks against openwind, bempp-cl and Hornresp: [CROSS_SOFTWARE_VALIDATION](docs/CROSS_SOFTWARE_VALIDATION.md)
- Validation records and CI jobs: [VALIDATION](docs/VALIDATION.md)

## Examples

- [Example site](https://timini.github.io/horn-simulation/)
- [examples/](examples/README.md): an 800 to 1600 Hz design and a size-capped 320 to 5000 Hz design, with reports, STEP files and data.
- [Worked example](docs/WORKED_EXAMPLE.md): one auto run from command to result.

## Documentation

The [documentation index](docs/index.md) lists every page. Items from the original brief that were not built are in [DESCOPED](docs/DESCOPED.md). Planned work, including a searchable database of driver and horn trade-offs, is in [ROADMAP](docs/ROADMAP.md).

## Contributing

Open an issue first to discuss a change. [DEVELOPMENT](docs/DEVELOPMENT.md) covers the layout, environments, tests and CI. [AGENTS.md](AGENTS.md) sets the delivery rules.

## License

Code is under the [MIT licence](LICENSE). External reference datasets keep the terms recorded in `data/validation/references.json`.
