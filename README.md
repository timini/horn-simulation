# Horn Loudspeaker Simulation Pipeline

> Design and simulate horn-loaded loudspeakers — from a target frequency band to ranked driver-horn recommendations — using FEM acoustics, a driver database, and automated optimization.

[![CI](https://github.com/timini/horn-simulation/actions/workflows/ci.yml/badge.svg)](https://github.com/timini/horn-simulation/actions/workflows/ci.yml)
![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)

## What is this?

An open-source tool for acoustic horn design. Give it a target frequency band and it will generate horn geometries, solve the Helmholtz equation with FEM, couple the results with real loudspeaker drivers via Thiele-Small parameters, and rank every driver-horn combination — producing a self-contained HTML report you can open in any browser.

Built for audio engineers, acousticians, DIY speaker builders, and researchers.

Browse the [ready-to-open design examples](examples/README.md) for a generated report, winning STEP file, response data and commands for your own band.

**Current status: experimental predictions.** The complete automatic workflow runs, including bounded refinement, failure checks and reports. Its driver/interface and listening-distance output models have not passed independent physical validation. Missing driver evidence is reported as `insufficient_evidence`; an unsuccessful search returns no feasible design. Acoustic STEP exports describe an air volume, not fabrication-ready hardware.

The solver container uses the conservative Nehalem OpenBLAS kernel on x86 to keep the pinned numerical environment consistent across host CPUs. This is a workaround under validation for intermittent CI geometry failures, not a confirmed diagnosis of the upstream cause. Independent basis/area checks and a gross mesh-to-CAD area guard catch corrupt geometry before reporting a response. An explicit `OPENBLAS_CORETYPE` override requires validating that alternative backend.

Use `just run-auto --target_f_low 500 --target_f_high 4000` for an isolated run. The launcher writes a manifest, source snapshot and container hashes under a new `results/<run-id>/` directory. Defaults are 2.83 V RMS, 1 m from the mouth, 6 dB maximum band ripple, a uniform baffled-piston observer, ten screened geometries and six additional refinement evaluations. Override these with `--voltage_rms`, `--observation_distance`, `--max_ripple_db`, `--lem_top_n` and `--refinement_budget`. Fix dimensions or set minimum/maximum search bounds when space is limited.

The [measured exponential-horn comparison](docs/POST_HIXSON_VALIDATION.md) passes its fixed input-impedance limits over `1 ≤ ka ≤ 5` (approximately 201–1,007 Hz at the simulated sound speed). This validates a bounded horn-loading comparison; complete driver/assembly qualification remains outstanding.

Optional thermoviscous losses and finite-flange pipe radiation now have independent impedance checks; see [the measured results and remaining failures](docs/LOSS_PHYSICS_VALIDATION.md). Legacy FEM–BEM horn coupling and directivity are disabled pending a correct exterior model.

See [acoustic assumptions](docs/ACOUSTIC_CONTRACT.md), [existing measurement sources and importer](docs/VALIDATION_DATA.md), and [implementation/validation status](docs/IMPLEMENTATION_STATUS.md). Resume only an unchanged source/data/container/Nextflow snapshot with `python scripts/run_pipeline.py --run-dir results/<run-id> -resume`.


## Scope and sister project

This project focuses on selecting a driver and sizing a conventional single-driver horn for a target frequency band, initially using rigid, axisymmetric, front-loaded horns with characterized driver-to-throat interfaces.

**Multiple-entry horn (MEH) design and simulation are outside this project's scope.** Work on multiple drivers feeding a shared horn through separate entries, their acoustic interaction, and crossover optimization belongs in the sister project: [MEH Design Studio](https://github.com/timini/meh-design-studio).

## Key features

- **7 horn profiles** — conical, exponential, hyperbolic, tractrix, oblate spheroid, Le Cléac’h, and constant directivity
- **FEM Helmholtz solver** — FEniCSx/dolfinx with adaptive meshing and radiation BC
- **Driver database** with Thiele-Small parameter coupling
- **3 operating modes** — single simulation, auto comparison, full-auto design exploration
- **Self-contained HTML reports** — rankings, plots, and driver tables, all base64-embedded
- **Containerized and parallelized** — Docker + Nextflow, frequency bands solved in parallel

## How it works

### Pipeline flow

```mermaid
flowchart LR
    A[Target frequency band] --> B[Pre-screen drivers]
    A --> C[Generate horn geometry]
    C --> D[FEM solver — Helmholtz equation]
    D --> E[Couple with drivers]
    E --> F[Score & rank combinations]
    F --> G[HTML report + plots]
```

### Package architecture

```mermaid
flowchart TB
    subgraph Orchestration
        NF[Nextflow pipeline]
    end
    subgraph Packages
        CORE[horn-core<br/>parameters & data structures]
        GEO[horn-geometry<br/>gmsh STEP generation]
        SOLVER[horn-solver<br/>FEniCSx FEM — Helmholtz]
        ANALYSIS[horn-analysis<br/>scoring, ranking, reports]
        DRIVERS[horn-drivers<br/>driver database]
    end
    NF --> GEO --> SOLVER --> ANALYSIS
    DRIVERS --> ANALYSIS
    CORE --> GEO & SOLVER & ANALYSIS & DRIVERS
```

### Frequency parallelization

```mermaid
flowchart TB
    G[generate_geometry<br/>gmsh/OCC → STEP file] --> S0[solve band 0]
    G --> S1[solve band 1]
    G --> SN[solve band N]
    S0 --> M[merge_results<br/>pandas concat]
    S1 --> M
    SN --> M
    M --> P[generate_plots<br/>SPL vs frequency]
```

## Operating modes

### Single mode (default)

Simulate one horn with explicit geometry parameters and get a frequency response plot.

For imported geometry, supply both `--step_file path/to/horn.step` and `--length` in metres. The supported orientation has the inlet at z = 0 and the outlet at z = length; imported files do not inherit the parametric 0.5 m default. Reports use actual CAD boundary areas and omit inferred circular radii; supply `--horn_3d_png` for your own geometry image, otherwise the report shows an explicit placeholder.

```bash
just run \
    --throat_radius 0.05 --mouth_radius 0.2 --length 0.5
```

### Auto mode

Fix any dimensions you know, derive the others from the band, screen seven profile families, evaluate a shortlist with FEM, and refine within a fixed budget. Ranking uses the requested band and records infeasible combinations and missing evidence.

```bash
just run-auto \
    --target_f_low 500 --target_f_high 4000 \
    --mouth_radius 0.2 --length 0.5 --top_n 10
```

### Fullauto mode

Specify **only** a target frequency band. The system derives horn geometry analytically (mouth radius from cutoff frequency, length from quarter-wave to half-wave), generates a grid of seven profiles and candidate dimensions, screens it analytically, runs FEM on a shortlist, and refines promising dimensions. Reports describe the best evaluated candidates, not a proven global optimum.

```bash
just run-fullauto \
    --target_f_low 500 --target_f_high 4000
```

## Quick start

### Prerequisites

- [Docker](https://www.docker.com/get-started)
- [just](https://github.com/casey/just) (task runner)
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
- **Java 17–22** (required by Nextflow; Java 25+ is not supported). On macOS: `brew install openjdk@21`

### Build

```bash
just build
```

### Run

```bash
just run-fullauto \
    --target_f_low 500 --target_f_high 4000
```

Find the latest completed run with `just latest-run`, then open `outputs/auto/report/auto_report.html` beneath that directory.

### Test

```bash
just test                        # all packages
just test-package horn-solver    # single package
```

## Parameters

<details>
<summary><strong>Single mode (default)</strong></summary>

| Parameter | Description | Default |
|-----------|-------------|---------|
| `throat_radius` | Horn throat radius (m) | `0.05` |
| `mouth_radius` | Horn mouth radius (m) | `0.2` |
| `length` | Horn length along Z-axis (m) | `0.5` |
| `profile` | Horn flare profile | `conical` |
| `min_freq` | Minimum simulation frequency (Hz) | `500` |
| `max_freq` | Maximum simulation frequency (Hz) | `8000` |
| `num_intervals` | Number of frequency steps | `100` |
| `mesh_size` | Target mesh element size (m) | `0.01` |
| `num_bands` | Parallel frequency band jobs | `8` |
| `outdir` | Direct Nextflow output directory | Required for direct invocation; launcher chooses an isolated run |

</details>

<details>
<summary><strong>Auto mode (<code>--mode auto</code>)</strong></summary>

| Parameter | Description | Default |
|-----------|-------------|---------|
| `target_f_low` | Target low frequency (Hz) | `500` |
| `target_f_high` | Target high frequency (Hz) | `4000` |
| `drivers_db` | Path to driver database JSON | `data/drivers.json` |
| `top_n` | Number of top results to return | `10` |


</details>

<details>
<summary><strong>Fullauto mode (<code>--mode fullauto</code>)</strong></summary>

| Parameter | Description | Default |
|-----------|-------------|---------|
| `target_f_low` | Target low frequency (Hz) | `500` |
| `target_f_high` | Target high frequency (Hz) | `4000` |
| `drivers_db` | Path to driver database JSON | `data/drivers.json` |
| `top_n` | Number of top results to return | `10` |
| `num_mouth_radii` | Mouth radius grid points | `3` |
| `num_lengths` | Length grid points | `3` |
| `num_intervals` | Frequency steps per simulation | `100` |
| `num_bands` | Parallel frequency band jobs | `8` |

</details>

## Packages

| Package | Purpose | Key Dependencies |
|---------|---------|-----------------|
| `horn-core` | Shared data structures (`HornParameters`, `DriverParameters`) | numpy |
| `horn-geometry` | Procedural horn geometry generation (STEP files) | gmsh (OpenCASCADE) |
| `horn-solver` | FEM acoustic solving + meshing | FEniCSx/dolfinx, gmsh |
| `horn-analysis` | Result merging, plotting, scoring, ranking, HTML reports | pandas, matplotlib, scipy |
| `horn-drivers` | Driver database loading and validation | numpy, horn-core |

## Acoustic modelling

### Governing equation

The solver computes the steady-state acoustic pressure field inside the horn by solving the **time-harmonic Helmholtz equation**:

```
∇²p + k²p = 0
```

where `p` is complex acoustic pressure, `k = 2πf / c₀` is the wave number, `f` is frequency, and `c₀ = 343 m/s` is the speed of sound in air.

This is discretised using the **Finite Element Method (FEM)** via FEniCSx/dolfinx. The weak (variational) form used is:

```
∫_Ω [∇p · ∇q − k²pq] dx − jk ∫_outlet pq ds = 0
```

where `q` is a test function from a first-order Lagrange (P1) finite element space on the tetrahedral volume mesh. The surface integral implements the first-order Sommerfeld radiation condition at the outlet. The linear system is solved with a direct LU factorisation via PETSc at each frequency.

### Boundary conditions

The horn mesh has three tagged boundary regions, identified automatically by the z-coordinate of each surface's centre of mass:

| Boundary | Tag | Location | Condition |
|----------|-----|----------|-----------|
| **Inlet** (throat) | 2 | z = 0 | Dirichlet: p = 1 Pa (unit driving pressure) |
| **Outlet** (mouth) | 3 | z = length | Robin: ∂p/∂n = −jkp (radiation impedance) |
| **Walls** | 4 | Remaining surfaces | Neumann: ∂p/∂n = 0 (sound-hard walls) |

The wall Neumann condition is natural (satisfied implicitly by the variational form). The inlet Dirichlet condition models a piston driver producing a uniform pressure at the throat. The outlet Robin condition (first-order Sommerfeld radiation condition) is implemented as `a -= jk ∫_outlet p·q ds`, which adds radiation damping and is accurate for ka < ~3.

### Meshing and adaptive element sizing

The STEP geometry is imported into gmsh's OpenCASCADE kernel, which generates a 3D tetrahedral mesh with uniform element size. To ensure accurate wave resolution, the solver enforces a **λ/6 rule**: the element size must not exceed one-sixth of the shortest wavelength being simulated:

```
h_adaptive = c₀ / (6 × f_max)
```

The actual mesh size used is the finer of the user-specified `mesh_size` and `h_adaptive`. Since each frequency band has its own `f_max`, lower bands automatically get coarser (faster) meshes while higher bands get finer meshes.

For example, at `f_max = 8000 Hz`: λ_min = 43 mm, so h_adaptive = 7.1 mm.

### SPL calculation

Sound Pressure Level is computed from the RMS pressure integrated over the **outlet surface** (horn mouth), giving a physically meaningful metric independent of mesh refinement or horn volume:

```
p_rms = √( ∫_outlet |p|² ds  /  A_outlet )

SPL = 20 × log₁₀(p_rms / p_ref)
```

where `p_ref = 20 µPa` is the standard acoustic reference pressure. The outlet area `A_outlet` is computed once before the frequency loop since the mesh is static.

### Frequency sweep

Frequencies are logarithmically spaced using `np.geomspace`, providing finer resolution at lower frequencies where acoustic behaviour changes more rapidly. The total range is split into `num_bands` independent sub-ranges that run in parallel as separate Nextflow processes, each in its own Docker container. Results are merged and sorted by frequency after all bands complete.

### Physical constants

| Constant | Value | Description |
|----------|-------|-------------|
| c₀ | 343.0 m/s | Speed of sound in air (~20 °C) |
| ρ₀ | 1.225 kg/m³ | Air density at sea level |
| p_ref | 20 µPa | SPL reference pressure |

### Limitations

- **Experimental driver recommendations**: driver/interface coupling and absolute listening-distance output still need matched assembly validation. Missing evidence is disclosed in reports.
- **Restricted radiation models**: the optional [modal aperture model](docs/MODAL_APERTURE.md) accounts for nonuniform axisymmetric mouth velocity and passes an independent infinite-baffle comparison for one worked horn. The default local termination remains an approximation. Arbitrary exterior diffraction and directivity are not validated; legacy FEM–BEM horn coupling is disabled.
- **Rigid structure**: optional boundary-layer losses cover a restricted pipe regime, not flexible walls, porous materials or general damping. Measured-reference failures remain visible.
- **Air and driver assumptions**: reference runs record air properties; direct motor coupling requires the default air properties. Nonlinear distortion and broadband thermal performance are not predicted.
- **Profile qualification**: seven profiles are implemented, but validation of simple pipes and conical/exponential search grids does not qualify every profile or frequency band.
- **Acoustic CAD only**: exported STEP files describe the air volume. Walls, mounting and the actual driver interface require mechanical design before fabrication.

## Architecture

This is a monorepo with each package in `packages/`. Each package has its own Dockerfile with multi-stage builds (base → production → test). Tests run inside Docker containers to ensure reproducibility.

The pipeline is orchestrated by Nextflow (`main.nf`), which maps each process to its corresponding Docker container via `nextflow.config`.

## Roadmap

The goal is a dependable **target band → driver and horn dimensions** workflow for conventional single-driver horns. [Issue #81](https://github.com/timini/horn-simulation/issues/81) tracks completion; [the current status and delivery roadmap](docs/SINGLE_HORN_ROADMAP.md) records the evidence, remaining issues and acceptance gates. MEH development belongs in [MEH Design Studio](https://github.com/timini/meh-design-studio).

The automated search, bounded refinement, seven profiles, axial STEP import, isolated runs and reports work. The independent interior/ideal-motor comparison, cavity test, modal radiation comparison, candidate-resolution studies and four-band search audit now have archived passing evidence. See the [worked default example](docs/WORKED_EXAMPLE.md) and [improved-radiation example](docs/MODAL_APERTURE.md). Predictions remain experimental while the following physical gates are outstanding:

1. Obtain complete source/interface and calibration data for two sufficiently documented reference assemblies, and qualify a small driver set (#81). Known FSN classification is corrected; the larger catalogue audit remains #74.
2. Model and independently check the interfaces/rear loads required by those references, extending the numerical domain only where needed (#81).
3. Compare frozen absolute response, impedance and ordering with both measured assemblies. Differences smaller than uncertainty remain ties (#81).
4. Repeat release acceptance with qualified inputs and publish the reference's mechanical/interface drawings alongside acoustic CAD (#81).

Output management (#75) is delivered. Interior field export (#45), flexible boundary tagging (#49), and folded/back-loaded or automatic phase-plug design (#51) remain later extensions unless a selected reference case requires them. Historical closure of exterior-solver issues does not certify the disabled legacy coupling.

## Contributing

Contributions welcome. Please open an issue first to discuss what you'd like to change.

## License

Project code is available under the [MIT licence](LICENSE). External reference datasets retain the source-specific terms recorded in the validation catalog.


### Resumable driver imports

Install the scraper extras with `pip install -e "packages/horn-drivers[scrape]"`.
Run `horn-scrape-drivers --db data/drivers --manufacturers Eminence --state results/scrape-state.json`.
Valid records with the required fields, source provenance and current scraper schema are skipped; use `--refresh` to explicitly refetch them. Older records are refreshed to remove provisional Mms values.
Writes replace individual files atomically, retaining existing permissions or respecting the process umask for new files, so interrupted or failed refreshes preserve the last good record.
Discovery failures, incomplete pagination and failed driver retrieval return a failure status; an all-current resume succeeds.
`--patience-hours` optionally waits through origin outages; exhausting that budget stops the batch instead of restarting the wait for the next driver. Throttling retries remain bounded.

Driver pages must identify the requested URL. Missing Mms is not replaced with dry Mmd, and program power is not converted into an assumed continuous rating. Records lacking essential known parameters are rejected; additional chamber/interface validation is still required before making physical driver recommendations.
Refreshes preserve known driver categories. New records have an unknown category until supported metadata is available; diaphragm area alone cannot distinguish a compression driver from a small cone driver.
Manually enriched interface data survives refresh. Omitted numeric fields survive only with explicit per-field provenance; stale unsourced values are removed. Source-derived nominal-size estimates are kept separate from verified nominal sizes.
During migration, an existing continuous-power rating must have an independent source in `parameter_sources.power_w`; otherwise it is removed because older scraper versions inferred it from program power. A manufacturer's progress is marked incomplete before discovery, so a failed refresh cannot leave an earlier completion flag in place.

See the [independent solver protocol](docs/INDEPENDENT_SOLVER_VALIDATION.md) and [physical assembly reference audit](docs/PHYSICAL_REFERENCE_AUDIT.md) for the current qualification work.
See [run locations, latest-run lookup, resume and cleanup](docs/RUN_MANAGEMENT.md).

Start with the [reproducible 800–1600 Hz worked example](docs/WORKED_EXAMPLE.md), then check a candidate using the [mesh, frequency and loft resolution procedure](docs/CANDIDATE_RESOLUTION.md). Read [cutoff and sweep-bound metric definitions](docs/REPORT_METRICS.md) when interpreting the report.

### September 2026 catalogue audit and 6.5-inch study

[Driver-source research](docs/research/driver-sources-2026-09.md) documents new sources, model revisions, conflicting size labels and the shared catalogue repair. [Source inventory](data/catalogue-sources.json) distinguishes imported records from discovery leads. The catalogue has 29 manufacturer-sourced records, 375 quarantined entries and 1,485 legacy-unverified records; it is not a fully independently verified database. The loader rejects quarantine and unknown inductance, and unverified nominal sizes cannot satisfy size-constrained searches.

```sh
python scripts/audit_driver_catalogue.py --apply
python scripts/refresh_redcatt_drivers.py --size 6.5
python scripts/study_6p5_mid.py --output results/6p5-mid-500-6500
```

The source audit annotates the shared directory and legacy subset without deleting historic records. REDCATT refresh imports selected engineering facts under its research-use terms; it does not infer removals from a filtered feed. Keep the failing secondary-source identity guard enabled. Secondary-source sweeps preserve manufacturer-verified records even with `--refresh`; refresh those through their primary source.

The [completed 500–6,500 Hz study](examples/6p5-mid-500-6500/README.md) uses 26 verified nominal 6.5-inch cone drivers from the shared catalogue, handles factory-sealed rear loads, and reserves a 64 mm central HF housing. It includes predictions and acoustic STEP files. It remains a reduced-order screen; annular FEM, a cone-following phase plug, HF crossover summation and the outer 15-inch horn are unqualified.
