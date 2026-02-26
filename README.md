# Horn Loudspeaker Simulation Pipeline

A containerized pipeline for simulating and analysing horn loudspeaker acoustic performance. Built with FEniCSx (FEM), gmsh, and orchestrated by Nextflow.

## Pipeline Overview

The pipeline takes horn geometry parameters and driver characteristics as input, then runs a full acoustic simulation and produces frequency response plots.

```
                         Horn Parameters
                     (throat, mouth, length)
                              |
                              v
                  +---------------------+
                  |  generate_geometry   |  horn-geometry container
                  |  (gmsh/OCC kernel)   |  Produces STEP file
                  +---------------------+
                              |
                         horn.step
                              |
              +---------------+---------------+
              |               |               |        Frequency band
              v               v               v        parallelisation
        +-----------+   +-----------+   +-----------+
        |  solve    |   |  solve    |   |  solve    |  horn-solver container
        |  band 0   |   |  band 1   |   |  band N   |  FEniCSx FEM
        +-----------+   +-----------+   +-----------+
              |               |               |
         results_0.csv   results_1.csv   results_N.csv
              |               |               |
              +---------------+---------------+
                              |
                              v
                  +---------------------+
                  |   merge_results     |  horn-analysis container
                  |   (pandas concat)   |  Combines band CSVs
                  +---------------------+
                              |
                      final_results.csv
                              |
                              v
                  +---------------------+
                  |   generate_plots    |  horn-analysis container
                  |   (matplotlib)      |  SPL vs Frequency plot
                  +---------------------+
                              |
                    frequency_response.png
```

### Packages

| Package | Purpose | Key Dependencies |
|---------|---------|-----------------|
| `horn-geometry` | Procedural horn geometry generation | gmsh (OpenCASCADE) |
| `horn-solver` | FEM acoustic solving + meshing | FEniCSx/dolfinx, gmsh |
| `horn-analysis` | Result merging, plotting, scoring, ranking | pandas, matplotlib, scipy |
| `horn-core` | Shared data structures (HornParameters, DriverParameters) | numpy |
| `horn-drivers` | Driver database loading and validation | numpy, horn-core |

### What it does

1. **Geometry generation**: Creates a 3D STEP model of a conical horn from parametric inputs (throat radius, mouth radius, length) using gmsh's OpenCASCADE kernel
2. **Meshing + Solving**: Imports the STEP file into gmsh, generates a tetrahedral volume mesh with tagged boundaries (inlet/outlet/walls), then solves the Helmholtz equation at each frequency using FEniCSx. The frequency range is split into N bands and solved in parallel across separate containers
3. **Result merging**: Collects per-band CSV files and concatenates them into a single sorted dataset
4. **Plotting**: Generates a frequency response plot (SPL vs frequency on a log scale)

## Getting Started

### Prerequisites

- [Docker](https://www.docker.com/get-started)
- [just](https://github.com/casey/just) (task runner)
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) (for running the full pipeline)

### Build

```bash
just build
```

### Test

Run all package tests in Docker:
```bash
just test
```

Test a single package:
```bash
just test-package horn-solver
```

### Run the Pipeline

```bash
nextflow run main.nf -profile docker
```

## Parameters

### Single mode (default)

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
| `outdir` | Output directory | `./results` |

### Auto mode (`--mode auto`)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `target_f_low` | Target low frequency (Hz) | `500` |
| `target_f_high` | Target high frequency (Hz) | `4000` |
| `drivers_db` | Path to driver database JSON | `data/drivers.json` |
| `top_n` | Number of top results to return | `10` |

## Example Usage

### Quick development run
```bash
nextflow run main.nf -profile docker --outdir ./results/dev --num_intervals 10
```

### Full resolution analysis
```bash
nextflow run main.nf -profile docker --outdir ./results/full \
    --min_freq 20 --max_freq 20000 --num_intervals 2000
```

### Compare two horn designs
```bash
# Run Horn A
nextflow run main.nf -profile docker --outdir ./results/horn_A \
    --throat_radius 0.05 --mouth_radius 0.2 --length 0.5

# Run Horn B with different length
nextflow run main.nf -profile docker --outdir ./results/horn_B \
    --throat_radius 0.05 --mouth_radius 0.2 --length 0.8

# Plot comparison
python3 -m horn_analysis.compare_horns \
    results/horn_A/final_results.csv "Horn A" \
    results/horn_B/final_results.csv "Horn B" \
    results/comparison.png
```

### Auto-select drivers and horn profiles
```bash
# Run the full auto pipeline: pre-screens drivers, simulates all 3 profiles
# (conical, exponential, hyperbolic), couples each driver post-hoc, and ranks
nextflow run main.nf -profile docker --mode auto \
    --target_f_low 500 --target_f_high 4000 \
    --mouth_radius 0.2 --length 0.5 --top_n 10
```

This runs only 3 FEM simulations (one per profile) and couples all pre-screened drivers in pure Python via the transfer function. Outputs: `auto_ranking.json`, `auto_comparison.png`, `auto_summary.txt`, and a self-contained `auto_report.html` (open in any browser — all plots are base64-embedded).

### CLI tools for individual steps
```bash
# Pre-screen drivers for a target spec
horn-prescreen --drivers-db data/drivers.json \
    --target-f-low 500 --target-f-high 4000 \
    --mouth-radius 0.2 --length 0.5

# Rank drivers against a solver result
horn-rank --solver-csv results/final_results.csv \
    --drivers-db data/drivers.json --throat-radius 0.025 \
    --target-f-low 500 --target-f-high 4000 --top-n 5
```

## Acoustic Modelling

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

- **First-order radiation BC**: The outlet uses a first-order Sommerfeld (Robin) condition, which is accurate for ka < ~3 but increasingly reflective at higher frequencies/larger apertures.
- **Sound-hard walls**: No absorption or damping. Walls are perfectly rigid.
- **Constant air properties**: Temperature and humidity dependence not modelled.
- **Three profiles**: Supports conical, exponential, and hyperbolic. Tractrix profile is not yet available.

## Architecture

This is a monorepo with each package in `packages/`. Each package has its own Dockerfile with multi-stage builds (base → production → test). Tests run inside Docker containers to ensure reproducibility.

The pipeline is orchestrated by Nextflow (`main.nf`), which maps each process to its corresponding Docker container via `nextflow.config`.

## Current Status

The pipeline runs end-to-end. SPL is computed from the outlet-surface RMS pressure with a Robin (radiation impedance) BC at the horn mouth. A validation suite verifies the solver against analytical solutions (straight tube, Webster equation) and provides infrastructure for cross-validation against published data.

## Roadmap

Prioritised capabilities for reaching feature parity with tools like AKABAK. See linked GitHub issues for details.

### Priority 1 — Near-term

- Interior field visualisation (VTK/ParaView export from dolfinx) — [#45](https://github.com/timini/horn-simulation/issues/45)
- Arbitrary STEP file import workflow (user-supplied geometry) — [#47](https://github.com/timini/horn-simulation/issues/47)

### Priority 2 — Medium-term

- Exterior radiation / directivity (Kirchhoff-Helmholtz integral post-processing) — [#48](https://github.com/timini/horn-simulation/issues/48)
- Flexible boundary tagging (replace z-coordinate heuristic with surface naming) — [#49](https://github.com/timini/horn-simulation/issues/49)

### Priority 3 — Longer-term

- Complex geometry support (folded horns, phase plugs, back-loaded horns) — [#51](https://github.com/timini/horn-simulation/issues/51)
- Tractrix horn profile
- Wall absorption / damping materials
- Second-order radiation BC for large ka

### Completed

- HTML report for auto-mode: single self-contained `auto_report.html` with rankings table, 4 embedded plots (coupled SPL, raw profile SPL, impedance, phase/group delay), driver T-S parameter table, and summary cards
- Driver coupling with T-S parameters (transfer function + auto-select pipeline) — [#50](https://github.com/timini/horn-simulation/issues/50)
- Analysis features: impedance plots, scoring, driver DB — [#35](https://github.com/timini/horn-simulation/issues/35)
- Profile diversity: conical, exponential, hyperbolic
