# Horn Loudspeaker Simulation Pipeline

A containerized pipeline for simulating and analysing horn loudspeaker acoustic performance. Built with FEniCSx (FEM), Bempp (BEM), gmsh, and orchestrated by Nextflow.

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
        |  band 0   |   |  band 1   |   |  band N   |  FEniCSx FEM + Bempp BEM
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
| `horn-solver` | FEM acoustic solving + meshing | FEniCSx/dolfinx, gmsh, bempp-cl |
| `horn-analysis` | Result merging, plotting, scoring | pandas, matplotlib, scipy |
| `horn-core` | Shared data structures (WIP) | — |

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

| Parameter | Description | Default |
|-----------|-------------|---------|
| `throat_radius` | Horn throat radius (m) | `0.05` |
| `mouth_radius` | Horn mouth radius (m) | `0.2` |
| `length` | Horn length along Z-axis (m) | `0.5` |
| `min_freq` | Minimum simulation frequency (Hz) | `500` |
| `max_freq` | Maximum simulation frequency (Hz) | `8000` |
| `num_intervals` | Number of frequency steps | `100` |
| `mesh_size` | Target mesh element size (m) | `0.01` |
| `num_bands` | Parallel frequency band jobs | `8` |
| `outdir` | Output directory | `./results` |

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

## Architecture

This is a monorepo with each package in `packages/`. Each package has its own Dockerfile with multi-stage builds (base → production → test). Tests run inside Docker containers to ensure reproducibility.

The pipeline is orchestrated by Nextflow (`main.nf`), which maps each process to its corresponding Docker container via `nextflow.config`.

## Current Status

The pipeline runs end-to-end for the FEM-only path. The BEM coupling (exterior radiation from the horn mouth) is not yet functional — `bempp.x.dolfinx` coupling layer is not available in any pip release of bempp-cl. The SPL metric currently uses an L2 norm over the full volume as a placeholder. See GitHub issues for the detailed roadmap.
