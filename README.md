# Horn Loudspeaker Simulation Pipeline

This project provides a complete, containerized pipeline for the simulation and analysis of horn loudspeakers. It includes packages for procedural geometry generation, FEM-based acoustic solving, and a Dagster-based orchestration pipeline.

## Architecture

This project is structured as a monorepo managed by `uv`, with a strong emphasis on containerization and isolated package testing.

- **Packages**: The core logic is divided into distinct Python packages located in the `packages/` directory.
  - `horn-core`: Common data structures and utilities shared across packages.
  - `horn-geometry`: Procedural generation of horn geometries using `gmsh` and `FreeCAD`.
  - `horn-solver`: FEM acoustic solver using `FEniCSx`.
  - `horn-dagster`: A Dagster pipeline for orchestrating the geometry and solver services.
- **Containerization**: Each package has its own `Dockerfile` with a multi-stage build process. This creates lean, optimized images for production and separate environments for testing.
- **Testing**: All tests are run via a top-level `Makefile`. This ensures that tests are executed in a clean, consistent, and isolated Docker environment, preventing "it works on my machine" issues.

## Getting Started

### Prerequisites

- [Docker](https://www.docker.com/get-started)
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)

### Setup

1.  **Build the Docker images**:
    The pipeline relies on Docker images for each of the core packages. Build them using the `Makefile`.

    ```bash
    make build
    ```

## Running the Pipeline

The simulation is orchestrated using Nextflow. The main script is `main.nf`. You can customize the simulation by passing parameters on the command line.

### Parameters

The following parameters can be set using the `--<param_name>` syntax (e.g., `--min_freq 20`).

| Parameter          | Description                                                              | Default Value |
| ------------------ | ------------------------------------------------------------------------ | ------------- |
| `outdir`           | The directory where all output files will be saved.                      | (required)    |
| **Horn Geometry**  |                                                                          |               |
| `throat_radius`    | Radius of the horn's throat (the narrow end) in meters.                  | `0.05`        |
| `mouth_radius`     | Radius of the horn's mouth (the wide end) in meters.                     | `0.2`         |
| `length`           | Length of the horn along the Z-axis in meters.                           | `0.5`         |
| **Simulation**     |                                                                          |               |
| `min_freq`         | Minimum frequency for the simulation sweep in Hz.                        | `500`         |
| `max_freq`         | Maximum frequency for the simulation sweep in Hz.                        | `8000`        |
| `num_intervals`    | The number of frequency steps to simulate in the sweep.                  | `100`         |
| `mesh_size`        | The target mesh element size in meters for the FEM simulation.           | `0.01`        |
| **Execution**      |                                                                          |               |
| `num_bands`        | Number of parallel jobs to split the frequency sweep into.               | `8`           |

### Example Usage

#### 1. Quick Development Run
This is useful for quickly testing the pipeline with low-resolution settings.

```bash
nextflow run main.nf --outdir ./results/dev_test --num_intervals 10
```

#### 2. Full Resolution Analysis
This runs a detailed simulation suitable for final analysis. This may take a significant amount of time.

```bash
nextflow run main.nf --outdir ./results/full_run \
    --min_freq 20 \
    --max_freq 20000 \
    --num_intervals 2000
```

#### 3. Comparing Two Horn Designs
You can run the pipeline multiple times with different geometric parameters and save the results to different directories.

First, run the simulation for "Horn A":
```bash
nextflow run main.nf --outdir ./results/horn_A \
    --throat_radius 0.05 --mouth_radius 0.2 --length 0.5
```

Next, run the simulation for "Horn B" with a different length:
```bash
nextflow run main.nf --outdir ./results/horn_B \
    --throat_radius 0.05 --mouth_radius 0.2 --length 0.8
```

After the runs are complete, you can use the `compare_horns.py` script to plot both results on the same graph:

```bash
python3 -m horn_analysis.compare_horns \
    results/horn_A/final_results.csv "Horn A" \
    results/horn_B/final_results.csv "Horn B" \
    results/comparison.png
```

This will generate `comparison.png` in the `results/` directory. 