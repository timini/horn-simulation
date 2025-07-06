# HORN: Horn Loudspeaker Simulation Pipeline

This project is a cloud-deployable horn loudspeaker simulation pipeline based on the architecture described in `docs/Revised Horn Simulation Pipeline Architecture.md`.

## Project Structure

- `docs/`: Contains project documentation.
- `src/horn/`: Contains the main source code for the simulation pipeline.
  - `data_ingestion/`: Handles fetching driver data and parsing datasheets.
  - `geometry/`: Generates 3D horn models.
  - `simulation/`: Runs the core FEM/BEM simulation.
  - `analysis/`: Analyzes simulation results and scores designs.
  - `main.py`: The main orchestrator for the pipeline.
- `tests/`: Contains tests for the source code, following a TDD approach.

## Getting Started

### Prerequisites

This project uses [Poetry](https://python-poetry.org/) for dependency management. Please install it first.

Many of the core simulation tools (`FEniCSx`, `Bempp`, `Gmsh`, `FreeCAD`) are complex system-level dependencies. The recommended way to run this project is via a `Dockerfile` (to be created), which will set up a consistent environment.

### Installation

1. Clone the repository.
2. Install dependencies:
   ```bash
   poetry install
   ```

### Running Tests

To run the test suite:
```bash
poetry run pytest
``` 