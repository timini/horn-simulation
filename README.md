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
- [Python 3.10+](https://www.python.org/downloads/)
- `uv` (can be installed via `pip install uv`)

### Setup

1.  **Activate the virtual environment**:
    This project uses `uv` to manage dependencies for the local development environment (e.g., for IDE support).

    ```bash
    uv venv
    source .venv/bin/activate
    ```

2.  **Install development dependencies**:
    This will install all packages and their dependencies, including `pytest`, into the virtual environment.

    ```bash
    uv pip install -e .
    ```

## Running Tests

All tests should be executed via the `Makefile`. This is the canonical way to ensure tests run in the correct, isolated container environment.

To run the tests for a specific package, use the `make test` command with the `package` variable.

**Example: Run tests for `horn-geometry`**
```bash
make test package=geometry
```

**Example: Run tests for `horn-solver`**
```bash
make test package=solver
``` 