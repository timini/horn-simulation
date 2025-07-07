# Horn Pipeline Orchestration Plan

This document outlines the current and future strategy for orchestrating the multi-container horn simulation pipeline.

## Phase 1: Local-First Orchestration (Current State)

The current pipeline is orchestrated by a single Python script, `src/horn/main.py`. This script acts as a self-contained workflow engine.

### How It Works

1.  **Sequential Execution**: The script defines the pipeline stages (Data Ingestion, Geometry, Meshing, Solving, Analysis) and runs them in a fixed order.
2.  **Containerization**: The computationally heavy and dependency-rich stages (Geometry, Meshing, Solving) are executed inside dedicated Docker containers (`horn-freecad-app`, `horn-solver-app`).
3.  **Orchestration via `subprocess`**: The main script uses Python's `subprocess.run` to execute `docker run` commands for each containerized stage.
4.  **Data Passing via Volume Mounts**: A temporary directory is created on the host machine for each pipeline run. This directory is volume-mounted (`-v /path/on/host:/data`) into each container, allowing each stage to read the output from the previous one and write its own output.

### Advantages of This Approach

-   **Simplicity**: It requires only Python and Docker to run, with no external platform dependencies.
-   **Portability**: Any user can clone the repository and execute the full pipeline locally.
-   **Testability**: The orchestration logic is easily tested by mocking the `subprocess.run` calls, as demonstrated in `tests/test_e2e.py`.

This model is ideal for development, debugging, and single-run executions.

## Phase 2: Production-Grade Orchestration (Future)

The local-first approach is not designed for large-scale, parallel execution. For that, we will migrate the pipeline to a dedicated workflow orchestration platform.

### The Migration Plan

The core logic of the simulation is encapsulated within the Docker containers. This makes the transition to a formal orchestrator straightforward.

1.  **Develop Core Logic First**: The immediate priority is to complete the implementation of the scientific code within `solver.py` and `horn_generator.py`. We will use the current local orchestrator to test and validate this logic.
2.  **"Lift and Shift" to an Orchestrator**: Once the core simulation is working reliably, we will replace the `main.py` script with a formal workflow definition. The `docker run` commands will be translated into steps in a workflow template.
3.  **Data Handling**: The volume-mounted temporary directory will be replaced by the orchestrator's native artifact-passing mechanism (e.g., S3 buckets, PVCs).

### Candidate Platforms

-   **Argo Workflows**: The primary candidate. It is Kubernetes-native and designed for container-based, multi-step computational workflows like ours.
-   **Kubeflow Pipelines**: An alternative, particularly if the "Analysis" stage evolves to include significant machine learning components.
-   **CI/CD Systems (e.g., GitHub Actions)**: A simpler alternative for event-driven workflows where jobs can pass artifacts to one another.

By containerizing the logic first, we have created portable, reusable components that can be deployed in a simple local environment today and scaled up to a production-grade, distributed environment tomorrow with minimal changes to the core application code. 