# Horn Loudspeaker Simulation and Optimization Project

## 1. Project Vision

The ultimate goal of this project is to create a powerful, automated simulation pipeline to optimize the design of horn loudspeakers. By combining a specific driver unit with a computationally generated horn geometry, we can analyze the acoustic performance and iterate on the design to achieve optimal results.

The long-term vision includes integrating a database of loudspeaker drivers to find the best possible driver/horn combination for a given set of performance criteria.

## 2. Simulation Pipeline Architecture

The project is designed as a multi-stage pipeline, where each stage is a distinct, containerized microservice responsible for one part of the simulation process. This architecture ensures separation of concerns and allows for swapping out implementations (e.g., different physics solvers) without affecting the entire system.

The pipeline consists of the following four stages:

![Pipeline Diagram](https://i.imgur.com/your_diagram_link_here.png)  <!-- Placeholder for a future diagram -->

### Stage 1: Geometry Generation (`horn-geometry`)

-   **Responsibility:** Generates a 3D model of the horn based on a set of input parameters (e.g., throat size, mouth size, flare rate, length).
-   **Inputs:** A configuration file (e.g., JSON or YAML) specifying the horn's parameters.
-   **Technology:** FreeCAD (or another CAD kernel) accessed programmatically.
-   **Output:** A standard 3D geometry file, such as a `.stp` (STEP) file.

### Stage 2: Meshing (`horn-solver`)

-   **Responsibility:** Takes the 3D model and divides its volume into a finite element mesh suitable for simulation. This is a critical prerequisite for the solver.
-   **Inputs:** A `.stp` file from the Geometry stage.
-   **Technology:** `gmsh` for mesh generation.
-   **Output:** A mesh object consumable by the physics solver.

### Stage 3: Solving (`horn-solver`)

-   **Responsibility:** Applies the principles of acoustics (using either the Boundary Element Method (BEM) or Finite Element Method (FEM)) to the mesh to calculate the horn's acoustic properties. It simulates how sound waves propagate through the horn.
-   **Inputs:**
    -   The mesh from the Meshing stage.
    -   Parameters for the driver unit (e.g., Thiele/Small parameters like `Bl`, `Re`).
    -   A frequency range to simulate.
-   **Technology:** `FEniCSx/dolfinx` for the FEM/BEM calculations.
-   **Output:** A dataset containing the raw simulation results, likely frequency vs. Sound Pressure Level (SPL) and other acoustic metrics. This will be saved to a `.csv` or similar format.

### Stage 4: Analysis & Visualization (`horn-analysis`)

-   **Responsibility:** Processes the raw data from the Solver to generate human-readable plots and summary statistics. This stage turns raw numbers into actionable insights about the horn's performance.
-   **Inputs:** The results dataset from the Solver stage.
-   **Technology:** `pandas` for data manipulation, `matplotlib` or `plotly` for plotting.
-   **Outputs:**
    -   Frequency response plot (SPL vs. Frequency).
    -   Polar plots (sound dispersion).
    -   Impedance plots.
    -   A summary report of key performance indicators (KPIs).

## 3. Package Overview

-   `horn-core`: A central package for shared data structures, constants, and utilities used across the pipeline (e.g., defining a standard `Driver` or `HornParameters` class).
-   `horn-geometry`: The implementation of the Geometry Generation stage.
-   `horn-solver`: The implementation of both the Meshing and Solving stages, as they are tightly coupled by the `dolfinx` and `gmsh` dependencies.
-   `horn-analysis`: (To be created) The implementation of the Analysis & Visualization stage.

## 4. Orchestration

A master script or orchestrator (`main.py`) will be responsible for executing the pipeline in sequence, passing the output of one stage as the input to the next. For development, this can be a simple Python script; for production, this could be a more robust workflow manager like Airflow or a simple shell script that runs the containers in order. 