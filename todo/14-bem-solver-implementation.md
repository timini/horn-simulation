# 14: BEM Solver Implementation

**Date:** 2025-07-23

**Status:** Not Started ⚪

## 1. Goal

Implement a Boundary Element Method (BEM) solver to accurately simulate the acoustic radiation from the horn mouth into an open, infinite domain. This approach is more computationally efficient and physically accurate for open-domain radiation problems than a pure FEM model with absorbing boundary conditions.

This task will focus on creating a standalone BEM solver that can be integrated into the existing Nextflow pipeline.

## 2. Detailed Steps

### Step 1: BEM Library and Environment Setup
- **Action:** Finalize the choice of a Python BEM library. `Bempp` is the primary candidate due to its features and potential for FEniCSx integration.
- **Action:** Create a new, dedicated Docker environment for the BEM solver to manage its specific dependencies. This will live in a new `horn-bem-solver` package.
- **Deliverable:** A `Dockerfile` that successfully builds with `bempp` and its dependencies installed.

### Step 2: Basic BEM Solver Structure
- **Action:** Create a new package `packages/horn-bem-solver`.
- **Action:** Inside the new package, create the main solver script, e.g., `bem_solver.py`.
- **Action:** The script will be responsible for:
    - Loading a surface mesh of the horn mouth.
    - Defining physical parameters (frequency, speed of sound, air density).
    - Setting up the Bempp function spaces (`DP`, `P1`, etc.).
- **Deliverable:** A runnable script that can load a mesh and initialize the BEM spaces.

### Step 3: Define Boundary Conditions
- **Action:** Implement the velocity boundary condition on the horn mouth surface. This will be a uniform velocity for initial tests, representing the output from the horn's throat.
- **Action:** The Sommerfeld radiation condition at infinity is implicitly handled by the BEM formulation, so no explicit implementation is needed, but this should be noted in the code documentation.
- **Deliverable:** A function within `bem_solver.py` that applies the velocity boundary condition to the BEM problem.

### Step 4: Helmholtz Equation Formulation
- **Action:** Formulate the Helmholtz integral equation using the chosen BEM library.
- **Action:** Assemble the boundary operators (single-layer, double-layer, and hypersingular operators).
- **Action:** Solve the resulting system of linear equations to find the pressure on the boundary.
- **Deliverable:** A working solver that computes the acoustic pressure on the horn mouth surface for a given frequency.

### Step 5: Post-Processing and Far-Field Calculation
- **Action:** Once the surface pressure is known, calculate the sound pressure at points in the far-field.
- **Action:** Implement a function to define a listening grid (e.g., a semi-sphere in front of the horn mouth).
- **Action:** Use the Kirchhoff-Helmholtz integral formula to evaluate the pressure at the far-field points.
- **Action:** Calculate and plot the on-axis and off-axis Sound Pressure Level (SPL).
- **Deliverable:** A script that can generate SPL plots from the BEM solution.

### Step 6: Integration into Nextflow Pipeline
- **Action:** Create a new Nextflow process in `main.nf` for the BEM solver.
- **Action:** This process will take the horn mouth mesh from the geometry stage as input.
- **Action:** The process will execute the `bem_solver.py` script within its Docker container.
- **Action:** The results (SPL data, plots) will be published to the `results` directory.
- **Deliverable:** An updated `main.nf` and `nextflow.config` that incorporates the BEM solver.

### Step 7: Verification and Validation
- **Action:** Create a verification test case by simulating a simple pulsating sphere and comparing the BEM results to the known analytical solution.
- **Action:** Compare the BEM results for the conical horn with the results from the FEM solver to understand the differences and validate the implementation.
- **Deliverable:** Test scripts and documentation of the verification results.

## 3. Key Challenges
- **Dependency Management:** Ensuring a stable and reproducible Docker environment for `bempp`.
- **Performance:** BEM can be memory-intensive. The solver will need to be optimized for larger meshes and higher frequencies.
- **Accuracy:** Ensuring the numerical implementation is correct and provides physically meaningful results.
