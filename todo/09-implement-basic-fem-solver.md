# 09: Implement a Basic FEM/BEM Solver

**Date:** 2024-07-26

**Status:** Completed

## 1. Goal

Replace the dummy placeholder `run_simulation` function in `packages/horn-solver/src/horn_solver/solver.py` with a basic, but physically meaningful, acoustic simulation using the Finite Element Method (FEM) with `dolfinx`.

This is the most critical task to achieve a meaningful end-to-end simulation pipeline.

## 2. Scope of Work

1.  **Physics Modeling:**
    -   Implement a solver for the time-harmonic acoustic wave equation (the Helmholtz equation).
    -   Define appropriate boundary conditions:
        -   An **inlet** surface (the "driver") with a specified normal velocity.
        -   A **radiation** or **absorbing** boundary condition at the outlet (the "mouth") to simulate sound radiating into open space.
        -   **Hard walls** for all other surfaces (zero normal velocity).

2.  **Gmsh Integration:**
    -   Update the meshing logic in `run_simulation_from_step` to tag the boundary surfaces (inlet, outlet, walls) with physical markers. This is essential for applying the correct boundary conditions in the solver.

3.  **Solver Implementation:**
    -   Define the FEM variational problem (weak form) for the Helmholtz equation.
    -   Write the `dolfinx` code to assemble and solve the resulting linear system for the acoustic pressure.
    -   The solver must handle complex numbers, as acoustic pressure is a complex quantity.

4.  **Frequency Sweep:**
    -   Loop through the input `freq_range` to solve the problem at multiple frequencies.

5.  **Output Generation:**
    -   Calculate a meaningful output metric, such as the average sound pressure at the horn's mouth for each frequency.
    -   Save the results (`frequency`, `spl`) to the output CSV file, maintaining the existing format.

## 3. Testing Strategy

A detailed `TESTING_STRATEGY.md` document will be created. The core approach will be **verification** against a known analytical solution (e.g., a plane wave in a straight tube) to confirm the solver's correctness before applying it to complex horn geometries.

## 4. Implementation Summary & Challenges

This task was successfully completed, but involved overcoming a major, unexpected hurdle related to the solver's environment.

-   **Boundary Tagging:** The meshing logic in `solver.py` was successfully refactored and extended. A TDD approach was used to implement physical group tagging for the "inlet," "outlet," and "wall" surfaces based on their coordinates, which was verified with a new test case.
-   **FEM Implementation:** The dummy `run_simulation` function was replaced with a `dolfinx` implementation of the time-harmonic Helmholtz equation, using the tagged boundaries to apply a source term and a radiation condition.
-   **The JIT Compilation Failure:** The primary challenge was a persistent `ValueError: Unexpected complex value in real expression`. This blocked progress for a significant time and was initially misdiagnosed as an API usage error.
-   **Resolution:** As detailed in `10-debugging-complex-solver.md`, the root cause was discovered to be an incorrect environment configuration in the Docker image. The issue was resolved by permanently setting the required environment variables for the complex-petsc build using the `ENV` instruction in the `Dockerfile`. After fixing the environment, remaining `dolfinx` API issues were minor and quickly resolved.

The final, working solver is now committed and verified by the test suite. 