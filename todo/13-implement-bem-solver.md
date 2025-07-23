# 13: Implement BEM for Acoustic Radiation

**Date:** 2025-07-22

**Status:** Not Started ⚪

## 1. Goal

To significantly improve the accuracy of the simulation by replacing the current, simplified radiation boundary condition (a Robin boundary condition) with a more physically accurate model using the Boundary Element Method (BEM). This involves coupling the existing `dolfinx` FEM solver with a BEM library like `Bempp`.

## 2. Background

The current FEM model simulates the acoustic waves inside the horn volume. However, to model how the sound radiates from the horn's mouth into open space, we need a method that can handle infinite domains. BEM is ideal for this. We will use a coupled FEM-BEM approach:

1.  **FEM (Interior):** Solves the Helmholtz equation inside the horn, as is currently done.
2.  **BEM (Exterior):** Solves for the sound pressure on the radiating surface (the horn's mouth) and the surrounding space.

The two solvers are coupled at the boundary between them (the "mouth" surface).

## 3. Scope of Work

1.  **Dependency and Environment Setup:**
    -   **Status:** Completed (with limitations) ✅
    -   **Summary:** The initial attempt to create a Docker environment on an ARM64 host failed due to the lack of a pre-built `bempp-cl` package for the `linux/arm64` architecture. A detailed analysis confirmed this is a known gap in the ecosystem.
    -   **Resolution:** A two-pronged approach was attempted:
        1.  **Short-Term Workaround:** A `Dockerfile.amd64` was created to build a `linux/amd64` image using pre-built `conda-forge` packages. This build was successful and provides a functional, albeit emulated, environment for development.
        2.  **Long-Term Solution:** A `Dockerfile` for native `linux/arm64` builds was created. However, this build failed because a critical dependency, `pocl-dev`, is not available in the `ubuntu:22.04` arm64 repositories. This makes a native build infeasible at this time.
    -   **Conclusion:** The `amd64` emulated environment is the only viable solution for now. All development will proceed using the `Dockerfile.amd64` image.
    -   **Next Steps:** With the environment now established, the next step is to implement the BEM solver itself.

2.  **BEM Solver Implementation:**
    -   Create a new Python module, e.g., `horn_solver/bem_solver.py`.
    -   Define the function spaces for the BEM problem on the boundary mesh (the horn mouth).
    -   Assemble the BEM operators (single-layer and double-layer potentials).
    -   Create a function that takes the FEM solution at the boundary as input and returns the BEM solution.

3.  **FEM-BEM Coupling:**
    -   Modify the main `run_simulation` function in `horn_solver/solver.py`.
    -   The variational form of the FEM problem will be modified to include the BEM operator, which acts as a non-local boundary condition at the horn mouth.
    -   The combined FEM-BEM system will be solved as a single, monolithic problem.

4.  **Testing Strategy:**
    -   **Verification:** Create a new test case to verify the coupled solver against a known analytical solution, such as the pulsating sphere. This is essential to ensure the coupling is implemented correctly.
    -   **Integration:** Update the existing end-to-end tests to use the new coupled solver.

## 4. Key Challenges

-   **Environment:** `Bempp` and `dolfinx` have complex dependencies. Creating a stable Docker environment where they coexist is the biggest anticipated risk.
-   **Mathematical Formulation:** The coupling formulation requires a solid understanding of the underlying physics and mathematics of both FEM and BEM.
-   **API Integration:** Ensuring the data structures (meshes, function spaces) from `dolfinx` and `Bempp` are compatible and correctly mapped at the coupling interface.
