# 10: Debugging the Complex-Valued Solver

**Date:** 2024-07-26

**Status:** Completed

## 1. Goal

To systematically debug and resolve the `ValueError: Unexpected complex value in real expression` that is blocking the implementation of the FEM solver. This requires a more fundamental approach than previous attempts.

## 2. Debugging Strategies

We will pursue the following strategies in order.

### Strategy 1: Isolate with a Minimal Reproducible Example

-   **Action:** Create a new, single, self-contained Python script (`debug_helmholtz.py`).
-   **Details:** This script will not use `gmsh` or any of our project's modules. It will use `dolfinx.mesh.create_unit_square` to generate a mesh and solve the simplest possible Helmholtz problem.
-   **Goal:** To determine if the `dolfinx v0.8.0` API for complex problems works as expected in the most basic case, completely isolated from our project's code.

### Strategy 2: Deconstruct the Canonical Demo

-   **Action:** Copy the official `dolfinx` Helmholtz demo source code into a local script (`debug_demo.py`).
-   **Details:** Run the exact, unmodified demo code inside our `dolfinx/dolfinx:v0.8.0` container.
-   **Goal:**
    -   If it fails, it indicates a potential problem with the Docker image's configuration.
    -   If it succeeds, incrementally modify the working demo code to resemble our `solver.py` implementation. The specific change that introduces the failure will reveal the root cause.

### Strategy 3: Manual Decomposition into Real and Imaginary Parts

-   **Action:** If native complex support remains elusive, rewrite the physics.
-   **Details:** The complex Helmholtz equation `∇²p + k²p = 0` can be split into a system of two real equations for `p_real` and `p_imag`. This involves creating a mixed `FunctionSpace` and solving a `2x2` system.
-   **Goal:** This is a workaround to unblock the project. It avoids the problematic JIT path for complex numbers entirely by solving a more complex, purely real-valued system. This is a last resort.

## 3. Current Action

Proceeding with **Strategy 2: Deconstruct the Canonical Demo**, as it provides the most direct comparison against a known-good implementation.

## 4. Resolution

The issue was successfully resolved. The root cause was not an error in the Python code's use of the UFL API, but a fundamental misunderstanding of the Docker environment.

-   **Discovery:** Running a canonical demo script inside the container revealed the true issue: the `dolfinx/dolfinx:v0.8.0` image, while containing a complex-enabled build of PETSc, **defaults to a real-only environment**.

-   **Initial Flawed Fix:** The initial attempt to fix this was to `RUN . /usr/local/bin/dolfinx-complex-mode` in the Dockerfile. This failed because each `RUN` command executes in a separate shell, and the environment variables set by the script did not persist.

-   **The Correct Solution:** The `dolfinx-complex-mode` script was inspected to identify the five environment variables it sets (`PKG_CONFIG_PATH`, `CMAKE_PREFIX_PATH`, `PETSC_ARCH`, `PYTHONPATH`, `LD_LIBRARY_PATH`). These were then permanently set in the `Dockerfile` using the `ENV` instruction.

This ensured that all subsequent build steps and the final container runtime were executed within the correct complex-valued environment, which immediately resolved the persistent `ValueError`. Subsequent minor API mismatches were then easily debugged and fixed. 