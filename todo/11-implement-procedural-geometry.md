# 11: Implement Procedural Horn Geometry

**Date:** 2024-07-26

**Status:** Paused

## 1. Goal

To replace the static `test_box.stp` file with a dynamic, parametric geometry generation system within the `horn-geometry` package. This is the next critical step towards a true end-to-end simulation pipeline.

## 2. Scope of Work

1.  **Create `horn-geometry` Package:** The directory `packages/horn-geometry` exists, but it is empty. It needs to be set up as a proper Python package with an `__init__.py` and a `src` directory.

2.  **Implement a Geometry Generator:**
    -   Create a new module, e.g., `packages/horn-geometry/src/horn_geometry/generator.py`.
    -   Implement a function, e.g., `create_conical_horn(throat_radius, mouth_radius, length)`.
    -   This function will use a CAD kernel (likely the `gmsh` built-in kernel for simplicity, to avoid adding new dependencies like FreeCAD for now) to construct the horn geometry.
    -   The function must export the final geometry as a STEP file.

3.  **Create an End-to-End Test:**
    -   Create a new integration test, e.g., `packages/horn-geometry/tests/test_generator.py`.
    -   This test will call `create_conical_horn`.
    -   It will then pass the generated STEP file path to the `horn-solver`'s `run_simulation_from_step` function.
    -   The test will assert that the simulation runs to completion and generates a results file. This will be our first true end-to-end test of the entire pipeline.

## 3. Key Decisions

-   **CAD Kernel:** We will use `gmsh`'s internal OpenCASCADE kernel (`gmsh.model.occ`) to build the geometry. A detailed analysis in `docs/GEOMETRY_KERNEL_EVALUATION.md` concluded this is the optimal path as it introduces no new dependencies and tightly integrates with our existing meshing workflow.
-   **Initial Profile:** A simple conical horn is the ideal first target. We can add more complex profiles (e.g., exponential) later. 