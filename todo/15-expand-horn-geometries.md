# 15: Expand Horn Geometry Generation

**Date:** 2025-07-23

**Status:** Not Started ⚪

## 1. Goal

To extend the `horn-geometry` package to support the procedural generation of various common horn profiles beyond the current conical horn. This will enable the simulation and analysis of more complex and acoustically diverse horn designs.

## 2. Detailed Steps

### Step 1: Refactor Geometry Generator
- **Action:** Review the existing `horn_geometry/generator.py`.
- **Action:** Refactor the code to create a clear, extensible API. This should involve creating a base `Horn` class with common methods (`generate_mesh`, `export_step`, etc.) and then creating subclasses for each horn type.
- **Deliverable:** A refactored `generator.py` with a class-based structure for defining horn geometries.

### Step 2: Implement Exponential Horn
- **Action:** Create a new subclass, `ExponentialHorn`.
- **Action:** Implement the mathematical formula for the exponential profile: `A(x) = A_t * exp(m*x)`, where `A(x)` is the cross-sectional area at distance `x`, `A_t` is the throat area, and `m` is the flare rate constant.
- **Action:** Add command-line arguments and parameters to the Nextflow pipeline to allow users to select the `exponential` horn type and specify its parameters (e.g., `flare_rate`, `mouth_diameter`, `length`).
- **Deliverable:** A new `ExponentialHorn` class and updated CLI/pipeline parameters.

### Step 3: Implement Tractrix Horn
- **Action:** Create a new subclass, `TractrixHorn`.
- **Action:** Implement the mathematical formula for the Tractrix curve, which defines the horn's profile. The Tractrix curve is defined by the property that the length of the tangent from the point of tangency to the asymptote is constant.
- **Action:** Add the necessary parameters for the `tractrix` horn type to the pipeline.
- **Deliverable:** A new `TractrixHorn` class.

### Step 4: Implement Spherical/Wave-Guide Horn
- **Action:** Create a new subclass, `SphericalHorn`.
- **Action:** This profile is often used in modern compression drivers and is defined by a spherical wavefront. The implementation will involve generating a profile that follows a segment of a circle.
- **Action:** Add the necessary parameters (e.g., `radius`, `angle`) for the `spherical` horn type.
- **Deliverable:** A new `SphericalHorn` class.

### Step 5: (Optional) Implement JMLC / Le Cléac'h Horn
- **Action:** Research the profile developed by Jean-Michel Le Cléac'h. This is a more complex, numerically-derived profile designed for a specific acoustic response.
- **Action:** Implement a generator for this profile, which may involve interpolating from published data or implementing the generating equations.
- **Action:** Add the `jmlc` horn type to the pipeline.
- **Deliverable:** A new `JMLCHorn` class.

### Step 6: Update Testing
- **Action:** For each new horn type, create a dedicated unit test in `packages/horn-geometry/tests/`.
- **Action:** The tests should:
    1. Generate a STEP file for the new geometry.
    2. Verify that the key dimensions (throat area, mouth area, length) are correct.
    3. (Optional) Verify that the generated mesh is valid and has no degenerate elements.
- **Deliverable:** New test files (e.g., `test_exponential_horn.py`) for each geometry.

### Step 7: Documentation
- **Action:** Update the `README.md` and any relevant documentation in `docs/` to describe the new horn types.
- **Action:** Clearly document the required parameters for each geometry.
- **Deliverable:** Updated documentation files.

## 3. Key Challenges
- **Mathematical Correctness:** Ensuring the formulas for each horn profile are implemented accurately.
- **Geometric Robustness:** The geometry kernel (`pygmsh` or `gmsh`) must be able to successfully mesh the more complex curves of the new profiles without errors. This may require careful tuning of meshing parameters.
- **Parameterization:** Defining a clear and intuitive set of parameters for each horn type.
