# Horn Solver Testing Strategy

## 1. Introduction

The results of a physics simulation are meaningless without a rigorous process to ensure their correctness. This document outlines the strategy for testing the `horn-solver` to build confidence in its outputs. Our approach distinguishes between two key concepts:

-   **Verification:** Are we solving the mathematical equations correctly? (i.e., is the code free of bugs and does it correctly implement the chosen numerical model?).
-   **Validation:** Are we solving the correct equations? (i.e., does our mathematical model accurately represent real-world acoustics?).

For the initial MVP, our focus is entirely on **Verification**. We must first prove that our solver correctly implements the physics before we can attempt to validate it against real-world measurements.

## 2. Verification Strategy: Method of Manufactured Solutions

The gold standard for verification is to test the solver against a problem with a known, exact analytical solution. For our case, we will use the **Method of Manufactured Solutions (MMS)**, for which we will solve the acoustic problem for a simple geometry where the outcome is known.

### The Test Case: Plane Wave in a Duct

-   **Geometry:** A simple, straight, rectangular duct. We will use the existing `test_box.stp` for this, treating it as a section of a duct.
-   **Physics:** We will simulate a simple plane wave traveling down the length of theduct.
-   **Analytical Solution:** For a plane wave in a duct, the acoustic pressure `p` at any point `x` along the duct is given by a simple formula. If we impose a velocity at the inlet (x=0), the pressure is known everywhere.
-   **Test Data:** No external test data is required. The test itself will contain:
    1.  The simple geometry (`test_box.stp`).
    2.  The analytical formula for the plane wave.
    3.  The code to run the `dolfinx` solver on the geometry.
    4.  The code to compare the `dolfinx` result to the result from the analytical formula.

### How the Test Will Work

The `test_solver.py` will be updated to perform these steps:

1.  **Run Simulation:** Execute the `run_simulation_from_step` on the `test_box.stp` file for a single, known frequency.
2.  **Extract Numerical Result:** Read the output `.csv` file produced by the solver.
3.  **Calculate Analytical Result:** Using the known physical parameters (frequency, air density, speed of sound) and the plane wave formula, calculate the exact expected pressure at the measurement point.
4.  **Compare and Assert:** Compare the numerical result from the solver with the analytical result. They should match within a small tolerance (to account for numerical discretization error).

Passing this test will give us high confidence that the solver code is correctly implementing the physics of the Helmholtz equation.

## 3. Validation Suite

The validation suite lives in `tests/validation/` and proves the solver produces physically correct results at three tiers of confidence:

| Case | Tier | Tolerance | Reference |
|------|------|-----------|-----------|
| V1: Straight tube | Tier 1 (exact) | 0.5 dB | Analytical forward-traveling wave |
| V2: Conical horn | Tier 1 (approximate) | 3 dB | Webster horn equation |
| V3: Community horn | Tier 2 (cross-validation) | 6 dB | Published measurements |

Run with: `pytest packages/horn-solver/tests/validation/ -v -m validation`

See `tests/validation/README.md` for full documentation of each case, tolerance rationale, and instructions for adding new reference data.

## 4. Validation Strategy (Future Work)

Once the solver is **verified**, the next stage is **validation** against real-world measurements. This involves:

1.  **Physical Measurements:** Measure a real driver and horn combination with a calibrated microphone.
2.  **Model Physical Horn:** Accurately model the exact geometry of the measured horn.
3.  **Model Physical Driver:** Use the measured Thiele/Small parameters of the driver as input to the simulation.
4.  **Compare Results:** Compare the simulated SPL curve to the measured SPL curve. Discrepancies will inform refinements to the physics model (e.g., adding air viscosity, more complex boundary conditions).

V3 in the validation suite provides the infrastructure for this: add digitized measurement data as CSV, and the test framework handles the comparison automatically.

By following this two-stage process—verification first, then validation—we can build a simulation tool that is not only functional but also trustworthy.