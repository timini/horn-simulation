# 12: Implement Frequency Sweep and Plotting

**Date:** 2024-07-26

**Status:** In Progress

## 1. Goal

To achieve a minimal but complete end-to-end simulation run that produces a visual output. This involves modifying the solver to run over a range of frequencies and creating a script to plot the results. This is the highest priority for the MVP.

## 2. Scope of Work

1.  **Implement Frequency Sweep in Solver:**
    -   Modify the `run_simulation` function in `packages/horn-solver/src/horn_solver/solver.py`.
    -   The function will iterate over a number of steps within the provided `freq_range`.
    -   For each frequency, it will solve the FEM problem and calculate the resulting SPL.
    -   All results (`frequency`, `spl`) will be collected and saved to a single CSV file.

2.  **Create Plotting Script:**
    -   Create a new script: `scripts/plot_results.py`.
    -   This script will take two command-line arguments: an input CSV file path and an output image file path.
    -   It will use `pandas` to read the data and `matplotlib` to generate a line plot of SPL vs. Frequency.
    -   Add `matplotlib` as a test dependency.

3.  **Update End-to-End Test:**
    -   Modify `test_e2e_meshing_and_solving` in `packages/horn-solver/tests/test_solver.py`.
    -   After the simulation runs, the test will read the output CSV and assert that it contains multiple data rows.
    -   The test will then execute the `plot_results.py` script using `subprocess`.
    -   Finally, it will assert that the plot image file was created and is not empty.

## 3. Key Decisions

-   **Plotting Library:** `matplotlib` is the standard, most well-supported choice for this task and is already available in the base `dolfinx` image.
-   **Workflow:** Keeping the plotting separate from the solver is a good design choice. The solver's job is to produce data; other tools are responsible for visualizing it. 