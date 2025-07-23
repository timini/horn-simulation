# 03: Simulation

**Status:** Completed ✅

This task covered the implementation of the core computational stages of the pipeline: meshing and solving.

## Key Accomplishments

- **Meshing**: Implemented a `create_mesh_from_step` function in `horn_solver/solver.py` that uses `gmsh` to generate a mesh from a STEP file and tag boundaries for simulation.
- **Solver**: Implemented a `run_simulation` function that uses `dolfinx` to solve the Helmholtz equation for acoustic pressure.
- **E2E Integration**: Integrated the simulation stage into the main Nextflow pipeline.