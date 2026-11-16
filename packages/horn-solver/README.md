# Horn Solver

This package contains the core simulation logic for the horn loudspeaker project.

## Architecture

- **Interior:** `dolfinx` solves the Helmholtz equation for the acoustic pressure inside the horn's air volume (`solver.py`).
- **Mouth radiation:** `modal_baffled` (default for generated axisymmetric horns) couples 16 axisymmetric aperture modes to a Rayleigh radiation impedance in an infinite baffle (`modal_boundary.py`, maths in `horn_core.modal_radiation`). The local models `flanged_piston`, `unflanged_piston`, `finite_flange`, `plane_wave` and `closed` are Robin or rigid boundaries.
- **Meshing:** `meshing.py` (gmsh only, no dolfinx) sets the degree-dependent element size, the throat refinement field, the 1 % throat-area guard and the node limit, and turns an unmeshable candidate into a per-candidate rejection.

The earlier FEM-BEM coupling (`bem_coupling.py`, bempp-cl) mapped the whole boundary instead of the mouth and is not reachable from the solver; only its standalone operator tests remain. See `docs/EXTERIOR_RADIATION_VALIDATION.md` and `docs/MODAL_APERTURE.md`.
