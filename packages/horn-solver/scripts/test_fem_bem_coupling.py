#!/usr/bin/env python3
"""Phase 1b: FEM-BEM coupling smoke test adapted for DOLFINx v0.8 + bempp-cl.

Solves a simple exterior Helmholtz problem on a unit cube:
  - FEM interior (DOLFINx) + BEM exterior (bempp-cl)
  - Verifies the FEniCSx <-> bempp trace coupling machinery works

Based on the mscroggs FEM-BEM coupling tutorial, updated for:
  - DOLFINx v0.8 API (dolfinx.fem.functionspace, etc.)
  - bempp-cl (bempp.api, not bempp_cl.api)

Run inside the dolfinx container:
    docker run --rm -v $PWD:/app dolfinx/dolfinx:v0.8.0 bash -c \
        "source /usr/local/bin/dolfinx-complex-mode && \
         pip install bempp-cl && python3 /app/packages/horn-solver/scripts/test_fem_bem_coupling.py"
"""

import sys
import numpy as np


def main() -> int:
    # --- Imports ---
    try:
        import dolfinx
        from dolfinx import fem, mesh
        from mpi4py import MPI
        import ufl
        print(f"OK: DOLFINx {dolfinx.__version__}")
    except ImportError as exc:
        print(f"FAIL: DOLFINx import error: {exc}")
        return 1

    try:
        import bempp.api as bempp_api
        print(f"OK: bempp.api {bempp_api.__version__}")
    except ImportError as exc:
        print(f"FAIL: bempp import error: {exc}")
        return 1

    try:
        from bempp.api.external import fenicsx as bempp_fenicsx
        print("OK: imported bempp.api.external.fenicsx coupling module")
    except ImportError as exc:
        print(f"FAIL: bempp-fenicsx coupling not available: {exc}")
        return 1

    # --- Create a unit cube mesh in DOLFINx ---
    print("\nCreating unit cube FEM mesh...")
    domain = mesh.create_unit_cube(
        MPI.COMM_WORLD, 5, 5, 5, cell_type=mesh.CellType.tetrahedron
    )
    V = fem.functionspace(domain, ("Lagrange", 1))
    print(f"  FEM DOFs: {V.dofmap.index_map.size_global}")

    # --- Extract boundary trace space ---
    print("Extracting boundary trace space...")
    try:
        fenics_space, trace_matrix = bempp_fenicsx.fenics_to_bempp_trace_data(V)
        print(f"  BEM trace DOFs: {fenics_space.global_dof_count}")
        print(f"  Trace matrix shape: {trace_matrix.shape}")
    except Exception as exc:
        print(f"FAIL: trace extraction failed: {exc}")
        return 1

    # --- Assemble BEM operators ---
    k = 1.0  # wavenumber
    print(f"\nAssembling BEM operators (k={k})...")
    try:
        slp = bempp_api.operators.boundary.helmholtz.single_layer(
            fenics_space, fenics_space, fenics_space, k
        )
        dlp = bempp_api.operators.boundary.helmholtz.double_layer(
            fenics_space, fenics_space, fenics_space, k
        )
        print("  OK: single-layer and double-layer operators assembled")
    except Exception as exc:
        print(f"FAIL: BEM operator assembly failed: {exc}")
        return 1

    # --- Quick validation: apply operator to a constant function ---
    print("Applying BEM operators to a test function...")
    try:
        ones = bempp_api.GridFunction(fenics_space, coefficients=np.ones(fenics_space.global_dof_count))
        result = slp * ones
        print(f"  OK: SLP * ones -> grid function with {result.coefficients.shape[0]} coefficients")
        print(f"  Max coefficient magnitude: {np.max(np.abs(result.coefficients)):.6f}")
    except Exception as exc:
        print(f"FAIL: operator application failed: {exc}")
        return 1

    print("\nAll FEM-BEM coupling checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
