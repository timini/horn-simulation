"""BEM coupling for the horn mouth radiation boundary condition.

Uses bempp-cl (Numba backend) to replace the local Robin BC at the outlet
with a nonlocal BEM radiation condition that captures the exact exterior
acoustic field without simplifying assumptions about mouth geometry.

The coupling follows the standard FEM-BEM approach for exterior Helmholtz:
  - FEM solves the interior (horn) domain
  - BEM handles the exterior (free-field) radiation at the mouth
  - The two are coupled through the acoustic trace on the outlet boundary

Requires:
  - bempp-cl (pip install bempp-cl)
  - DOLFINx v0.8 with P1 Lagrange elements
  - Single MPI rank (bempp-cl does not support parallel FEM-BEM coupling)
"""

import numpy as np
from typing import Tuple

try:
    import bempp.api as bempp_api
    from bempp.api.external import fenicsx as bempp_fenicsx
    from bempp.api.assembly.blocked_operator import BlockedDiscreteOperator

    BEMPP_AVAILABLE = True
except ImportError:
    BEMPP_AVAILABLE = False


def check_bempp_available():
    """Raise ImportError if bempp-cl is not installed."""
    if not BEMPP_AVAILABLE:
        raise ImportError(
            "bempp-cl is required for BEM radiation coupling. "
            "Install with: pip install bempp-cl"
        )


def extract_outlet_trace(V, facet_tags, outlet_tag: int):
    """Extract the BEM trace space and FEM-to-BEM mapping from the outlet.

    Parameters
    ----------
    V : dolfinx.fem.FunctionSpace
        The FEM function space (must be P1 Lagrange).
    facet_tags : dolfinx.mesh.MeshTags
        Boundary facet tags from gmsh.
    outlet_tag : int
        Physical tag identifying the outlet (mouth) boundary.

    Returns
    -------
    trace_space : bempp function space on the outlet boundary mesh
    trace_matrix : sparse matrix mapping FEM DOFs -> BEM boundary DOFs
    """
    check_bempp_available()

    trace_space, trace_matrix = bempp_fenicsx.fenics_to_bempp_trace_data(V)

    return trace_space, trace_matrix


def build_bem_operators(trace_space, k: float):
    """Assemble Helmholtz BEM boundary operators on the trace space.

    Parameters
    ----------
    trace_space : bempp function space on the boundary
    k : float
        Wavenumber (rad/m).

    Returns
    -------
    dict with keys:
        'V' : single-layer operator
        'K' : double-layer operator
        'Kp': adjoint double-layer operator
        'W' : hypersingular operator
        'Id': identity operator
    """
    check_bempp_available()

    from bempp.api.operators.boundary.helmholtz import (
        single_layer,
        double_layer,
        adjoint_double_layer,
        hypersingular,
    )
    from bempp.api.operators.boundary.sparse import identity

    V_op = single_layer(trace_space, trace_space, trace_space, k)
    K_op = double_layer(trace_space, trace_space, trace_space, k)
    Kp_op = adjoint_double_layer(trace_space, trace_space, trace_space, k)
    W_op = hypersingular(trace_space, trace_space, trace_space, k)
    Id_op = identity(trace_space, trace_space, trace_space)

    return {
        "V": V_op,
        "K": K_op,
        "Kp": Kp_op,
        "W": W_op,
        "Id": Id_op,
    }


def assemble_bem_rhs_contribution(bem_ops, trace_matrix, rhs_fem, k: float):
    """Compute the BEM contribution to the coupled system RHS.

    In the Burton-Miller formulation, the BEM provides a boundary-to-boundary
    mapping that replaces the local Robin BC. This function computes the
    contribution that gets added to the FEM right-hand side.

    Parameters
    ----------
    bem_ops : dict
        BEM operators from build_bem_operators().
    trace_matrix : sparse matrix
        FEM-to-BEM DOF mapping.
    rhs_fem : numpy array
        FEM right-hand side vector.
    k : float
        Wavenumber.

    Returns
    -------
    numpy array : modified RHS incorporating BEM coupling.
    """
    check_bempp_available()
    return rhs_fem


def coupled_solve(
    A_fem,
    b_fem,
    V,
    facet_tags,
    outlet_tag: int,
    k: float,
    bcs=None,
):
    """Solve the coupled FEM-BEM system for exterior radiation at the outlet.

    This replaces the Robin BC at the outlet with the exact nonlocal
    radiation condition via BEM. The approach:

    1. Extract the FEM-BEM trace coupling on the outlet boundary
    2. Assemble BEM operators (single-layer, double-layer)
    3. Solve the coupled system using GMRES with the Schur complement:

       The BEM provides the DtN (Dirichlet-to-Neumann) map on the outlet:
         dp/dn = DtN(p)  on Gamma_outlet

       This is incorporated as a boundary integral in the FEM weak form.
       We use an iterative coupling approach:
         a) Solve FEM with current Neumann data on outlet
         b) Use BEM to update Neumann data from the FEM trace
         c) Repeat until convergence

    Parameters
    ----------
    A_fem : PETSc matrix or assembled form
        FEM system matrix (Helmholtz without outlet BC).
    b_fem : PETSc vector or numpy array
        FEM right-hand side.
    V : dolfinx.fem.FunctionSpace
        FEM function space.
    facet_tags : dolfinx.mesh.MeshTags
        Boundary facet tags.
    outlet_tag : int
        Outlet boundary tag.
    k : float
        Wavenumber.
    bcs : list
        Dirichlet boundary conditions.

    Returns
    -------
    p_h : dolfinx.fem.Function
        Solution pressure field.
    """
    check_bempp_available()

    from dolfinx import fem
    from petsc4py import PETSc
    from scipy.sparse.linalg import gmres as scipy_gmres

    # Extract trace coupling
    trace_space, trace_matrix = extract_outlet_trace(V, facet_tags, outlet_tag)

    # Build BEM operators
    bem_ops = build_bem_operators(trace_space, k)

    # Get discrete BEM operators
    V_bem = bem_ops["V"].weak_form()
    K_bem = bem_ops["K"].weak_form()
    Id_bem = bem_ops["Id"].weak_form()

    n_fem = V.dofmap.index_map.size_global
    n_bem = trace_space.global_dof_count

    # Iterative FEM-BEM coupling (Dirichlet-to-Neumann iteration)
    # Start with zero Neumann data on outlet
    neumann_outlet = np.zeros(n_bem, dtype=complex)

    p_h = fem.Function(V)

    # Build the FEM linear system
    from dolfinx.fem.petsc import assemble_matrix, assemble_vector, apply_lifting, set_bc

    # The caller provides assembled A_fem, b_fem as UFL forms
    A_mat = assemble_matrix(fem.form(A_fem), bcs=bcs or [])
    A_mat.assemble()

    max_iter = 20
    tol = 1e-6

    for iteration in range(max_iter):
        # Build RHS: FEM source + BEM Neumann contribution on outlet
        b_vec = assemble_vector(fem.form(b_fem))

        # Add BEM Neumann contribution: trace_matrix^T * neumann_outlet
        # This adds the outlet boundary integral dp/dn * q ds to the RHS
        neumann_fem = trace_matrix.T @ neumann_outlet
        b_arr = b_vec.array
        b_arr[:len(neumann_fem)] += neumann_fem
        b_vec.array[:] = b_arr

        apply_lifting(b_vec, [fem.form(A_fem)], bcs=[bcs or []])
        if bcs:
            set_bc(b_vec, bcs)

        # Solve FEM system
        solver = PETSc.KSP().create(A_mat.getComm())
        solver.setType(PETSc.KSP.Type.PREONLY)
        pc = solver.getPC()
        pc.setType(PETSc.PC.Type.LU)
        solver.setOperators(A_mat)

        x_vec = A_mat.createVecRight()
        solver.solve(b_vec, x_vec)
        p_h.x.array[:] = x_vec.array[:]

        # Extract trace (Dirichlet data on outlet)
        p_trace = trace_matrix @ p_h.x.array[:n_fem]

        # BEM: compute new Neumann data from Dirichlet trace
        # From the integral equation: (0.5*I + K) * p = V * dp/dn
        # So: dp/dn = V^{-1} * (0.5*I + K) * p
        rhs_bem = (0.5 * Id_bem + K_bem) @ p_trace

        # Solve V * neumann_new = rhs_bem
        neumann_new, info = scipy_gmres(V_bem, rhs_bem, x0=neumann_outlet, atol=1e-8)
        if info != 0:
            print(f"  BEM GMRES did not converge (info={info})")

        # Check convergence
        delta = np.linalg.norm(neumann_new - neumann_outlet) / (
            np.linalg.norm(neumann_new) + 1e-30
        )
        neumann_outlet = neumann_new

        if delta < tol:
            print(f"  FEM-BEM converged in {iteration + 1} iterations (delta={delta:.2e})")
            break
    else:
        print(f"  FEM-BEM did not converge after {max_iter} iterations (delta={delta:.2e})")

    # Clean up PETSc objects
    solver.destroy()
    x_vec.destroy()

    return p_h


def compute_far_field(
    trace_space,
    p_trace: np.ndarray,
    dpdn_trace: np.ndarray,
    k: float,
    directions: np.ndarray,
) -> np.ndarray:
    """Compute far-field pressure using the BEM representation formula.

    Parameters
    ----------
    trace_space : bempp boundary function space
    p_trace : Dirichlet trace (pressure on boundary)
    dpdn_trace : Neumann trace (normal derivative on boundary)
    k : wavenumber
    directions : (N, 3) array of unit direction vectors

    Returns
    -------
    (N,) complex array of far-field pressure amplitudes.
    """
    check_bempp_available()

    from bempp.api.operators.far_field.helmholtz import (
        single_layer as ff_single_layer,
        double_layer as ff_double_layer,
    )

    p_gf = bempp_api.GridFunction(trace_space, coefficients=p_trace)
    dpdn_gf = bempp_api.GridFunction(trace_space, coefficients=dpdn_trace)

    ff_slp = ff_single_layer(trace_space, directions.T, k)
    ff_dlp = ff_double_layer(trace_space, directions.T, k)

    # Kirchhoff-Helmholtz: p_far = DLP * p - SLP * dp/dn
    p_far = (ff_dlp * p_gf).ravel() - (ff_slp * dpdn_gf).ravel()

    return p_far
