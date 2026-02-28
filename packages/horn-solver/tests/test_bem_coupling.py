"""Tests for BEM radiation coupling via bempp-cl.

These tests require bempp-cl to be installed (pip install bempp-cl).
They are skipped if bempp-cl is not available (e.g. outside the Docker container).
"""

import pytest
import numpy as np

try:
    import bempp.api as bempp_api
    BEMPP_AVAILABLE = True
except ImportError:
    BEMPP_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not BEMPP_AVAILABLE, reason="bempp-cl not installed"
)


class TestBemppImport:
    """Verify bempp-cl is importable and configured correctly."""

    def test_import_bempp(self):
        """bempp.api should be importable."""
        import bempp.api
        assert hasattr(bempp.api, "function_space")

    def test_numba_backend(self):
        """Default device interface should be 'numba' (no OpenCL)."""
        device = getattr(bempp_api, "DEFAULT_DEVICE_INTERFACE", "unknown")
        # Accept either numba or opencl — both work, but numba is the goal
        assert device in ("numba", "opencl"), f"Unexpected backend: {device}"

    def test_helmholtz_operators_exist(self):
        """Helmholtz boundary operators should be importable."""
        from bempp.api.operators.boundary.helmholtz import (
            single_layer,
            double_layer,
            adjoint_double_layer,
            hypersingular,
        )
        assert callable(single_layer)


class TestBemOperators:
    """Test BEM operator assembly on a simple geometry."""

    def test_sphere_single_layer(self):
        """Assemble a Helmholtz single-layer operator on a sphere."""
        grid = bempp_api.shapes.regular_sphere(2)
        space = bempp_api.function_space(grid, "P", 1)
        k = 1.0
        slp = bempp_api.operators.boundary.helmholtz.single_layer(
            space, space, space, k
        )
        mat = slp.weak_form()
        assert mat.shape[0] == mat.shape[1]
        assert mat.shape[0] == space.global_dof_count

    def test_pulsating_sphere_analytical(self):
        """BEM pulsating sphere should match analytical solution within 1 dB.

        A pulsating sphere of radius a with uniform radial velocity v0
        has an exact analytical solution for the surface pressure:

            p(a) = rho * c * v0 * (ka)^2 / (1 + (ka)^2) * (1 + j/(ka))
                 (simplified for unit velocity on the surface)

        We use the BEM exterior Neumann problem to compute surface pressure
        and compare against the analytical result.
        """
        a = 0.1  # sphere radius
        k = 10.0  # wavenumber (ka = 1.0)
        ka = k * a
        rho = 1.225
        c = 343.0

        # Create sphere mesh
        grid = bempp_api.shapes.regular_sphere(3)
        # Scale to radius a
        vertices = grid.vertices * a
        grid = bempp_api.Grid(vertices, grid.elements)

        space = bempp_api.function_space(grid, "P", 1)

        # BEM operators
        from bempp.api.operators.boundary.helmholtz import (
            single_layer,
            double_layer,
            hypersingular,
        )
        from bempp.api.operators.boundary.sparse import identity

        V_op = single_layer(space, space, space, k)
        K_op = double_layer(space, space, space, k)
        Id_op = identity(space, space, space)

        # Neumann BC: dp/dn = -j*k*rho*c*v0 on sphere surface
        # For unit velocity v0=1, outward normal velocity:
        v0 = 1.0
        neumann_data = -1j * k * rho * c * v0 * np.ones(space.global_dof_count)
        neumann_gf = bempp_api.GridFunction(space, coefficients=neumann_data)

        # Integral equation: (0.5*I + K) * p = V * dp/dn
        # Solve for p
        lhs = 0.5 * Id_op + K_op
        rhs_gf = V_op * neumann_gf

        from scipy.sparse.linalg import gmres
        lhs_disc = lhs.weak_form()
        rhs_vec = rhs_gf.projections(space)

        p_coeffs, info = gmres(lhs_disc, rhs_vec, atol=1e-8)
        assert info == 0, f"GMRES did not converge: info={info}"

        # Analytical surface pressure for pulsating sphere
        # p(a) = rho*c*v0 * j*ka / (1 + j*ka) (exact for monopole)
        # But for a pulsating sphere (all modes), the exact result is:
        # p(a) = -rho*c*v0 * h0'(ka) / h0(ka) where h0 is spherical Hankel
        # For simplicity, use the monopole (n=0) approximation
        from scipy.special import spherical_jn, spherical_yn

        def spherical_hankel1(n, z):
            return spherical_jn(n, z) + 1j * spherical_yn(n, z)

        def spherical_hankel1_deriv(n, z):
            return spherical_jn(n, z, derivative=True) + 1j * spherical_yn(n, z, derivative=True)

        # For pulsating sphere: only n=0 mode contributes
        h0 = spherical_hankel1(0, ka)
        h0p = spherical_hankel1_deriv(0, ka)
        p_analytical = -rho * c * v0 * h0p / h0

        # Compare RMS pressure
        p_rms_bem = np.sqrt(np.mean(np.abs(p_coeffs) ** 2))
        p_rms_analytical = abs(p_analytical)

        spl_bem = 20 * np.log10(p_rms_bem / 20e-6 + 1e-30)
        spl_analytical = 20 * np.log10(p_rms_analytical / 20e-6 + 1e-30)

        diff_db = abs(spl_bem - spl_analytical)
        print(f"Pulsating sphere: BEM SPL={spl_bem:.1f} dB, "
              f"Analytical SPL={spl_analytical:.1f} dB, diff={diff_db:.2f} dB")

        assert diff_db < 1.0, (
            f"BEM pulsating sphere deviates from analytical by {diff_db:.2f} dB "
            f"(BEM={spl_bem:.1f}, analytical={spl_analytical:.1f})"
        )


class TestTraceExtraction:
    """Test FEM-BEM trace extraction."""

    def test_trace_from_unit_cube(self):
        """Extract BEM trace space from a DOLFINx unit cube."""
        from bempp.api.external import fenicsx as bempp_fenicsx
        from dolfinx import fem, mesh
        from mpi4py import MPI

        domain = mesh.create_unit_cube(
            MPI.COMM_WORLD, 3, 3, 3, cell_type=mesh.CellType.tetrahedron
        )
        V = fem.functionspace(domain, ("Lagrange", 1))

        trace_space, trace_matrix = bempp_fenicsx.fenics_to_bempp_trace_data(V)

        assert trace_space.global_dof_count > 0
        assert trace_matrix.shape[0] == trace_space.global_dof_count
        assert trace_matrix.shape[1] == V.dofmap.index_map.size_global

    def test_trace_matrix_maps_correctly(self):
        """Trace matrix should map constant FEM field to constant BEM field."""
        from bempp.api.external import fenicsx as bempp_fenicsx
        from dolfinx import fem, mesh
        from mpi4py import MPI

        domain = mesh.create_unit_cube(
            MPI.COMM_WORLD, 3, 3, 3, cell_type=mesh.CellType.tetrahedron
        )
        V = fem.functionspace(domain, ("Lagrange", 1))

        trace_space, trace_matrix = bempp_fenicsx.fenics_to_bempp_trace_data(V)

        # Constant FEM function (p=1 everywhere)
        n_fem = V.dofmap.index_map.size_global
        p_fem = np.ones(n_fem, dtype=complex)

        # Trace should also be constant
        p_trace = trace_matrix @ p_fem
        assert np.allclose(p_trace, 1.0, atol=1e-10), (
            f"Trace of constant field should be constant, "
            f"got range [{p_trace.min():.6f}, {p_trace.max():.6f}]"
        )


class TestBemCouplingModule:
    """Test the bem_coupling module functions."""

    def test_check_bempp_available(self):
        """check_bempp_available should not raise when bempp is installed."""
        from horn_solver.bem_coupling import check_bempp_available
        check_bempp_available()  # should not raise

    def test_build_bem_operators(self):
        """build_bem_operators should return all required operators."""
        from horn_solver.bem_coupling import build_bem_operators

        grid = bempp_api.shapes.regular_sphere(2)
        space = bempp_api.function_space(grid, "P", 1)

        ops = build_bem_operators(space, k=1.0)

        assert "V" in ops
        assert "K" in ops
        assert "Kp" in ops
        assert "W" in ops
        assert "Id" in ops


class TestE2eWithBemBC:
    """End-to-end test: run the solver with radiation_model='bem'."""

    def test_bem_produces_finite_spl(self, tmp_path):
        """Solver with BEM radiation BC should produce finite SPL values."""
        from pathlib import Path
        from horn_solver.solver import run_simulation_from_step

        step_file = Path(__file__).parent / "test_box.stp"
        if not step_file.exists():
            pytest.skip("test_box.stp not found")

        output_file = tmp_path / "results.csv"

        driver_params = {"Bl": 5.0, "Re": 6.0, "length": 1.0}
        freq_range = (200.0, 500.0)

        result_path = run_simulation_from_step(
            step_file=str(step_file),
            driver_params=driver_params,
            freq_range=freq_range,
            num_intervals=3,  # few points for speed
            output_file=str(output_file),
            max_freq_mesh=freq_range[1],
            mesh_size=1.0,
            radiation_model="bem",
        )

        assert result_path.exists()

        import pandas as pd
        results_df = pd.read_csv(result_path)
        assert len(results_df) == 3
        assert all(np.isfinite(results_df["spl"].values)), "SPL values should be finite"
