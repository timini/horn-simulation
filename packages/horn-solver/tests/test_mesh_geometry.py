"""Independent checks of the basis and geometry used by the numerical backend."""
import numpy as np
import pytest
import basix
import ufl
from dolfinx import fem
from petsc4py import PETSc
from horn_geometry.generator import create_conical_horn
from horn_solver.solver import create_mesh_from_step, run_simulation, OUTLET_TAG


@pytest.mark.parametrize('degree', [1, 2])
def test_lagrange_basis_interpolates_its_nodes_and_preserves_constants(degree):
    element = basix.create_element(basix.ElementFamily.P, basix.CellType.tetrahedron,
                                   degree, basix.LagrangeVariant.equispaced)
    table = element.tabulate(1, element.points)
    np.testing.assert_allclose(table[0, :, :, 0], np.eye(len(element.points)), atol=1e-12)
    np.testing.assert_allclose(table[1:, :, :, 0].sum(axis=2), 0, atol=1e-12)


def test_mixed_geometry_orders_preserve_small_pipe_area(tmp_path):
    step = tmp_path/'pipe.step'
    create_conical_horn(.007, .007, .18, step)
    expected = np.pi * .007**2
    for degree in (1, 2, 1, 2):
        domain, tags = create_mesh_from_step(str(step), .006, .18, geometry_order=degree)
        ds = ufl.Measure('ds', domain=domain, subdomain_data=tags)
        one = fem.Constant(domain, PETSc.ScalarType(1))
        measured = fem.assemble_scalar(fem.form(one * ds(OUTLET_TAG))).real
        np.testing.assert_allclose(measured, expected, rtol=.11 if degree == 1 else .002)


def test_solver_rejects_gross_mesh_cad_area_disagreement(tmp_path):
    step = tmp_path/'pipe.step'
    create_conical_horn(.007, .007, .18, step)
    domain, tags = create_mesh_from_step(str(step), .006, .18, geometry_order=2)
    domain.horn_boundary_areas['mouth'] *= 1000
    with pytest.raises(RuntimeError, match='inconsistent with CAD'):
        run_simulation(domain, tags, (200, 400), 2, {}, str(tmp_path/'bad.csv'), element_degree=2)
    assert not (tmp_path/'bad.csv').exists()
