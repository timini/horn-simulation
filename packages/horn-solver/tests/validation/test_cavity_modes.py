"""Verify the production volume operator independently of horn/TMM references."""
import itertools
import pytest
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
from mpi4py import MPI
from dolfinx import fem, mesh
import ufl
from horn_solver.solver import bulk_helmholtz_form


pytestmark = pytest.mark.validation

def matrix(form):
    assembled = fem.petsc.assemble_matrix(fem.form(form))
    assembled.assemble()
    indptr, indices, values = assembled.getValuesCSR()
    assert np.max(np.abs(values.imag)) < 1e-14
    result = csr_matrix((values.real.copy(), indices.copy(), indptr.copy()), shape=assembled.getSize())
    assembled.destroy()
    return result


def test_rectangular_cavity_modes_converge():
    """First twelve nonzero rigid-box modes: fine P2 error <1%, decreasing."""
    assert MPI.COMM_WORLD.size == 1, 'Modal verification requires a single rank'
    lengths = np.array([.47, .33, .22])
    expected = sorted(343./2*np.linalg.norm(np.array(n)/lengths)
                      for n in itertools.product(range(6), repeat=3) if any(n))[:12]
    errors = []
    for counts in [(4,3,2), (8,6,4), (12,9,6)]:
        domain = mesh.create_box(MPI.COMM_WORLD, [np.zeros(3), lengths], counts,
                                 cell_type=mesh.CellType.tetrahedron)
        space = fem.functionspace(domain, ('Lagrange', 2))
        p, q = ufl.TrialFunction(space), ufl.TestFunction(space)
        # Recover both matrices from the exact operator used by production.
        zero = bulk_helmholtz_form(p, q, 0.)
        unit = bulk_helmholtz_form(p, q, 1.)
        stiffness, mass = matrix(zero), matrix(zero-unit)
        eigenvalues = np.sort(eigsh(stiffness, k=13, M=mass, sigma=-1., which='LM',
                                   return_eigenvectors=False, tol=1e-10))
        assert abs(eigenvalues[0]) < 1e-7, 'Rigid cavity must have a constant zero mode'
        observed = 343./(2*np.pi)*np.sqrt(eigenvalues[1:])
        error = np.abs(observed/expected-1)
        errors.append(error)
        print(counts, 'maximum relative modal frequency error:', max(error))
    assert max(errors[-1]) < .01
    assert np.all(errors[1] < errors[0])
    assert np.all(errors[2] < errors[1])
