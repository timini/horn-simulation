"""Experimental Galerkin coupling to axisymmetric Rayleigh aperture modes.

The exterior is an infinite rigid baffle. This module is currently a numerical
qualification path, not an unvalidated replacement for the production default.
"""
import numpy as np
from scipy.sparse import bmat, csr_matrix
from scipy.sparse.linalg import splu
from scipy.special import j0
from dolfinx import fem
from dolfinx.fem import petsc as fem_petsc
from petsc4py import PETSc
import ufl
from horn_core.modal_radiation import modal_radiation_impedance, radial_roots


def aperture_projection(V, ds, outlet_tag, radius, modes):
    if V.mesh.comm.size!=1:
        raise ValueError('Modal aperture qualification currently requires one MPI rank')
    x=ufl.SpatialCoordinate(V.mesh)
    radial=ufl.sqrt(x[0]**2+x[1]**2)
    test=ufl.TestFunction(V)
    measure=ds(outlet_tag,metadata={'quadrature_degree':16})
    columns=[]
    for root in radial_roots(modes):
        shape=1. if root==0 else ufl.bessel_J(0,float(root)*radial/radius)/float(j0(root))
        vector=fem_petsc.assemble_vector(fem.form(ufl.inner(shape,test)*measure))
        vector.ghostUpdate(addv=PETSc.InsertMode.ADD,mode=PETSc.ScatterMode.REVERSE)
        columns.append(vector.array.copy())
        vector.destroy()
    return csr_matrix(np.column_stack(columns))


def solve_modal_aperture(V, a, L, bcs, projection, radius, k, *, rho=1.225,c=343.):
    """Couple FEM pressure and modal rho*c*velocity in one linear system."""
    if V.mesh.comm.size!=1:
        raise ValueError('Modal aperture qualification currently requires one MPI rank')
    form=fem.form(a)
    matrix=fem_petsc.assemble_matrix(form,bcs=bcs);matrix.assemble()
    indptr,indices,values=matrix.getValuesCSR()
    bulk=csr_matrix((values.copy(),indices.copy(),indptr.copy()),shape=matrix.getSize())
    vector=fem_petsc.assemble_vector(fem.form(L))
    fem_petsc.apply_lifting(vector,[form],bcs=[bcs])
    vector.ghostUpdate(addv=PETSc.InsertMode.ADD,mode=PETSc.ScatterMode.REVERSE)
    fem_petsc.set_bc(vector,bcs)
    rhs=vector.array.copy()
    matrix.destroy();vector.destroy()
    impedance=modal_radiation_impedance(k*radius,projection.shape[1])
    weak_impedance=np.pi*radius**2*impedance
    # Multiplying the modal pressure-continuity equation by i*k makes
    # the assembled block matrix complex symmetric, with consistent units.
    system=bmat([[bulk,1j*k*projection],
                 [1j*k*projection.T,csr_matrix(-1j*k*weak_impedance)]],format='csc')
    forcing=np.r_[rhs,np.zeros(projection.shape[1],dtype=complex)]
    result=splu(system).solve(forcing)
    residual=system@result-forcing
    relative=float(np.linalg.norm(residual)/max(np.linalg.norm(forcing),1e-30))
    n=V.dofmap.index_map.size_local*V.dofmap.index_map_bs
    pressure=fem.Function(V);pressure.x.array[:]=result[:n];pressure.x.scatter_forward()
    scaled_velocity=result[n:]
    modal_pressure=projection.T@result[:n]
    expected=weak_impedance@scaled_velocity
    closure=float(np.linalg.norm(modal_pressure-expected)/max(np.linalg.norm(modal_pressure),np.linalg.norm(expected),1e-30))
    if not np.isfinite(result).all() or relative>1e-8 or closure>1e-8:
        raise RuntimeError(f'Unreliable modal solve: residual={relative:g}, interface={closure:g}')
    velocity=scaled_velocity/(rho*c)
    power=float(np.real(np.vdot(velocity,weak_impedance@(rho*c*velocity))))
    return pressure,velocity,dict(relative_residual=relative,interface_relative_error=closure,radiated_power_w=power)
