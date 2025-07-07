from typing import Tuple, Dict, Any
from pathlib import Path
import os
import pandas as pd
import numpy as np
import gmsh
from dolfinx.io import gmshio

# Define physical group tags for boundaries
INLET_TAG, OUTLET_TAG, WALL_TAG = 2, 3, 4

# Physical constants
C0 = 343.0  # Speed of sound in air (m/s)
RHO0 = 1.225 # Density of air (kg/m^3)

try:
    from dolfinx import mesh, fem
    from dolfinx.fem.petsc import LinearProblem
    from mpi4py import MPI
    import ufl
    from petsc4py.PETSc import ScalarType
except ImportError as e:
    # This is expected if running outside the solver container (e.g., for local tests)
    print(f"DEBUG: Caught ImportError: {e}")
    dolfinx = None
    fem = None
    MPI = None
    mesh = None

def create_mesh_from_step(step_file: str, mesh_size: float) -> Tuple["mesh.Mesh", "mesh.MeshTags"]:
    """
    Generates a mesh from a STEP file and tags the boundaries.
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("horn")
    gmsh.model.occ.importShapes(step_file)
    gmsh.model.occ.synchronize()

    # Tag the physical volume for meshing
    volumes = gmsh.model.occ.getEntities(dim=3)
    if not volumes:
        gmsh.model.occ.synchronize()
        volumes = gmsh.model.occ.getEntities(dim=3)
        if not volumes:
            raise RuntimeError("No 3D volumes found in the STEP file after import.")
    
    physical_volume = gmsh.model.addPhysicalGroup(3, [v[1] for v in volumes], tag=1)
    gmsh.model.setPhysicalName(3, physical_volume, "main_volume")

    # Tag the boundary surfaces for applying boundary conditions
    inlet_surfaces = []
    outlet_surfaces = []
    wall_surfaces = []

    surfaces = gmsh.model.occ.getEntities(dim=2)
    for surface in surfaces:
        com = gmsh.model.occ.getCenterOfMass(surface[0], surface[1])
        # Note: This assumes a specific orientation of the test_box.stp file
        # Y-min is inlet, Y-max is outlet.
        if np.allclose(com, [0.5, 0.0, 0.5]):
            inlet_surfaces.append(surface[1])
        elif np.allclose(com, [0.5, 1.0, 0.5]):
            outlet_surfaces.append(surface[1])
        else:
            wall_surfaces.append(surface[1])
    
    gmsh.model.addPhysicalGroup(2, inlet_surfaces, INLET_TAG)
    gmsh.model.setPhysicalName(2, INLET_TAG, "inlet")
    gmsh.model.addPhysicalGroup(2, outlet_surfaces, OUTLET_TAG)
    gmsh.model.setPhysicalName(2, OUTLET_TAG, "outlet")
    gmsh.model.addPhysicalGroup(2, wall_surfaces, WALL_TAG)
    gmsh.model.setPhysicalName(2, WALL_TAG, "walls")

    gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size)
    gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)
    gmsh.model.mesh.generate(3)
    
    # Note the change here: we are now capturing all three return values.
    domain, cell_tags, facet_tags = gmshio.model_to_mesh(gmsh.model, MPI.COMM_WORLD, 0, gdim=3)
    gmsh.finalize()

    return domain, facet_tags

def run_simulation_from_step(
    step_file: str,
    freq_range: Tuple[float, float],
    driver_params: Dict[str, Any],
    output_file: str,
    max_freq_mesh: float,
    mesh_size: float = 0.01,  # Default to a fine mesh for production
) -> Path:
    """
    Generates a mesh from a STEP file and then runs the simulation.
    """
    if not fem:
        raise ImportError("FEniCSx (dolfinx) is required for the simulation.")

    # --- Meshing ---
    domain, facet_tags = create_mesh_from_step(step_file, mesh_size)
    print(f"Successfully loaded mesh: {domain.name} with {domain.topology.index_map(domain.topology.dim).size_global} cells.")

    # --- Solving ---
    return run_simulation(domain, facet_tags, freq_range, driver_params, output_file)

def run_simulation(
    domain: "mesh.Mesh",
    facet_tags: "mesh.MeshTags",
    freq_range: Tuple[float, float],
    driver_params: Dict[str, Any],
    output_file: str,
) -> Path:
    """
    Runs the FEM simulation for the Helmholtz equation using a dolfinx mesh object.
    """
    # For now, solve for a single frequency, then we will implement the sweep
    min_freq, max_freq = freq_range
    frequency = (min_freq + max_freq) / 2.0

    # --- FEM Problem Setup ---
    # Define function space for complex pressure
    V = fem.functionspace(domain, ("Lagrange", 1))

    # Define trial and test functions
    p = ufl.TrialFunction(V)
    q = ufl.TestFunction(V)

    # Physical parameters
    omega = 2 * np.pi * frequency # angular frequency
    k = omega / C0 # wave number

    # Define the variational problem (weak form of Helmholtz equation)
    # LHS: Integral over the domain volume. Use ufl.inner for complex problems.
    a = ufl.inner(ufl.grad(p), ufl.grad(q)) * ufl.dx - k**2 * ufl.inner(p, q) * ufl.dx
    # Add Robin boundary condition for the outlet (radiation condition)
    a -= 1j * k * ufl.inner(p, q) * ufl.ds(OUTLET_TAG)

    # RHS: Define the driver velocity at the inlet. Use ufl.inner.
    u_n = fem.Constant(domain, ScalarType(1.0))
    L = -1j * RHO0 * omega * ufl.inner(u_n, q) * ufl.ds(INLET_TAG)

    # --- Solve ---
    # Set up and solve the linear problem. The JIT compilation will happen automatically.
    problem = LinearProblem(a, L, bcs=[], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    p_h = problem.solve()

    # --- Post-processing ---
    # For now, we compute a placeholder SPL from the L2-norm of the solution
    # In the future, we will integrate pressure over the outlet surface
    p_ref = 20e-6 # Reference pressure for SPL
    # Calculate L2 norm of the complex pressure field. ufl.inner(p,p) is equivalent
    # to p*conj(p) for complex fields.
    norm_p_L2_squared = fem.assemble_scalar(fem.form(ufl.inner(p_h, p_h) * ufl.dx))
    norm_p_L2 = np.sqrt(domain.comm.allreduce(norm_p_L2_squared, op=MPI.SUM))
    spl = 20 * np.log10(norm_p_L2 / p_ref)

    # --- Output Generation ---
    output_path = Path(output_file)
    results_df = pd.DataFrame({
        "frequency": [frequency],
        "spl": [spl]
    })
    results_df.to_csv(output_path, index=False)

    print(f"Successfully ran simulation and generated results: {output_path}")
    
    return output_path
