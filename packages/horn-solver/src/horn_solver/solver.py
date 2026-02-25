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

from dolfinx import mesh, fem
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
import ufl
from petsc4py.PETSc import ScalarType

def create_mesh_from_step(step_file: str, mesh_size: float, horn_length: float) -> Tuple["mesh.Mesh", "mesh.MeshTags"]:
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
        # Identify surfaces by their Z-coordinate.
        # Inlet (throat) is at z=0, Outlet (mouth) is at z=horn_length.
        if np.isclose(com[2], 0.0):
            inlet_surfaces.append(surface[1])
        elif np.isclose(com[2], horn_length):
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
    num_intervals: int,
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
    # We need the horn length to correctly identify the outlet boundary.
    # The length is a parameter of the geometry generation. For now, let's
    # extract it from the Nextflow parameters. This is a bit of a hack.
    # A better solution would be to save it as a metadata file.
    # We will assume it's available in the arguments for now.
    from os import environ
    # A bit of a hack: read the length from the params passed to the script.
    # This is not ideal. A better way would be to pass it in the channel.
    # For now, we'll rely on the `--length` argument being available.
    # Since we are not in main, we will just use a fixed value.
    # This will be passed from main.nf
    horn_length = driver_params.get("length", 0.4) # Default for safety

    # Frequency-adaptive meshing: element size < λ/6 = c/(6*f_max)
    adaptive_size = C0 / (6.0 * max_freq_mesh)
    final_mesh_size = min(mesh_size, adaptive_size)
    print(f"Mesh size: user={mesh_size}, adaptive={adaptive_size:.4f}, using={final_mesh_size:.4f}")

    domain, facet_tags = create_mesh_from_step(step_file, final_mesh_size, horn_length)
    print(f"Successfully loaded mesh: {domain.name} with {domain.topology.index_map(domain.topology.dim).size_global} cells.")

    # --- Solving ---
    return run_simulation(domain, facet_tags, freq_range, num_intervals, driver_params, output_file)

def run_simulation(
    domain: "mesh.Mesh",
    facet_tags: "mesh.MeshTags",
    freq_range: Tuple[float, float],
    num_intervals: int,
    driver_params: Dict[str, Any],
    output_file: str,
) -> Path:
    """
    Runs the FEM simulation for the Helmholtz equation using a dolfinx mesh object.
    """
    min_freq, max_freq = freq_range
    # Let's use a logarithmic frequency sweep for more detail at lower frequencies
    frequencies = np.geomspace(min_freq, max_freq, num_intervals)
    
    # --- FEM Problem Setup ---
    # Define function space for complex pressure
    V = fem.functionspace(domain, ("Lagrange", 1))

    # --- Boundary measure and outlet area (computed once, mesh doesn't change) ---
    ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags)
    one = fem.Constant(domain, ScalarType(1.0))
    outlet_area = fem.assemble_scalar(fem.form(one * ds(OUTLET_TAG)))
    outlet_area = domain.comm.allreduce(outlet_area, op=MPI.SUM).real
    print(f"Outlet surface area: {outlet_area:.6f} m^2")

    # --- Store results ---
    results = []

    print(f"\nRunning simulation for {num_intervals} frequencies from {min_freq:.1f}Hz to {max_freq:.1f}Hz...")
    
    for frequency in frequencies:
        print(f"    Solving for {frequency:.2f} Hz...")
        
        # Define trial and test functions
        p = ufl.TrialFunction(V)
        q = ufl.TestFunction(V)

        # Physical parameters
        omega = 2 * np.pi * frequency # angular frequency
        k = omega / C0 # wave number

        # Define the variational problem (weak form of Helmholtz equation)
        # LHS: Integral over the domain volume. Use ufl.inner for complex problems.
        a = ufl.inner(ufl.grad(p), ufl.grad(q)) * ufl.dx - k**2 * ufl.inner(p, q) * ufl.dx

        # Radiation impedance BC at outlet (first-order Sommerfeld condition: dp/dn = -jkp)
        a -= 1j * k * ufl.inner(p, q) * ds(OUTLET_TAG)

        # RHS: Define the driver velocity at the inlet.
        # We now use a strong Dirichlet condition to enforce p=1 at the inlet,
        # mirroring the successful test case.
        L = ufl.inner(fem.Constant(domain, ScalarType(0.0)), q) * ufl.dx
        
        # --- Boundary Conditions ---
        # Enforce p=1 at the inlet
        inlet_pressure = fem.Function(V)
        inlet_pressure.x.array[:] = 1.0 + 0j
        
        # Locate facets for inlet boundary using the facet_tags from GMSH
        inlet_facets = facet_tags.find(INLET_TAG)
        inlet_dofs = fem.locate_dofs_topological(V, domain.topology.dim - 1, inlet_facets)
        bc = fem.dirichletbc(inlet_pressure, inlet_dofs)

        # --- Solve ---
        # Set up and solve the linear problem. The JIT compilation will happen automatically.
        problem = LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        p_h = problem.solve()

        # Robin BC (radiation impedance) is applied in the bilinear form above

        # --- Post-processing ---
        # Integrate |p|^2 over the outlet surface for physically meaningful SPL
        p_outlet_sq = fem.assemble_scalar(fem.form(ufl.inner(p_h, p_h) * ds(OUTLET_TAG)))
        p_outlet_sq = domain.comm.allreduce(p_outlet_sq, op=MPI.SUM).real

        # RMS pressure at outlet
        p_rms = np.sqrt(p_outlet_sq / outlet_area) if outlet_area > 0 else 0.0

        # SPL relative to 20 µPa reference
        p_ref = 20e-6
        spl = 20 * np.log10(p_rms / p_ref + 1e-12)

        print(f"      -> p_rms = {p_rms:.4f}, SPL = {spl:.2f} dB")
        results.append({"frequency": frequency, "spl": spl})

    # --- Output Generation ---
    output_path = Path(output_file)
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_path, index=False)

    print(f"\nSuccessfully ran simulation and generated results: {output_path}")
    
    return output_path

def main():
    """Command-line interface for the solver."""
    import argparse
    parser = argparse.ArgumentParser(description="Run FEM simulation on a horn STEP file.")
    parser.add_argument("--step-file", type=str, required=True, help="Path to the input STEP file.")
    parser.add_argument("--output-file", type=str, required=True, help="Path to save the output CSV results.")
    parser.add_argument("--min-freq", type=float, default=100.0, help="Minimum frequency in Hz.")
    parser.add_argument("--max-freq", type=float, default=1000.0, help="Maximum frequency in Hz.")
    parser.add_argument("--num-intervals", type=int, default=100, help="Number of frequency steps.")
    parser.add_argument("--mesh-size", type=float, default=0.01, help="Mesh element size.")
    parser.add_argument("--length", type=float, required=True, help="Length of the horn for boundary tagging.")
    args = parser.parse_args()

    run_simulation_from_step(
        step_file=args.step_file,
        freq_range=(args.min_freq, args.max_freq),
        num_intervals=args.num_intervals,
        driver_params={"length": args.length},  # Pass length to solver
        output_file=args.output_file,
        max_freq_mesh=args.max_freq, # Use max_freq for this for now
        mesh_size=args.mesh_size,
    )

if __name__ == "__main__":
    main()
