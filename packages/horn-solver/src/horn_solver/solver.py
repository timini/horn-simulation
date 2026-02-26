from typing import Tuple, Dict, Any, Optional
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


def _compute_neumann_velocity(
    frequency: float,
    driver: "DriverParameters",
    throat_area: float,
    z_horn_real: float,
    z_horn_imag: float,
    v_g: float = 2.83,
) -> complex:
    """Compute normal velocity at the inlet for Neumann BC mode.

    Uses the driver electro-mechanical model and the horn throat impedance
    (from a prior Phase A Dirichlet simulation) to compute the diaphragm
    velocity at a single frequency.

    Returns:
        Complex normal velocity v_n (m/s) at the throat.
    """
    omega = 2.0 * np.pi * frequency

    # Electrical impedance
    z_e = driver.re_ohm + 1j * omega * driver.le_h

    # Mechanical impedance of the driver suspension
    rms = driver.rms_kg_per_s if driver.rms_kg_per_s is not None else 0.0
    z_mech_driver = rms + 1j * omega * driver.mms_kg + 1.0 / (1j * omega * driver.cms_m_per_n)

    # Horn throat impedance (specific acoustic → mechanical)
    z_horn = z_horn_real + 1j * z_horn_imag
    z_mech_load = z_horn * (driver.sd_m2 ** 2) / throat_area

    z_mech_total = z_mech_driver + z_mech_load
    z_mot = (driver.bl_tm ** 2) / z_mech_total

    current = v_g / (z_e + z_mot)
    velocity = driver.bl_tm * current / z_mech_total

    # Volume velocity at throat = v_diaphragm * Sd, normal velocity = U / S_throat
    v_n = velocity * driver.sd_m2 / throat_area
    return v_n

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
    mesh_size: float = 0.01,
    **kwargs,
) -> Path:
    """Generate a mesh from a STEP file and run the simulation.

    Extra keyword arguments (``bc_mode``, ``driver``, ``throat_area``,
    ``z_horn_initial``) are forwarded to :func:`run_simulation`.
    """
    if not fem:
        raise ImportError("FEniCSx (dolfinx) is required for the simulation.")

    horn_length = driver_params.get("length", 0.4)

    # Frequency-adaptive meshing: element size < λ/6 = c/(6*f_max)
    adaptive_size = C0 / (6.0 * max_freq_mesh)
    final_mesh_size = min(mesh_size, adaptive_size)
    print(f"Mesh size: user={mesh_size}, adaptive={adaptive_size:.4f}, using={final_mesh_size:.4f}")

    domain, facet_tags = create_mesh_from_step(step_file, final_mesh_size, horn_length)
    print(f"Successfully loaded mesh: {domain.name} with "
          f"{domain.topology.index_map(domain.topology.dim).size_global} cells.")

    return run_simulation(
        domain, facet_tags, freq_range, num_intervals,
        driver_params, output_file, **kwargs,
    )

def run_simulation(
    domain: "mesh.Mesh",
    facet_tags: "mesh.MeshTags",
    freq_range: Tuple[float, float],
    num_intervals: int,
    driver_params: Dict[str, Any],
    output_file: str,
    bc_mode: str = "dirichlet",
    driver: Optional[Any] = None,
    throat_area: Optional[float] = None,
    z_horn_initial: Optional[Dict[str, np.ndarray]] = None,
) -> Path:
    """Run the FEM simulation for the Helmholtz equation.

    Args:
        domain: FEniCSx mesh object.
        facet_tags: Boundary tags from gmsh.
        freq_range: (min_freq, max_freq) in Hz.
        num_intervals: Number of frequency points.
        driver_params: Dict with at least ``length`` key.
        output_file: Path for output CSV.
        bc_mode: ``"dirichlet"`` (Phase A, p=1 at inlet) or ``"neumann"``
                 (Phase B, prescribed velocity from driver model).
        driver: DriverParameters instance (required for neumann mode).
        throat_area: Physical throat cross-section in m² (required for neumann).
        z_horn_initial: Dict with ``frequencies``, ``z_real``, ``z_imag`` arrays
                        from a prior Dirichlet run (required for neumann).

    Returns:
        Path to the output CSV file.
    """
    if bc_mode == "neumann":
        if driver is None or throat_area is None or z_horn_initial is None:
            raise ValueError(
                "Neumann BC mode requires driver, throat_area, and z_horn_initial"
            )
        from scipy.interpolate import interp1d
        _z_real_interp = interp1d(
            z_horn_initial["frequencies"], z_horn_initial["z_real"],
            kind="linear", fill_value="extrapolate",
        )
        _z_imag_interp = interp1d(
            z_horn_initial["frequencies"], z_horn_initial["z_imag"],
            kind="linear", fill_value="extrapolate",
        )

    min_freq, max_freq = freq_range
    frequencies = np.geomspace(min_freq, max_freq, num_intervals)

    # --- FEM Problem Setup ---
    V = fem.functionspace(domain, ("Lagrange", 1))

    ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags)
    one = fem.Constant(domain, ScalarType(1.0))
    outlet_area = fem.assemble_scalar(fem.form(one * ds(OUTLET_TAG)))
    outlet_area = domain.comm.allreduce(outlet_area, op=MPI.SUM).real
    print(f"Outlet surface area: {outlet_area:.6f} m^2")

    results = []

    mode_label = f"({bc_mode} BC)"
    print(f"\nRunning simulation {mode_label} for {num_intervals} frequencies "
          f"from {min_freq:.1f}Hz to {max_freq:.1f}Hz...")

    for frequency in frequencies:
        print(f"    Solving for {frequency:.2f} Hz...")

        p = ufl.TrialFunction(V)
        q = ufl.TestFunction(V)

        omega = 2 * np.pi * frequency
        k = omega / C0

        # Helmholtz weak form
        a = (ufl.inner(ufl.grad(p), ufl.grad(q)) * ufl.dx
             - k**2 * ufl.inner(p, q) * ufl.dx)

        # Robin BC at outlet (radiation impedance)
        a -= 1j * k * ufl.inner(p, q) * ds(OUTLET_TAG)

        bcs = []

        if bc_mode == "neumann":
            # --- Neumann mode (Phase B): prescribed velocity at inlet ---
            # Compute v_n from driver model using Phase A impedance
            z_r = float(_z_real_interp(frequency))
            z_i = float(_z_imag_interp(frequency))
            v_n = _compute_neumann_velocity(
                frequency, driver, throat_area, z_r, z_i,
            )

            # Neumann source: dp/dn = -jωρ₀·v_n
            # In the weak form, the boundary integral becomes:
            # L += -jωρ₀·v_n · q · ds(INLET)
            neumann_val = -1j * omega * RHO0 * v_n
            L = (ufl.inner(fem.Constant(domain, ScalarType(0.0)), q) * ufl.dx
                 + fem.Constant(domain, ScalarType(neumann_val)) * q * ds(INLET_TAG))
        else:
            # --- Dirichlet mode (Phase A): p=1 at inlet ---
            L = ufl.inner(fem.Constant(domain, ScalarType(0.0)), q) * ufl.dx

            inlet_pressure = fem.Function(V)
            inlet_pressure.x.array[:] = 1.0 + 0j
            inlet_facets = facet_tags.find(INLET_TAG)
            inlet_dofs = fem.locate_dofs_topological(
                V, domain.topology.dim - 1, inlet_facets,
            )
            bcs.append(fem.dirichletbc(inlet_pressure, inlet_dofs))

        # --- Solve ---
        problem = LinearProblem(
            a, L, bcs=bcs,
            petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
        )
        p_h = problem.solve()

        # --- Post-processing ---
        p_outlet_sq = fem.assemble_scalar(
            fem.form(ufl.inner(p_h, p_h) * ds(OUTLET_TAG))
        )
        p_outlet_sq = domain.comm.allreduce(p_outlet_sq, op=MPI.SUM).real
        p_rms = np.sqrt(p_outlet_sq / outlet_area) if outlet_area > 0 else 0.0

        p_ref = 20e-6
        spl = 20 * np.log10(p_rms / p_ref + 1e-12)

        # Phase
        p_outlet_integral = fem.assemble_scalar(fem.form(p_h * ds(OUTLET_TAG)))
        p_outlet_integral = domain.comm.allreduce(p_outlet_integral, op=MPI.SUM)
        p_avg = p_outlet_integral / outlet_area if outlet_area > 0 else 0.0
        phase_deg = float(np.degrees(np.angle(p_avg)))

        # Throat impedance
        n = ufl.FacetNormal(domain)
        dp_dn_form = fem.form(ufl.dot(ufl.grad(p_h), n) * ds(INLET_TAG))
        dp_dn_integral = fem.assemble_scalar(dp_dn_form)
        dp_dn_integral = domain.comm.allreduce(dp_dn_integral, op=MPI.SUM)
        inlet_area_val = fem.assemble_scalar(fem.form(one * ds(INLET_TAG)))
        inlet_area_val = domain.comm.allreduce(inlet_area_val, op=MPI.SUM).real
        dp_dn_avg = dp_dn_integral / inlet_area_val if inlet_area_val > 0 else 0.0

        if bc_mode == "dirichlet":
            # Z = p_inlet / v_n, with p_inlet = 1
            if abs(dp_dn_avg) > 1e-30:
                z_throat = (1j * omega * RHO0) / (-dp_dn_avg)
            else:
                z_throat = 0.0 + 0.0j
        else:
            # In Neumann mode, compute impedance from average inlet pressure
            p_inlet_integral = fem.assemble_scalar(fem.form(p_h * ds(INLET_TAG)))
            p_inlet_integral = domain.comm.allreduce(p_inlet_integral, op=MPI.SUM)
            p_inlet_avg = p_inlet_integral / inlet_area_val if inlet_area_val > 0 else 0.0
            if abs(v_n) > 1e-30:
                z_throat = p_inlet_avg / v_n
            else:
                z_throat = 0.0 + 0.0j

        z_real = float(np.real(z_throat))
        z_imag = float(np.imag(z_throat))

        print(f"      -> p_rms = {p_rms:.4f}, SPL = {spl:.2f} dB, "
              f"phase = {phase_deg:.1f}, Z = {z_real:.1f} + {z_imag:.1f}j")
        results.append({
            "frequency": frequency,
            "spl": spl,
            "phase_deg": phase_deg,
            "z_real": z_real,
            "z_imag": z_imag,
        })

    # --- Output Generation ---
    output_path = Path(output_file)
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_path, index=False)

    print(f"\nSuccessfully ran simulation and generated results: {output_path}")

    return output_path

def main():
    """Command-line interface for the solver."""
    import argparse
    import json as _json

    parser = argparse.ArgumentParser(description="Run FEM simulation on a horn STEP file.")
    parser.add_argument("--step-file", type=str, required=True, help="Path to the input STEP file.")
    parser.add_argument("--output-file", type=str, required=True, help="Path to save the output CSV results.")
    parser.add_argument("--min-freq", type=float, default=100.0, help="Minimum frequency in Hz.")
    parser.add_argument("--max-freq", type=float, default=1000.0, help="Maximum frequency in Hz.")
    parser.add_argument("--num-intervals", type=int, default=100, help="Number of frequency steps.")
    parser.add_argument("--mesh-size", type=float, default=0.01, help="Mesh element size.")
    parser.add_argument("--length", type=float, required=True, help="Length of the horn for boundary tagging.")
    parser.add_argument("--bc-mode", type=str, default="dirichlet",
                        choices=["dirichlet", "neumann"],
                        help="Boundary condition mode (default: dirichlet).")
    parser.add_argument("--driver-json", type=str, default=None,
                        help="Driver database JSON (required for neumann mode).")
    parser.add_argument("--driver-id", type=str, default=None,
                        help="Driver ID in database (required for neumann mode).")
    parser.add_argument("--throat-area", type=float, default=None,
                        help="Throat cross-section area in m² (required for neumann mode).")
    parser.add_argument("--phase-a-csv", type=str, default=None,
                        help="Phase A solver CSV with Z_horn data (required for neumann mode).")
    args = parser.parse_args()

    # Build extra kwargs for neumann mode
    extra_kwargs = {"bc_mode": args.bc_mode}

    if args.bc_mode == "neumann":
        if not all([args.driver_json, args.driver_id, args.throat_area, args.phase_a_csv]):
            parser.error(
                "Neumann mode requires --driver-json, --driver-id, "
                "--throat-area, and --phase-a-csv"
            )
        from horn_drivers.loader import load_driver
        extra_kwargs["driver"] = load_driver(args.driver_json, args.driver_id)
        extra_kwargs["throat_area"] = args.throat_area

        phase_a_df = pd.read_csv(args.phase_a_csv)
        extra_kwargs["z_horn_initial"] = {
            "frequencies": phase_a_df["frequency"].values,
            "z_real": phase_a_df["z_real"].values,
            "z_imag": phase_a_df["z_imag"].values,
        }

    run_simulation_from_step(
        step_file=args.step_file,
        freq_range=(args.min_freq, args.max_freq),
        num_intervals=args.num_intervals,
        driver_params={"length": args.length},
        output_file=args.output_file,
        max_freq_mesh=args.max_freq,
        mesh_size=args.mesh_size,
        **extra_kwargs,
    )

if __name__ == "__main__":
    main()
