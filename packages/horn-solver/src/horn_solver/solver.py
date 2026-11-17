from typing import Tuple, Dict, Any, Optional
from pathlib import Path
import pandas as pd
import numpy as np
import gmsh
from dolfinx.io import gmshio

# Define physical group tags for boundaries
INLET_TAG, OUTLET_TAG, WALL_TAG = 2, 3, 4

# Physical constants
from horn_core.duct import DEFAULT_AIR, AirProperties, boundary_layer_depths, circular_pipe_radiation

C0 = 343.0  # Speed of sound in air (m/s)
RHO0 = 1.225 # Density of air (kg/m^3)

from dolfinx import mesh, fem
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
import ufl
from petsc4py.PETSc import ScalarType
from petsc4py import PETSc

from horn_core.acoustics import validate_band, SCHEMA_VERSION
from horn_core.cad_units import import_step_metres

from horn_solver.radiation import (
    piston_radiation_impedance,
    unflanged_radiation_impedance,
)

from horn_solver.meshing import generate_qualified_mesh

RADIATION_MODELS = ("plane_wave", "flanged_piston", "modal_baffled", "unflanged_piston", "finite_flange", "closed")


def bulk_helmholtz_form(pressure, test, wavenumber):
    """Production rigid-air volume operator, shared with modal verification."""
    return (ufl.inner(ufl.grad(pressure), ufl.grad(test))
            - wavenumber**2 * ufl.inner(pressure, test)) * ufl.dx


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
    # Same motor and chamber network as the analysis coupling
    # (horn_core.driver_interface). Normal velocity = U_throat / S_throat.
    from horn_core.driver_interface import driver_operating_point
    point = driver_operating_point(driver, np.array([frequency], dtype=float),
                                   np.array([z_horn_real + 1j * z_horn_imag]), throat_area, v_g)
    return complex(point["throat_volume_velocity"][0] / throat_area)

def _mesh_columns(domain) -> Dict[str, Any]:
    """Mesh actually used, recorded on every response row."""
    info = getattr(domain, "horn_mesh_info", None)
    if not info:
        return {}
    columns = {
        "mesh_size_m": info["mesh_size_m"],
        "mesh_throat_size_m": info["throat_mesh_size_m"],
        "mesh_nodes": info["mesh_nodes"],
        "mesh_throat_area_ratio": info["throat_area_ratio"],
        "mesh_elements_per_wavelength": info["elements_per_wavelength"],
        "mesh_target_elements_per_wavelength": info["target_elements_per_wavelength"],
        "mesh_coarsened_for_node_limit": info["mesh_coarsened_for_node_limit"],
    }
    return {key: value for key, value in columns.items()
            if isinstance(value, (bool, np.bool_)) or np.isfinite(value)}


def create_mesh_from_step(step_file: str, mesh_size: float, horn_length: float, geometry_order: int = 1, *,
                          require_axisymmetric_disk=False, max_freq: Optional[float] = None,
                          sound_speed: float = DEFAULT_AIR.c, max_nodes: Optional[int] = None) -> Tuple["mesh.Mesh", "mesh.MeshTags"]:
    """Generate a qualified mesh from a STEP file and tag the boundaries.

    With ``max_freq`` the element size follows the degree's wavelength rule
    (see horn_solver.meshing); without it ``mesh_size`` is used as given.
    Both apply the throat refinement field, the 1 % throat-area guard and the
    node limit. Raises a MeshRejection subclass for a per-candidate rejection.
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("horn")
    # Honour the unit declared in the STEP file and work in metres.
    import_step_metres(gmsh, step_file, horn_length)

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
        bbox = gmsh.model.getBoundingBox(surface[0], surface[1])
        # Identify surfaces by their Z-coordinate.
        # Inlet (throat) is at z=0, Outlet (mouth) is at z=horn_length.
        if np.allclose([bbox[2], bbox[5]], 0.0, atol=1e-6):
            inlet_surfaces.append(surface[1])
        elif np.allclose([bbox[2], bbox[5]], horn_length, atol=1e-6, rtol=0):
            outlet_surfaces.append(surface[1])
        else:
            wall_surfaces.append(surface[1])
    
    if not inlet_surfaces or not outlet_surfaces or not wall_surfaces:
        gmsh.finalize()
        raise ValueError("STEP must contain an acoustic volume with planar inlet at z=0 and mouth at z=length")
    cad_areas = {
        "inlet": sum(gmsh.model.occ.getMass(2, t) for t in inlet_surfaces),
        "mouth": sum(gmsh.model.occ.getMass(2, t) for t in outlet_surfaces),
    }
    if require_axisymmetric_disk:
        from horn_solver.modal_geometry import verify_axisymmetric_disk_horn
        try:
            verify_axisymmetric_disk_horn(volumes, inlet_surfaces, outlet_surfaces)
        except Exception:
            gmsh.finalize()
            raise
    gmsh.model.addPhysicalGroup(2, inlet_surfaces, INLET_TAG)
    gmsh.model.setPhysicalName(2, INLET_TAG, "inlet")
    gmsh.model.addPhysicalGroup(2, outlet_surfaces, OUTLET_TAG)
    gmsh.model.setPhysicalName(2, OUTLET_TAG, "outlet")
    gmsh.model.addPhysicalGroup(2, wall_surfaces, WALL_TAG)
    gmsh.model.setPhysicalName(2, WALL_TAG, "walls")

    try:
        # A fixed size is expressed as an equivalent frequency for the rule.
        rule_freq = max_freq if max_freq is not None else sound_speed / (10.0 * mesh_size)
        mesh_info = generate_qualified_mesh(
            gmsh, inlet_surfaces, cad_areas["inlet"], user_size=mesh_size, c=sound_speed,
            f_max=rule_freq, degree=geometry_order, max_nodes=max_nodes)
    except Exception:
        gmsh.finalize()
        raise
    if max_freq is None:
        mesh_info["elements_per_wavelength"] = mesh_info["target_elements_per_wavelength"] = float("nan")
    print(f"Mesh: {mesh_info['mesh_nodes']} nodes (limit {mesh_info['max_mesh_nodes']}), "
          f"size {mesh_info['mesh_size_m']:.4g} m, throat size {mesh_info['throat_mesh_size_m']:.4g} m, "
          f"{mesh_info['elements_per_wavelength']:.3g} elements per wavelength, "
          f"throat area ratio {mesh_info['throat_area_ratio']:.5f}")

    # Note the change here: we are now capturing all three return values.
    domain, cell_tags, facet_tags = gmshio.model_to_mesh(gmsh.model, MPI.COMM_WORLD, 0, gdim=3)
    gmsh.finalize()

    domain.horn_boundary_areas = cad_areas
    domain.horn_axisymmetric_disk = require_axisymmetric_disk
    domain.horn_mesh_info = mesh_info
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
    if kwargs.get("radiation_model") == "bem":
        raise NotImplementedError("Legacy FEM-BEM horn coupling was removed: its whole-boundary trace does not represent a mouth-only exterior problem")

    horn_length = driver_params.get("length", 0.4)
    max_nodes = kwargs.pop("max_mesh_nodes", None)

    # Degree-dependent wavelength rule, throat grading and node limit:
    # see horn_solver.meshing.
    domain, facet_tags = create_mesh_from_step(step_file, mesh_size, horn_length,
                                               geometry_order=kwargs.get("element_degree", 1),
                                               require_axisymmetric_disk=kwargs.get("radiation_model")=="modal_baffled",
                                               max_freq=max_freq_mesh,
                                               sound_speed=kwargs.get("air", DEFAULT_AIR).c,
                                               max_nodes=max_nodes)
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
    radiation_model: str = "plane_wave",
    loss_model: str = "lossless",
    flange_width: float = 0.,
    air: AirProperties = DEFAULT_AIR,
    element_degree: int = 1,
    minimum_wall_scale: Optional[float] = None,
    inlet_velocity_rms: complex = 1.0,
) -> Path:
    """Run the FEM simulation for the Helmholtz equation.

    Args:
        domain: FEniCSx mesh object.
        facet_tags: Boundary tags from gmsh.
        freq_range: (min_freq, max_freq) in Hz.
        num_intervals: Number of frequency points.
        driver_params: Dict with at least ``length`` key.
        output_file: Path for output CSV.
        bc_mode: ``"dirichlet"`` (p=1 at inlet), ``"neumann"``
                 (velocity from driver model), or ``"velocity"``
                 (uniform prescribed inward RMS velocity).
        inlet_velocity_rms: Complex inward velocity in m/s for velocity mode.
        driver: DriverParameters instance (required for neumann mode).
        throat_area: Physical throat cross-section in m² (required for neumann).
        z_horn_initial: Dict with ``frequencies``, ``z_real``, ``z_imag`` arrays
                        from a prior Dirichlet run (required for neumann).
        radiation_model: Radiation impedance model at the outlet.
                         ``"plane_wave"`` (default, Z=rho*c),
                         ``"flanged_piston"`` (analytical piston in infinite baffle),
                         ``"modal_baffled"`` (16 axisymmetric aperture modes in
                         an infinite baffle), ``"unflanged_piston"``
                         (Levine-Schwinger approximation), ``"finite_flange"``
                         or ``"closed"``.

    Returns:
        Path to the output CSV file.
    """
    if loss_model not in {"lossless", "boundary_layer"}:
        raise ValueError("FEM supports lossless or boundary_layer walls")
    if radiation_model == "modal_baffled":
        if loss_model != "lossless" or element_degree != 1 or domain.comm.size != 1:
            raise ValueError("Experimental modal aperture requires lossless P1 and one MPI rank")
        if not getattr(domain, "horn_axisymmetric_disk", False):
            raise ValueError("Modal aperture requires a verified axisymmetric CAD volume with disk ports")
    if element_degree not in (1, 2):
        raise ValueError("Supported finite-element degrees are 1 and 2")
    if not np.isfinite(flange_width) or flange_width < 0:
        raise ValueError("Flange width must be finite and nonnegative")
    if bc_mode == "neumann" and air != DEFAULT_AIR:
        raise ValueError("Driver Neumann coupling with custom air is not supported")
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

    if bc_mode not in {"dirichlet", "neumann", "velocity"}:
        raise ValueError("Unknown inlet boundary condition")
    if bc_mode == "velocity":
        inlet_velocity_rms = complex(inlet_velocity_rms)
        if not np.isfinite(inlet_velocity_rms) or abs(inlet_velocity_rms) <= 1e-30:
            raise ValueError("Prescribed inlet velocity must be finite and nonzero")
    if radiation_model not in RADIATION_MODELS:
        raise ValueError("Unknown radiation model")
    min_freq, max_freq = freq_range
    validate_band(min_freq, max_freq)
    if num_intervals < 2:
        raise ValueError("At least two frequency points are required")
    frequencies = np.geomspace(min_freq, max_freq, num_intervals)

    # --- FEM Problem Setup ---
    V = fem.functionspace(domain, ("Lagrange", element_degree))

    ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags)
    one = fem.Constant(domain, ScalarType(1.0))
    outlet_area = fem.assemble_scalar(fem.form(one * ds(OUTLET_TAG)))
    outlet_area = domain.comm.allreduce(outlet_area, op=MPI.SUM).real
    if not np.isfinite(outlet_area) or outlet_area <= 0:
        raise ValueError("Mesh has no valid mouth boundary")
    print(f"Outlet surface area: {outlet_area:.6f} m^2")

    # Equivalent circular mouth radius for radiation impedance models
    cad_areas = getattr(domain, "horn_boundary_areas", {})
    physical_mouth_area = cad_areas.get("mouth", outlet_area)
    if (not np.isfinite(physical_mouth_area) or physical_mouth_area <= 0
            or not .5 <= outlet_area / physical_mouth_area <= 2.):
        raise RuntimeError(
            f"Integrated mesh mouth area {outlet_area:g} m² is inconsistent with "
            f"CAD area {physical_mouth_area:g} m²; check mesh and numerical backend"
        )
    a_mouth = np.sqrt(physical_mouth_area / np.pi)
    if radiation_model == "modal_baffled":
        if 2*np.pi*max_freq*a_mouth/air.c > 30:
            raise ValueError("Modal aperture quadrature requires ka <= 30")
        from horn_solver.modal_boundary import aperture_projection, solve_modal_aperture
        modal_projection = aperture_projection(V, ds, OUTLET_TAG, a_mouth, 16)
    print(f"Equivalent mouth radius: {a_mouth:.4f} m (radiation_model={radiation_model})")

    if loss_model == "boundary_layer":
        if minimum_wall_scale is None or not np.isfinite(minimum_wall_scale) or minimum_wall_scale <= 0:
            raise ValueError("Boundary-layer model requires a declared positive minimum wall curvature/gap scale")
        depths = boundary_layer_depths(min_freq, air)
        if max(depths)/minimum_wall_scale > .1:
            raise ValueError("Boundary-layer depth / minimum wall scale must be <= 0.1")
    results = []

    mode_label = f"({bc_mode} BC)"
    print(f"\nRunning simulation {mode_label} for {num_intervals} frequencies "
          f"from {min_freq:.1f}Hz to {max_freq:.1f}Hz...")

    for frequency in frequencies:
        print(f"    Solving for {frequency:.2f} Hz...")

        p = ufl.TrialFunction(V)
        q = ufl.TestFunction(V)

        omega = 2 * np.pi * frequency
        k = omega / air.c

        # Helmholtz weak form
        a = bulk_helmholtz_form(p, q, k)

        if loss_model == "boundary_layer":
            dv, dt = boundary_layer_depths(frequency, air)
            normal = ufl.FacetNormal(domain)
            tangent_p = ufl.grad(p)-ufl.dot(ufl.grad(p), normal)*normal
            tangent_q = ufl.grad(q)-ufl.dot(ufl.grad(q), normal)*normal
            # Berggren et al. Eq.31, conjugated test functions, exp(+i wt).
            a += ((1j-1)*dv/2 * ufl.inner(tangent_p, tangent_q)
                  +(1j-1)*(air.gamma-1)*dt*k**2/2 * ufl.inner(p,q))*ds(WALL_TAG)

        # Robin BC at outlet (radiation impedance); the modal model couples
        # its nonlocal aperture impedance in the solve below instead.
        if radiation_model == "modal_baffled":
            pass
        elif radiation_model == "closed":
            z_specific = None  # rigid termination; exactly zero mouth flow
        elif radiation_model == "finite_flange":
            z_specific = complex(circular_pipe_radiation(k, a_mouth, flange_width))
            a += (1j * k / z_specific) * ufl.inner(p, q) * ds(OUTLET_TAG)
        elif radiation_model == "flanged_piston":
            z_specific = piston_radiation_impedance(k, a_mouth)
            a += (1j * k / z_specific) * ufl.inner(p, q) * ds(OUTLET_TAG)
        elif radiation_model == "unflanged_piston":
            if k * a_mouth >= 1.5:
                raise ValueError("Unflanged approximation requires ka < 1.5")
            z_specific = unflanged_radiation_impedance(k, a_mouth)
            a += (1j * k / z_specific) * ufl.inner(p, q) * ds(OUTLET_TAG)
        else:
            z_specific = 1.0 + 0.0j  # plane_wave: Z = rho*c
            a += (1j * k / z_specific) * ufl.inner(p, q) * ds(OUTLET_TAG)

        bcs = []

        if bc_mode in {"neumann", "velocity"}:
            # --- Neumann mode (Phase B): prescribed velocity at inlet ---
            # Compute v_n from driver model using Phase A impedance
            if bc_mode == "velocity":
                v_n = inlet_velocity_rms
            else:
                z_r = float(_z_real_interp(frequency))
                z_i = float(_z_imag_interp(frequency))
                v_n = _compute_neumann_velocity(
                    frequency, driver, throat_area, z_r, z_i,
                )

            # Neumann source: dp/dn = +jωρ₀·v_n
            # In the weak form, the boundary integral becomes:
            # L += +jωρ₀·v_n · q · ds(INLET)
            neumann_val = 1j * omega * air.rho * v_n
            L = (ufl.inner(fem.Constant(domain, ScalarType(0.0)), q) * ufl.dx
                 + ufl.inner(fem.Constant(domain, ScalarType(neumann_val)), q) * ds(INLET_TAG))
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
        if radiation_model == "modal_baffled":
            p_h, modal_velocity, modal_health = solve_modal_aperture(
                V, a, L, bcs, modal_projection, a_mouth, k, rho=air.rho, c=air.c)
            relative_residual = modal_health['relative_residual']
            # The direct solve has an explicitly checked residual.
            converged_reason = 4
        else:
            problem = LinearProblem(
                a, L, bcs=bcs,
                petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
            )
            p_h = problem.solve()
            residual_vector = problem.b.duplicate()
            try:
                problem.A.mult(p_h.x.petsc_vec, residual_vector)
                residual_vector.axpy(-1.0, problem.b)
                relative_residual = float(residual_vector.norm()/max(problem.b.norm(),1e-30))
            finally:
                residual_vector.destroy()
            converged_reason = int(problem.solver.getConvergedReason())
            if not np.isfinite(relative_residual) or relative_residual > 1e-8 or converged_reason <= 0:
                raise RuntimeError(f"Unreliable solve at {frequency:g} Hz: residual={relative_residual:g}, PETSc reason={converged_reason}")

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
        if not np.isfinite(inlet_area_val) or inlet_area_val <= 0:
            raise ValueError("Mesh has no valid inlet boundary")
        if bc_mode == "dirichlet":
            # Recover integrated Dirichlet flux from the variational reaction.
            # Direct P1 boundary gradients are only first-order and introduce
            # large impedance errors near a throat-velocity minimum.
            reaction = fem.petsc.assemble_vector(fem.form(ufl.action(a, p_h) - L))
            reaction.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
            owned = V.dofmap.index_map.size_local * V.dofmap.index_map_bs
            inlet_owned = inlet_dofs[inlet_dofs < owned]
            dp_dn_integral = domain.comm.allreduce(np.sum(reaction.array[inlet_owned]), op=MPI.SUM)
            reaction.destroy()
        physical_inlet_area = cad_areas.get("inlet", inlet_area_val)
        # Specific impedance uses the physical CAD area; integrated FEM flux
        # is retained. This keeps p/U identical across impedance domains.
        dp_dn_avg = dp_dn_integral / physical_inlet_area if physical_inlet_area > 0 else 0.0

        if bc_mode == "dirichlet":
            if inlet_area_val <= 0:
                raise ValueError("Mesh has no valid inlet boundary")
            # Inward velocity = dp/dn / (i omega rho), exp(+i wt).
            if abs(dp_dn_avg) > 1e-30:
                z_throat = (1j * omega * air.rho) / dp_dn_avg
            else:
                raise FloatingPointError("Cannot recover impedance from zero inlet flow")
        else:
            # In Neumann mode, compute impedance from average inlet pressure
            p_inlet_integral = fem.assemble_scalar(fem.form(p_h * ds(INLET_TAG)))
            p_inlet_integral = domain.comm.allreduce(p_inlet_integral, op=MPI.SUM)
            p_inlet_avg = p_inlet_integral / inlet_area_val if inlet_area_val > 0 else 0.0
            if abs(v_n) > 1e-30:
                z_throat = p_inlet_avg * physical_inlet_area / (v_n * inlet_area_val)
            else:
                raise FloatingPointError("Cannot recover impedance from zero inlet flow")

        z_real = float(np.real(z_throat))
        z_imag = float(np.imag(z_throat))

        print(f"      -> p_rms = {p_rms:.4f}, SPL = {spl:.2f} dB, "
              f"phase = {phase_deg:.1f}, Z = {z_real:.1f} + {z_imag:.1f}j")
        # Local Robin condition gives outward mouth volume velocity exactly
        # in its own model.
        if radiation_model == "modal_baffled":
            mouth_u = physical_mouth_area * modal_velocity[0]
        elif radiation_model == "closed":
            mouth_u = 0j
        else:
            mouth_u = p_outlet_integral / (air.rho * air.c * z_specific)
        input_u = (dp_dn_integral/(1j*omega*air.rho) if bc_mode == "dirichlet"
                   else v_n*inlet_area_val)
        input_p = 1.0 if bc_mode == "dirichlet" else p_inlet_avg
        input_power = float(np.real(input_p*np.conjugate(input_u)))
        output_power = (modal_health['radiated_power_w'] if radiation_model == "modal_baffled" else
                        0. if radiation_model == "closed" else
                        float(p_outlet_sq*np.real(1/z_specific)/(air.rho*air.c)))
        viscous_power = thermal_power = 0.
        if loss_model == "boundary_layer":
            tangent_h = ufl.grad(p_h)-ufl.dot(ufl.grad(p_h), normal)*normal
            gv = fem.assemble_scalar(fem.form(ufl.inner(tangent_h,tangent_h)*ds(WALL_TAG)))
            gt = fem.assemble_scalar(fem.form(ufl.inner(p_h,p_h)*ds(WALL_TAG)))
            viscous_power = float(domain.comm.allreduce(gv, op=MPI.SUM).real*dv/(2*omega*air.rho))
            thermal_power = float(domain.comm.allreduce(gt, op=MPI.SUM).real*(air.gamma-1)*dt*omega/(2*air.rho*air.c**2))
        results.append({
            "frequency": frequency,
            "spl": spl,
            "phase_deg": phase_deg,
            "z_real": z_real,
            "z_imag": z_imag,
            "schema_version": SCHEMA_VERSION,
            "inlet_p_real": float(np.real(input_p)),
            "inlet_p_imag": float(np.imag(input_p)),
            "inlet_u_real": float(np.real(input_u)),
            "inlet_u_imag": float(np.imag(input_u)),
            "relative_residual": relative_residual,
            "converged_reason": converged_reason,
            "inlet_area_m2": float(physical_inlet_area),
            "mesh_inlet_area_m2": float(inlet_area_val),
            "mouth_area_m2": float(physical_mouth_area),
            "mesh_mouth_area_m2": float(outlet_area),
            "mouth_u_real": float(mouth_u.real),
            "mouth_u_imag": float(mouth_u.imag),
            "mouth_p_real": float(p_avg.real),
            "mouth_p_imag": float(p_avg.imag),
            "radiation_model": radiation_model,
            "loss_model": loss_model,
            "flange_width_m": flange_width,
            "air_c_m_s": air.c,
            "air_rho_kg_m3": air.rho,
            "air_gamma": air.gamma,
            "air_viscosity_pa_s": air.viscosity,
            "air_conductivity_w_m_k": air.conductivity,
            "air_heat_capacity_j_kg_k": air.heat_capacity,
            "element_degree": element_degree,
            "input_acoustic_power_w": input_power,
            "mouth_acoustic_power_w": output_power,
            "viscous_wall_power_w": viscous_power,
            "thermal_wall_power_w": thermal_power,
            "power_balance_available": True,
            "bc_mode": bc_mode,
            "phasor_convention": "exp(+iwt)_rms",
            "mesh_cells": domain.topology.index_map(domain.topology.dim).size_global,
            **_mesh_columns(domain),
        })
        if radiation_model == "modal_baffled":
            results[-1]['modal_mode_count'] = len(modal_velocity)
            results[-1]['modal_interface_relative_error'] = modal_health['interface_relative_error']
            for index, value in enumerate(modal_velocity):
                results[-1][f'modal_v_{index}_real'] = float(value.real)
                results[-1][f'modal_v_{index}_imag'] = float(value.imag)

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
    parser.add_argument("--mesh-size", type=float, default=0.01,
                        help="Largest element size (m); the finer of this and the degree's wavelength rule is used.")
    parser.add_argument("--max-mesh-nodes", type=int, default=None,
                        help="Node limit (default: HORN_MAX_MESH_NODES or 100000).")
    parser.add_argument("--length", type=float, required=True, help="Length of the horn for boundary tagging.")
    parser.add_argument("--bc-mode", type=str, default="dirichlet",
                        choices=["dirichlet", "neumann", "velocity"],
                        help="Boundary condition mode (default: dirichlet).")
    parser.add_argument("--inlet-velocity-rms", type=complex, default=1.0,
                        help="Uniform inward RMS velocity in m/s for velocity mode (complex values accepted).")
    parser.add_argument("--driver-json", type=str, default=None,
                        help="Driver database JSON (required for neumann mode).")
    parser.add_argument("--driver-id", type=str, default=None,
                        help="Driver ID in database (required for neumann mode).")
    parser.add_argument("--throat-area", type=float, default=None,
                        help="Throat cross-section area in m² (required for neumann mode).")
    parser.add_argument("--phase-a-csv", type=str, default=None,
                        help="Phase A solver CSV with Z_horn data (required for neumann mode).")
    parser.add_argument("--radiation-model", type=str, default="plane_wave",
                        choices=list(RADIATION_MODELS),
                        help="Radiation impedance model at the outlet (default: plane_wave).")
    parser.add_argument("--rejection-file", type=str, default=None,
                        help="On a per-candidate mesh or CAD rejection, write the reason to this JSON "
                             "file and exit 0 instead of failing.")
    parser.add_argument("--loss-model", choices=["lossless", "boundary_layer"], default="lossless")
    parser.add_argument("--flange-width", type=float, default=0.)
    parser.add_argument("--minimum-wall-scale", type=float)
    parser.add_argument("--element-degree", type=int, default=1)
    args = parser.parse_args()

    # Build extra kwargs for neumann mode
    extra_kwargs = {
        "loss_model": args.loss_model,
        "flange_width": args.flange_width,
        "minimum_wall_scale": args.minimum_wall_scale,
        "element_degree": args.element_degree,
        "bc_mode": args.bc_mode,
        "inlet_velocity_rms": args.inlet_velocity_rms,
        "radiation_model": args.radiation_model,
        "max_mesh_nodes": args.max_mesh_nodes,
    }

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

    try:
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
    except Exception as error:
        from horn_solver.meshing import rejection_reason, write_rejection
        if args.rejection_file is None or rejection_reason(error) is None:
            raise
        record = write_rejection(args.rejection_file, error, min_freq_hz=args.min_freq,
                                 max_freq_hz=args.max_freq)
        print(f"Candidate rejected ({record['rejection_reasons'][0]}): {error}")

if __name__ == "__main__":
    main()
