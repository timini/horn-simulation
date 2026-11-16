"""Mesh sizing and qualification for the horn FEM solve.

This module needs gmsh but not dolfinx, so its rules can be tested outside
the solver container.

Element size
------------
The size target depends on the Lagrange degree:

* P1: lambda/10 at the band's highest frequency.
* P2: lambda/6.

The basis is the dispersion (pollution) error of a Helmholtz finite element
solve. For degree p the relative phase-speed error per wavelength is about
(kh)^(2p) times 1/24 (p=1) or 1/1440 (p=2). That gives:

=============  =====  ======================
degree, size   kh     phase error per lambda
=============  =====  ======================
P1, lambda/6   1.05   4.6 %
P1, lambda/10  0.63   1.6 %
P2, lambda/6   1.05   0.08 %
=============  =====  ======================

On the default single-mode horn (conical, 50/200/500 mm) at 8 kHz, gmsh gives
about 61k nodes at P1 lambda/6, 266k at P1 lambda/10 and 475k P2 nodes at
lambda/6. P2 at lambda/6 is far more accurate per wavelength but costs about
1.8 times the unknowns of P1 at lambda/10. The default stays P1 at lambda/10
because the modal radiation model (the default for axisymmetric CAD)
supports P1 only, and all executed resolution evidence is P1
(docs/CANDIDATE_RESOLUTION.md, docs/MODAL_CANDIDATE_RESOLUTION.md).

When the lambda/10 mesh would exceed the node limit, the P1 size is coarsened
just enough to fit, but never beyond lambda/6, the previous production rule.
The solver records the elements per wavelength actually used so reports can
show it. If even lambda/6 does not fit, the candidate is rejected with reason
``mesh_node_limit_exceeded`` instead of stopping the run.

Throat
------
A gmsh Distance/Threshold field caps the element size at r_t/4 within one
throat radius of the throat disk, grading to the bulk size over two more
throat radii. The faceted mesh area of the throat disk must then be within
1 % of the CAD area, because the Dirichlet throat flux is integrated over the
mesh disk. A failure is retried once with a finer throat cap, then rejected
with reason ``throat_area_guard_failed``.
"""
import math
import os

import numpy as np

ELEMENTS_PER_WAVELENGTH = {1: 10.0, 2: 6.0}
# Coarsest size allowed when the node limit forces P1 coarsening.
MINIMUM_ELEMENTS_PER_WAVELENGTH = 6.0
THROAT_SIZE_FRACTION = 0.25
THROAT_AREA_TOLERANCE = 0.01
DEFAULT_MAX_MESH_NODES = 100_000
_DISPERSION_CONSTANT = {1: 1 / 24, 2: 1 / 1440}


class MeshRejection(RuntimeError):
    """A candidate cannot be meshed within the qualification rules."""
    reason = "mesh_rejected"


class MeshTooLargeError(MeshRejection):
    reason = "mesh_node_limit_exceeded"


class ThroatAreaGuardError(MeshRejection):
    reason = "throat_area_guard_failed"


def _check_degree(degree):
    if degree not in ELEMENTS_PER_WAVELENGTH:
        raise ValueError("Supported finite-element degrees are 1 and 2")


def wavelength_mesh_size(c, f_max, degree):
    """Target element size c / (n f_max) for this degree."""
    _check_degree(degree)
    if not (np.isfinite(c) and np.isfinite(f_max) and c > 0 and f_max > 0):
        raise ValueError("Sound speed and maximum frequency must be positive")
    return c / (ELEMENTS_PER_WAVELENGTH[degree] * f_max)


def target_mesh_size(user_size, c, f_max, degree):
    """Finer of the user size and the degree's wavelength rule."""
    if not (np.isfinite(user_size) and user_size > 0):
        raise ValueError("Mesh size must be positive")
    return min(user_size, wavelength_mesh_size(c, f_max, degree))


def coarsest_mesh_size(user_size, c, f_max, degree):
    """Coarsest size the node limit may force: lambda/6, never above the user size."""
    target = target_mesh_size(user_size, c, f_max, degree)
    return max(target, min(user_size, c / (MINIMUM_ELEMENTS_PER_WAVELENGTH * f_max)))


def dispersion_phase_error(degree, elements_per_wavelength):
    """Leading-order relative phase-speed error per wavelength."""
    _check_degree(degree)
    kh = 2 * math.pi / elements_per_wavelength
    return _DISPERSION_CONSTANT[degree] * kh ** (2 * degree)


def max_mesh_nodes(environ=None):
    environ = os.environ if environ is None else environ
    value = int(environ.get("HORN_MAX_MESH_NODES", DEFAULT_MAX_MESH_NODES))
    if value < 1:
        raise ValueError("HORN_MAX_MESH_NODES must be positive")
    return value


def throat_mesh_size(throat_radius, bulk_size, fraction=THROAT_SIZE_FRACTION):
    if not (throat_radius > 0 and bulk_size > 0):
        raise ValueError("Throat radius and bulk size must be positive")
    return min(bulk_size, fraction * throat_radius)


def set_throat_refinement(gmsh, inlet_surfaces, throat_radius, bulk_size, fraction=THROAT_SIZE_FRACTION):
    """Install the throat Distance/Threshold background field. Returns the throat size."""
    size = throat_mesh_size(throat_radius, bulk_size, fraction)
    field = gmsh.model.mesh.field
    for tag in field.list():
        field.remove(tag)
    distance = field.add("Distance")
    field.setNumbers(distance, "SurfacesList", list(inlet_surfaces))
    field.setNumber(distance, "Sampling", 40)
    threshold = field.add("Threshold")
    field.setNumber(threshold, "InField", distance)
    field.setNumber(threshold, "SizeMin", size)
    field.setNumber(threshold, "SizeMax", bulk_size)
    field.setNumber(threshold, "DistMin", throat_radius)
    field.setNumber(threshold, "DistMax", 3 * throat_radius)
    field.setAsBackgroundMesh(threshold)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMax", bulk_size)
    return size


def surface_mesh_area(gmsh, surfaces):
    """Area of the straight-sided triangles on the given surfaces."""
    tags, coords, _ = gmsh.model.mesh.getNodes()
    lookup = np.zeros(int(tags.max()) + 1, dtype=np.int64)
    lookup[tags.astype(np.int64)] = np.arange(len(tags))
    xyz = coords.reshape(-1, 3)
    area = 0.0
    for surface in surfaces:
        types, _, nodes = gmsh.model.mesh.getElements(2, surface)
        for element_type, connectivity in zip(types, nodes):
            name, _, _, per, _, primary = gmsh.model.mesh.getElementProperties(element_type)
            if primary != 3:
                raise ValueError(f"Only triangular surface meshes are supported, got {name}")
            triangles = lookup[np.asarray(connectivity, dtype=np.int64)].reshape(-1, per)[:, :3]
            a, b, c = (xyz[triangles[:, i]] for i in range(3))
            area += 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum()
    return float(area)


def node_count(gmsh):
    return len(gmsh.model.mesh.getNodes()[0])


def generate_qualified_mesh(gmsh, inlet_surfaces, cad_inlet_area, *, user_size, c, f_max,
                            degree, max_nodes=None):
    """Mesh the current gmsh model under the sizing, throat and node rules.

    Returns a dict describing the mesh actually generated. Raises
    MeshTooLargeError or ThroatAreaGuardError for a per-candidate rejection.
    """
    _check_degree(degree)
    if not (cad_inlet_area > 0):
        raise ValueError("CAD throat area must be positive")
    max_nodes = max_mesh_nodes() if max_nodes is None else int(max_nodes)
    throat_radius = math.sqrt(cad_inlet_area / math.pi)
    target = target_mesh_size(user_size, c, f_max, degree)
    coarsest = coarsest_mesh_size(user_size, c, f_max, degree)
    size, fraction = target, THROAT_SIZE_FRACTION
    area_retried = False
    for _ in range(8):
        gmsh.model.mesh.clear()
        throat_size = set_throat_refinement(gmsh, inlet_surfaces, throat_radius, size, fraction)
        gmsh.model.mesh.generate(3)
        ratio = surface_mesh_area(gmsh, inlet_surfaces) / cad_inlet_area
        if abs(ratio - 1) > THROAT_AREA_TOLERANCE:
            if area_retried:
                raise ThroatAreaGuardError(
                    f"Faceted throat area is {100 * (ratio - 1):+.2f} % from CAD "
                    f"(limit {100 * THROAT_AREA_TOLERANCE:g} %) with throat size {throat_size:.4g} m")
            area_retried, fraction = True, fraction / 2
            continue
        if degree == 2:
            gmsh.model.mesh.setOrder(2)
        nodes = node_count(gmsh)
        if nodes <= max_nodes:
            return dict(mesh_size_m=size, target_mesh_size_m=target, throat_mesh_size_m=throat_size,
                        mesh_nodes=nodes, max_mesh_nodes=max_nodes,
                        elements_per_wavelength=c / (f_max * size),
                        target_elements_per_wavelength=c / (f_max * target),
                        mesh_coarsened_for_node_limit=bool(size > target * (1 + 1e-9)),
                        throat_area_ratio=ratio)
        if size >= coarsest * (1 - 1e-9):
            raise MeshTooLargeError(
                f"Mesh needs {nodes} nodes at {c / (f_max * size):.1f} elements per wavelength, "
                f"above the limit of {max_nodes} (HORN_MAX_MESH_NODES)")
        size = min(coarsest, size * 1.05 * (nodes / max_nodes) ** (1 / 3))
    raise MeshTooLargeError(f"Mesh did not fit {max_nodes} nodes after repeated coarsening")


def rejection_reason(error):
    """Ranking rejection reason for a per-candidate failure, or None if it is a real fault."""
    from horn_solver.modal_geometry import ModalGeometryError
    if isinstance(error, MeshRejection):
        return error.reason
    if isinstance(error, ModalGeometryError):
        return "cad_geometry_rejected"
    return None


def write_rejection(path, error, **context):
    """Record a per-candidate rejection so the workflow can continue without it."""
    import json
    from pathlib import Path
    reason = rejection_reason(error)
    if reason is None:
        raise ValueError("Only per-candidate mesh or CAD rejections can be recorded")
    record = {"status": "rejected", "rejection_reasons": [reason], "rejection_detail": str(error), **context}
    Path(path).write_text(json.dumps(record, indent=2))
    return record
