"""The STEP wall must follow the profile radius within 1% everywhere along the horn.

main.nf and refine.py ask for 20 profile points. The old CD profile, with a
slope kink at the join, bulged 23% at the throat at that setting. A lofted
wall checked on a coarse 41-point grid was accepted while it deviated 2-3%
between the samples (cd 17.5/250/600 mm, os 12.5/300/100 mm), so these
tests re-import the written STEP and sample it densely along the profile.
"""
import gmsh
import numpy as np
import pytest

from horn_core.profiles import PROFILE_NAMES, get_radius_func, profile_geometry_error
from horn_geometry.generator import create_horn

PIPELINE_DEFAULT_SECTIONS = 20
GEOMETRIES = [
    (0.0254, 0.182, 0.4),
    (0.0125, 0.2, 0.5),
    (0.032498, 0.088709, 0.117906),
    # Accepted by the old coarse check while 2.3% and 3.1% off between samples.
    (0.0175, 0.25, 0.6),
    (0.0125, 0.3, 0.1),
    (0.0125, 0.3, 1.0),
]
DENSE_SAMPLES = 4001


def _dense_wall_error(step, radius_func):
    """Re-import the STEP in metres; sample the wall at 4001 profile x 12 angle points."""
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setString("Geometry.OCCTargetUnit", "M")
        gmsh.model.occ.importShapes(str(step))
        gmsh.model.occ.synchronize()
        worst, walls = 0.0, 0
        for dim, tag in gmsh.model.getEntities(2):
            if gmsh.model.getType(dim, tag) == "Plane":
                continue
            walls += 1
            lower, upper = gmsh.model.getParametrizationBounds(dim, tag)
            u, v = np.meshgrid(np.linspace(lower[0], upper[0], 12),
                               np.linspace(lower[1], upper[1], DENSE_SAMPLES))
            xyz = np.asarray(gmsh.model.getValue(dim, tag, np.column_stack([u.ravel(), v.ravel()]).ravel()))
            xyz = xyz.reshape(-1, 3)
            radius = np.hypot(xyz[:, 0], xyz[:, 1])
            target = np.array([radius_func(float(z)) for z in xyz[:, 2]])
            worst = max(worst, float(np.max(np.abs(radius - target) / target)))
        assert walls >= 1
        return worst
    finally:
        gmsh.finalize()


@pytest.mark.parametrize("geometry", GEOMETRIES, ids=lambda g: "x".join(f"{v*1000:.1f}" for v in g))
@pytest.mark.parametrize("profile", PROFILE_NAMES)
def test_wall_within_one_percent_everywhere(tmp_path, profile, geometry):
    r_t, r_m, length = geometry
    if profile_geometry_error(profile, r_t, r_m, length):
        pytest.skip("tractrix shorter than its natural length")
    func = get_radius_func(profile, r_t, r_m, length)
    step = create_horn(profile, r_t, r_m, length, tmp_path / "horn.step",
                       num_sections=PIPELINE_DEFAULT_SECTIONS)
    assert _dense_wall_error(step, func) <= 0.01


def test_wall_error_falls_as_points_are_added(tmp_path):
    """More profile points must not make the wall worse (the loft oscillated: 0.01% at 80, 11.6% at 160)."""
    r_t, r_m, length = 0.0125, 0.3, 1.0
    func = get_radius_func("cd", r_t, r_m, length)
    errors = []
    for n in (20, 40, 80, 160):
        step = create_horn("cd", r_t, r_m, length, tmp_path / f"cd{n}.step", num_sections=n)
        errors.append(_dense_wall_error(step, func))
    assert all(later <= earlier * 1.05 + 1e-6 for earlier, later in zip(errors, errors[1:])), errors
    assert errors[-1] < 0.001


def test_written_step_is_one_volume_with_two_port_disks(tmp_path):
    step = create_horn("os", 0.0125, 0.3, 0.1, tmp_path / "horn.step")
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setString("Geometry.OCCTargetUnit", "M")
        gmsh.model.occ.importShapes(str(step))
        gmsh.model.occ.synchronize()
        assert len(gmsh.model.getEntities(3)) == 1
        planes = [t for d, t in gmsh.model.getEntities(2) if gmsh.model.getType(d, t) == "Plane"]
        z = sorted(round(gmsh.model.occ.getCenterOfMass(2, t)[2], 9) for t in planes)
        assert z == [0.0, 0.1]
        areas = sorted(gmsh.model.occ.getMass(2, t) for t in planes)
        assert areas == pytest.approx([np.pi * 0.0125**2, np.pi * 0.3**2], rel=1e-6)
    finally:
        gmsh.finalize()
