"""The lofted STEP wall must follow the profile radius within 1% at the pipeline default.

main.nf and refine.py ask for 20 loft sections. The old CD profile, with a
slope kink at the join, bulged 23% at the throat at that setting.
"""
import gmsh
import numpy as np
import pytest

from horn_core.profiles import PROFILE_NAMES, get_radius_func
from horn_geometry.generator import create_horn

PIPELINE_DEFAULT_SECTIONS = 20
GEOMETRIES = [(0.0254, 0.182, 0.4), (0.0125, 0.2, 0.5), (0.032498, 0.088709, 0.117906)]


def _sampled_wall_error(step, radius_func, samples=61):
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
            u, v = np.meshgrid(np.linspace(lower[0], upper[0], samples),
                               np.linspace(lower[1], upper[1], samples))
            xyz = np.asarray(gmsh.model.getValue(dim, tag, np.column_stack([u.ravel(), v.ravel()]).ravel()))
            xyz = xyz.reshape(-1, 3)
            radius = np.hypot(xyz[:, 0], xyz[:, 1])
            target = np.array([radius_func(float(z)) for z in xyz[:, 2]])
            worst = max(worst, float(np.max(np.abs(radius - target) / target)))
        assert walls >= 1
        return worst
    finally:
        gmsh.finalize()


@pytest.mark.parametrize("geometry", GEOMETRIES, ids=lambda g: "x".join(f"{v*1000:.0f}" for v in g))
@pytest.mark.parametrize("profile", PROFILE_NAMES)
def test_lofted_wall_within_one_percent(tmp_path, profile, geometry):
    r_t, r_m, length = geometry
    func = get_radius_func(profile, r_t, r_m, length)
    step = create_horn(profile, r_t, r_m, length, tmp_path / "horn.step",
                       num_sections=PIPELINE_DEFAULT_SECTIONS)
    assert _sampled_wall_error(step, func) <= 0.01
