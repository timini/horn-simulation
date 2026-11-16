"""Generate horn STEP files by revolving a flare profile about the horn axis.

All profile mathematics lives in horn_core.profiles. The STEP files are
written in true size: the model is built in metres, gmsh's OCC target
unit is set to metres, and OCC writes millimetre coordinates under a
millimetre unit declaration.
"""
import argparse
from pathlib import Path
from typing import Callable, Optional

import gmsh
import numpy as np

from horn_core.cad_units import set_step_units
from horn_core.profiles import canonical_profile, get_radius_func


#: Largest permitted radial deviation of the STEP wall from the profile,
#: as a fraction of the local radius.
MAX_LOFT_RADIUS_ERROR = 0.01
#: Upper limit on profile points when refining to meet MAX_LOFT_RADIUS_ERROR.
MAX_LOFT_SECTIONS = 320
#: Minimum wall samples along the profile, and per profile span, in the check.
_MIN_WALL_SAMPLES = 2001
_SAMPLES_PER_SPAN = 20


class LoftAccuracyError(RuntimeError):
    """The STEP wall cannot follow the profile within MAX_LOFT_RADIUS_ERROR."""


def _profile_z_positions(
    radius_func: Callable[[float], float],
    length: float,
    num_points: int,
) -> np.ndarray:
    """Place profile points by arc length, packed where the profile bends.

    The density along the (z, r) curve is proportional to
    1 + sqrt(curvature * total arc length), which is close to the
    optimal spacing for spline interpolation and reduces to uniform arc
    length for a straight cone.
    """
    z_fine = np.linspace(0.0, length, 20001)
    r_fine = np.array([radius_func(z) for z in z_fine])
    ds = np.hypot(np.diff(z_fine), np.diff(r_fine))
    total = ds.sum()
    angle = np.arctan2(np.diff(r_fine), np.diff(z_fine))
    curvature = np.abs(np.gradient(angle)) / np.maximum(ds, 1e-300)
    density = ds * (1.0 + np.sqrt(curvature * total))
    cumulative = np.concatenate(([0.0], np.cumsum(density)))
    z_positions = np.interp(np.linspace(0.0, cumulative[-1], num_points), cumulative, z_fine)
    z_positions[0] = 0.0
    z_positions[-1] = length
    return z_positions


def wall_radius_error(radius_func: Callable[[float], float], num_points: int = 2) -> float:
    """Maximum relative radius error of the horn wall in the current gmsh model.

    Samples every non-planar surface densely along the profile (at least
    _MIN_WALL_SAMPLES points, and _SAMPLES_PER_SPAN per span between the
    num_points profile points) and at eight angles, and compares the
    distance from the z axis with radius_func(z) at the same z.
    """
    along = max(_MIN_WALL_SAMPLES, _SAMPLES_PER_SPAN * max(1, num_points - 1) + 1)
    worst = 0.0
    for dim, tag in gmsh.model.getEntities(2):
        if gmsh.model.getType(dim, tag) == "Plane":
            continue  # throat or mouth disk
        lower, upper = gmsh.model.getParametrizationBounds(dim, tag)
        # v runs along the profile on both a surface of revolution and a
        # lofted wall; u runs around the axis.
        grids = [np.linspace(lower[0], upper[0], 8, endpoint=False),
                 np.linspace(lower[1], upper[1], along)]
        uu, vv = np.meshgrid(grids[0], grids[1])
        params = np.column_stack([uu.ravel(), vv.ravel()]).ravel()
        points = np.asarray(gmsh.model.getValue(dim, tag, params)).reshape(-1, 3)
        radius = np.hypot(points[:, 0], points[:, 1])
        target = np.array([radius_func(float(z)) for z in points[:, 2]])
        worst = max(worst, float(np.max(np.abs(radius - target) / target)))
    return worst


def _build_wall(radius_func, length, num_points, model_name):
    """Revolve the interpolated (r, z) profile about the z axis into a solid.

    The result is one volume bounded by a planar throat disk, a planar
    mouth disk and a single surface of revolution. The revolution is exact,
    so the only geometric error is the profile spline's, which falls
    steadily as points are added.
    """
    occ = gmsh.model.occ
    gmsh.model.add(model_name)
    z_positions = _profile_z_positions(radius_func, length, num_points)
    points = [occ.addPoint(radius_func(z), 0.0, z) for z in z_positions]
    axis_throat = occ.addPoint(0.0, 0.0, 0.0)
    axis_mouth = occ.addPoint(0.0, 0.0, length)
    wall = occ.addSpline(points) if len(points) > 2 else occ.addLine(points[0], points[1])
    loop = occ.addCurveLoop([
        occ.addLine(axis_throat, points[0]),
        wall,
        occ.addLine(points[-1], axis_mouth),
        occ.addLine(axis_mouth, axis_throat),
    ])
    half_section = occ.addPlaneSurface([loop])
    occ.revolve([(2, half_section)], 0, 0, 0, 0, 0, 1, 2 * np.pi)
    occ.remove([(2, half_section)], recursive=True)
    occ.synchronize()


def _loft_horn_profile(
    radius_func: Callable[[float], float],
    length: float,
    num_sections: int,
    output_file: Path,
    model_name: str,
) -> Path:
    """Revolve the horn profile about the z-axis and write a STEP file.

    The wall is a surface of revolution of a spline through num_sections
    profile points. If it deviates from radius_func by more than
    MAX_LOFT_RADIUS_ERROR of the local radius anywhere along the profile,
    the point count is doubled until it does not (up to MAX_LOFT_SECTIONS).

    Args:
        radius_func: Function mapping z-position [0, length] to radius.
        length: Horn length along z-axis (m).
        num_sections: Minimum number of profile points (including throat and mouth).
        output_file: Path for the output STEP file.
        model_name: Name for the gmsh model.

    Returns:
        Path to the generated STEP file.

    Raises:
        LoftAccuracyError: the wall is still outside the limit at MAX_LOFT_SECTIONS.
    """
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1)
        set_step_units(gmsh)
        sections = max(2, int(num_sections))
        while True:
            _build_wall(radius_func, length, sections, model_name)
            error = wall_radius_error(radius_func, sections)
            if error <= MAX_LOFT_RADIUS_ERROR:
                break
            if sections >= MAX_LOFT_SECTIONS:
                raise LoftAccuracyError(
                    f"{model_name}: horn wall deviates {100 * error:.2f}% from the profile "
                    f"at {sections} profile points (limit {100 * MAX_LOFT_RADIUS_ERROR:.0f}%)")
            gmsh.model.remove()
            sections = min(2 * sections, MAX_LOFT_SECTIONS)
        if sections != num_sections:
            print(f"{model_name}: used {sections} profile points (requested {num_sections}); "
                  f"wall error {100 * error:.2f}%")
        gmsh.write(str(output_file))
    finally:
        gmsh.finalize()

    print(str(output_file))
    return output_file


def _create_profile_horn(profile, throat_radius, mouth_radius, length, output_file,
                         num_sections, **profile_kwargs):
    radius_func = get_radius_func(profile, throat_radius, mouth_radius, length, **profile_kwargs)
    return _loft_horn_profile(radius_func, length, num_sections, output_file,
                              f"{canonical_profile(profile)}_horn")


def create_conical_horn(throat_radius, mouth_radius, length, output_file, num_sections: int = 2) -> Path:
    """Generate a conical horn STEP file (linear radius profile)."""
    return _create_profile_horn("conical", throat_radius, mouth_radius, length, output_file, num_sections)


def create_exponential_horn(throat_radius, mouth_radius, length, output_file, num_sections: int = 20) -> Path:
    """Generate an exponential horn STEP file: r = r_t exp(ln(r_m/r_t) z / L)."""
    return _create_profile_horn("exponential", throat_radius, mouth_radius, length, output_file, num_sections)


def create_hyperbolic_horn(throat_radius, mouth_radius, length, output_file, num_sections: int = 20,
                           hypex_t: Optional[float] = None) -> Path:
    """Generate a hypex horn STEP file; the default T = 0 is the catenoidal horn."""
    kwargs = {} if hypex_t is None else {"hypex_t": hypex_t}
    return _create_profile_horn("hyperbolic", throat_radius, mouth_radius, length, output_file,
                                num_sections, **kwargs)


def create_tractrix_horn(throat_radius, mouth_radius, length, output_file, num_sections: int = 20) -> Path:
    """Generate a true tractrix horn STEP file (see horn_core.profiles.tractrix_radius)."""
    return _create_profile_horn("tractrix", throat_radius, mouth_radius, length, output_file, num_sections)


def create_os_horn(throat_radius, mouth_radius, length, output_file, num_sections: int = 20) -> Path:
    """Generate an oblate spheroidal (OS / Geddes) horn STEP file."""
    return _create_profile_horn("os", throat_radius, mouth_radius, length, output_file, num_sections)


def create_scaled_tractrix_horn(throat_radius, mouth_radius, length, output_file,
                                num_sections: int = 20) -> Path:
    """Generate an axially scaled tractrix horn STEP file (formerly "lecleach")."""
    return _create_profile_horn("scaled_tractrix", throat_radius, mouth_radius, length, output_file,
                                num_sections)


def create_lecleach_horn(throat_radius, mouth_radius, length, output_file, num_sections: int = 20) -> Path:
    """Deprecated alias: the old "lecleach" profile is the scaled tractrix."""
    return _create_profile_horn("lecleach", throat_radius, mouth_radius, length, output_file, num_sections)


def create_cd_horn(throat_radius, mouth_radius, length, output_file, num_sections: int = 30) -> Path:
    """Generate an exponential-throat, tangent-cone ("cd") horn STEP file."""
    return _create_profile_horn("cd", throat_radius, mouth_radius, length, output_file, num_sections)


_PROFILE_DISPATCH = {
    "conical": create_conical_horn,
    "exponential": create_exponential_horn,
    "hyperbolic": create_hyperbolic_horn,
    "tractrix": create_tractrix_horn,
    "os": create_os_horn,
    "scaled_tractrix": create_scaled_tractrix_horn,
    "lecleach": create_lecleach_horn,
    "cd": create_cd_horn,
}


def create_horn(
    profile: str,
    throat_radius: float,
    mouth_radius: float,
    length: float,
    output_file: Path,
    num_sections: int = 20,
    hypex_t: Optional[float] = None,
) -> Path:
    """Dispatch to the appropriate horn profile generator.

    Args:
        profile: A name in horn_core.profiles.PROFILE_NAMES ("lecleach" is a
            deprecated alias of "scaled_tractrix").
        throat_radius: Throat radius (m).
        mouth_radius: Mouth radius (m).
        length: Horn length (m).
        output_file: Output STEP file path.
        num_sections: Minimum number of profile points in the revolved wall.
        hypex_t: Hypex T for the "hyperbolic" profile (default 0, catenoidal).

    Returns:
        Path to the generated STEP file.

    Raises:
        horn_core.profiles.ProfileGeometryError: the profile cannot meet the
            requested dimensions (a tractrix shorter than its natural length).
    """
    if profile not in _PROFILE_DISPATCH:
        raise ValueError(f"Unknown profile '{profile}'. Choose from: {list(_PROFILE_DISPATCH.keys())}")
    kwargs = {}
    if hypex_t is not None:
        if profile != "hyperbolic":
            raise ValueError("hypex_t applies only to the 'hyperbolic' profile")
        kwargs["hypex_t"] = hypex_t
    return _PROFILE_DISPATCH[profile](
        throat_radius=throat_radius,
        mouth_radius=mouth_radius,
        length=length,
        output_file=output_file,
        num_sections=num_sections,
        **kwargs,
    )


def main():
    """Command-line interface for the horn geometry generator."""
    parser = argparse.ArgumentParser(description="Generate a horn STEP file.")
    parser.add_argument("--throat-radius", type=float, required=True)
    parser.add_argument("--mouth-radius", type=float, required=True)
    parser.add_argument("--length", type=float, required=True)
    parser.add_argument("--output-file", type=str, required=True)
    parser.add_argument(
        "--profile",
        type=str,
        choices=list(_PROFILE_DISPATCH.keys()),
        default="conical",
        help="Horn flare profile (default: conical). 'lecleach' is a deprecated alias of "
             "'scaled_tractrix'.",
    )
    parser.add_argument(
        "--num-sections",
        type=int,
        default=20,
        help="Minimum number of profile points in the revolved wall (default: 20). More are "
             "used when needed to keep the wall within 1%% of the profile radius.",
    )
    parser.add_argument(
        "--hypex-t",
        type=float,
        default=None,
        help="Hypex T for the hyperbolic profile (default 0: catenoidal; 1: exponential).",
    )
    args = parser.parse_args()

    output_path = Path(args.output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    create_horn(
        profile=args.profile,
        throat_radius=args.throat_radius,
        mouth_radius=args.mouth_radius,
        length=args.length,
        output_file=output_path,
        num_sections=args.num_sections,
        hypex_t=args.hypex_t,
    )


if __name__ == "__main__":
    main()
