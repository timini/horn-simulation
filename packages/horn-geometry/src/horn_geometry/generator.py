"""Generate horn STEP files by lofting circular sections through a flare profile.

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


def _adaptive_z_positions(
    radius_func: Callable[[float], float],
    length: float,
    num_sections: int,
) -> np.ndarray:
    """Compute z-positions that concentrate sections where the radius changes fast.

    Uses a blended metric: 50 % uniform arc-length + 50 % radius-change.  This
    ensures adequate section density everywhere (no huge gaps in flat regions)
    while still packing extra sections into steep regions (e.g. the mouth of a
    tractrix).  Degenerates to uniform spacing for constant-gradient profiles.
    """
    n_samples = 5000
    z_fine = np.linspace(0, length, n_samples)
    r_fine = np.array([radius_func(z) for z in z_fine])

    dr = np.abs(np.diff(r_fine))
    dz = np.diff(z_fine)

    # Blend: uniform component (dz) + curvature component (dr)
    total_dr = dr.sum()
    total_dz = dz.sum()

    if total_dr == 0:
        return np.linspace(0, length, num_sections)

    # Normalise both to the same scale, then blend 50/50
    metric = 0.5 * (dz / total_dz) + 0.5 * (dr / total_dr)
    cum_metric = np.concatenate(([0.0], np.cumsum(metric)))

    targets = np.linspace(0, cum_metric[-1], num_sections)
    z_positions = np.interp(targets, cum_metric, z_fine)

    # Guarantee exact endpoints
    z_positions[0] = 0.0
    z_positions[-1] = length
    return z_positions


#: Largest permitted radial deviation of the lofted wall from the profile,
#: as a fraction of the local radius.
MAX_LOFT_RADIUS_ERROR = 0.01
#: Upper limit on loft sections when refining to meet MAX_LOFT_RADIUS_ERROR.
MAX_LOFT_SECTIONS = 320


def wall_radius_error(radius_func: Callable[[float], float], samples: int = 41) -> float:
    """Maximum relative radius error of the lofted wall in the current gmsh model.

    Samples every non-planar surface on a (u, v) grid and compares the
    distance from the z axis with radius_func(z) at the same z.
    """
    worst = 0.0
    for dim, tag in gmsh.model.getEntities(2):
        if gmsh.model.getType(dim, tag) == "Plane":
            continue  # throat or mouth disk
        lower, upper = gmsh.model.getParametrizationBounds(dim, tag)
        uu, vv = np.meshgrid(np.linspace(lower[0], upper[0], samples),
                             np.linspace(lower[1], upper[1], samples))
        params = np.column_stack([uu.ravel(), vv.ravel()]).ravel()
        points = np.asarray(gmsh.model.getValue(dim, tag, params)).reshape(-1, 3)
        radius = np.hypot(points[:, 0], points[:, 1])
        target = np.array([radius_func(float(z)) for z in points[:, 2]])
        worst = max(worst, float(np.max(np.abs(radius - target) / target)))
    return worst


def _build_loft(radius_func, length, num_sections, model_name):
    gmsh.model.add(model_name)
    z_positions = _adaptive_z_positions(radius_func, length, num_sections)
    curve_loops = []
    for z in z_positions:
        wire = gmsh.model.occ.addCircle(0, 0, z, radius_func(z))
        curve_loops.append(gmsh.model.occ.addCurveLoop([wire]))
    # End-cap plane surfaces (throat and mouth), then the lofted wall.
    gmsh.model.occ.addPlaneSurface([curve_loops[0]])
    gmsh.model.occ.addPlaneSurface([curve_loops[-1]])
    gmsh.model.occ.addThruSections(curve_loops)
    gmsh.model.occ.synchronize()


def _loft_horn_profile(
    radius_func: Callable[[float], float],
    length: float,
    num_sections: int,
    output_file: Path,
    model_name: str,
) -> Path:
    """Loft circular cross-sections along the z-axis and write a STEP file.

    If the lofted wall deviates from radius_func by more than
    MAX_LOFT_RADIUS_ERROR of the local radius, the section count is
    doubled until it does not (up to MAX_LOFT_SECTIONS).

    Args:
        radius_func: Function mapping z-position [0, length] to radius.
        length: Horn length along z-axis (m).
        num_sections: Minimum number of cross-sections (including throat and mouth).
        output_file: Path for the output STEP file.
        model_name: Name for the gmsh model.

    Returns:
        Path to the generated STEP file.
    """
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1)
        set_step_units(gmsh)
        sections = max(2, int(num_sections))
        while True:
            _build_loft(radius_func, length, sections, model_name)
            error = wall_radius_error(radius_func)
            if error <= MAX_LOFT_RADIUS_ERROR:
                break
            if sections >= MAX_LOFT_SECTIONS:
                raise RuntimeError(
                    f"{model_name}: lofted wall deviates {100 * error:.2f}% from the profile "
                    f"at {sections} sections (limit {100 * MAX_LOFT_RADIUS_ERROR:.0f}%)")
            gmsh.model.remove()
            sections = min(2 * sections, MAX_LOFT_SECTIONS)
        if sections != num_sections:
            print(f"{model_name}: used {sections} loft sections (requested {num_sections}); "
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
        num_sections: Minimum number of cross-sections for lofting.
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
        help="Minimum number of cross-sections for lofting (default: 20). More are used "
             "when needed to keep the wall within 1%% of the profile radius.",
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
