import math

import gmsh
import numpy as np
from pathlib import Path
from typing import Callable
import argparse


def _loft_horn_profile(
    radius_func: Callable[[float], float],
    length: float,
    num_sections: int,
    output_file: Path,
    model_name: str,
) -> Path:
    """Shared helper that lofts through circular cross-sections along the z-axis.

    Args:
        radius_func: Function mapping z-position [0, length] to radius.
        length: Horn length along z-axis (m).
        num_sections: Number of intermediate cross-sections (including throat and mouth).
        output_file: Path for the output STEP file.
        model_name: Name for the gmsh model.

    Returns:
        Path to the generated STEP file.
    """
    gmsh.initialize()
    gmsh.model.add(model_name)
    gmsh.option.setNumber("General.Terminal", 1)

    z_positions = np.linspace(0, length, num_sections)
    curve_loops = []

    for z in z_positions:
        r = radius_func(z)
        wire = gmsh.model.occ.addCircle(0, 0, z, r)
        loop = gmsh.model.occ.addCurveLoop([wire])
        curve_loops.append(loop)

    # Create end-cap plane surfaces (throat and mouth)
    gmsh.model.occ.addPlaneSurface([curve_loops[0]])
    gmsh.model.occ.addPlaneSurface([curve_loops[-1]])

    # Loft through all sections
    gmsh.model.occ.addThruSections(curve_loops)
    gmsh.model.occ.synchronize()

    gmsh.write(str(output_file))
    gmsh.finalize()

    print(str(output_file))
    return output_file


def create_conical_horn(
    throat_radius: float,
    mouth_radius: float,
    length: float,
    output_file: Path,
    num_sections: int = 2,
) -> Path:
    """Generate a conical horn STEP file.

    A conical horn has a linear radius profile: r(z) = r_throat + (r_mouth - r_throat) * z / L.
    Only 2 sections are needed (throat + mouth), but more can be used for consistency.
    """

    def radius_func(z: float) -> float:
        return throat_radius + (mouth_radius - throat_radius) * z / length

    return _loft_horn_profile(radius_func, length, num_sections, output_file, "conical_horn")


def create_exponential_horn(
    throat_radius: float,
    mouth_radius: float,
    length: float,
    output_file: Path,
    num_sections: int = 20,
) -> Path:
    """Generate an exponential horn STEP file.

    Radius profile: r(z) = r_throat * exp(m * z / L)
    where m = ln(r_mouth / r_throat).
    """
    m = np.log(mouth_radius / throat_radius)

    def radius_func(z: float) -> float:
        return throat_radius * np.exp(m * z / length)

    return _loft_horn_profile(radius_func, length, num_sections, output_file, "exponential_horn")


def create_hyperbolic_horn(
    throat_radius: float,
    mouth_radius: float,
    length: float,
    output_file: Path,
    num_sections: int = 20,
) -> Path:
    """Generate a hyperbolic (hypex) horn STEP file.

    Radius profile: r(z) = r_throat * cosh(m * z / L)
    where m = acosh(r_mouth / r_throat).
    """
    m = np.arccosh(mouth_radius / throat_radius)

    def radius_func(z: float) -> float:
        return throat_radius * np.cosh(m * z / length)

    return _loft_horn_profile(radius_func, length, num_sections, output_file, "hyperbolic_horn")


def create_tractrix_horn(
    throat_radius: float,
    mouth_radius: float,
    length: float,
    output_file: Path,
    num_sections: int = 20,
) -> Path:
    """Generate a tractrix horn STEP file.

    Scaled tractrix parametric curve: rapid initial expansion, decelerating toward mouth.
    """
    t = np.linspace(np.pi - 1e-6, np.pi / 2, 500)
    y, x = np.sin(t), np.log(np.tan(t / 2)) + np.cos(t)
    x -= x[0]
    x_n, y_n = x / x[-1], y / y[-1]

    def radius_func(z: float) -> float:
        return throat_radius + (mouth_radius - throat_radius) * float(
            np.interp(z, x_n * length, y_n)
        )

    return _loft_horn_profile(radius_func, length, num_sections, output_file, "tractrix_horn")


def create_os_horn(
    throat_radius: float,
    mouth_radius: float,
    length: float,
    output_file: Path,
    num_sections: int = 20,
) -> Path:
    """Generate an oblate spheroidal (OS / Geddes) horn STEP file.

    Simple closed-form: r(z) = sqrt(r_t^2 + (z * tan(theta))^2).
    """
    theta = math.atan2(math.sqrt(mouth_radius**2 - throat_radius**2), length)

    def radius_func(z: float) -> float:
        return math.sqrt(throat_radius**2 + (z * math.tan(theta)) ** 2)

    return _loft_horn_profile(radius_func, length, num_sections, output_file, "os_horn")


def create_lecleach_horn(
    throat_radius: float,
    mouth_radius: float,
    length: float,
    output_file: Path,
    num_sections: int = 20,
) -> Path:
    """Generate a Le Cléac'h horn STEP file.

    Tractrix curve with constant a = mouth_radius, clipped where radius >= throat_radius.
    Gentle expansion at throat, aggressive flare at mouth.
    """
    t = np.linspace(np.pi - 1e-6, np.pi / 2, 500)
    y, x = np.sin(t), np.log(np.tan(t / 2)) + np.cos(t)
    x -= x[0]
    idx = np.searchsorted(y, throat_radius / mouth_radius)
    x_c, y_c = x[idx:] - x[idx], y[idx:]

    def radius_func(z: float) -> float:
        return float(np.interp(z, x_c / x_c[-1] * length, y_c / y_c[-1] * mouth_radius))

    return _loft_horn_profile(radius_func, length, num_sections, output_file, "lecleach_horn")


def create_cd_horn(
    throat_radius: float,
    mouth_radius: float,
    length: float,
    output_file: Path,
    num_sections: int = 30,
) -> Path:
    """Generate a constant directivity (CD) horn STEP file.

    Compound: exponential throat (30% of length) transitioning to conical body.
    """
    frac = 0.3
    z_t = frac * length
    r_trans = throat_radius * (mouth_radius / throat_radius) ** frac

    def radius_func(z: float) -> float:
        if z <= z_t:
            return throat_radius * math.exp(math.log(r_trans / throat_radius) * z / z_t)
        else:
            return r_trans + (mouth_radius - r_trans) * (z - z_t) / (length - z_t)

    return _loft_horn_profile(radius_func, length, num_sections, output_file, "cd_horn")


_PROFILE_DISPATCH = {
    "conical": create_conical_horn,
    "exponential": create_exponential_horn,
    "hyperbolic": create_hyperbolic_horn,
    "tractrix": create_tractrix_horn,
    "os": create_os_horn,
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
) -> Path:
    """Dispatch to the appropriate horn profile generator.

    Args:
        profile: One of "conical", "exponential", "hyperbolic".
        throat_radius: Throat radius (m).
        mouth_radius: Mouth radius (m).
        length: Horn length (m).
        output_file: Output STEP file path.
        num_sections: Number of cross-sections for lofting.

    Returns:
        Path to the generated STEP file.
    """
    if profile not in _PROFILE_DISPATCH:
        raise ValueError(f"Unknown profile '{profile}'. Choose from: {list(_PROFILE_DISPATCH.keys())}")
    return _PROFILE_DISPATCH[profile](
        throat_radius=throat_radius,
        mouth_radius=mouth_radius,
        length=length,
        output_file=output_file,
        num_sections=num_sections,
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
        help="Horn flare profile (default: conical).",
    )
    parser.add_argument(
        "--num-sections",
        type=int,
        default=20,
        help="Number of cross-sections for lofting (default: 20).",
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
    )


if __name__ == "__main__":
    main()
