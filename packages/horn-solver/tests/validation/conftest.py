"""Shared fixtures and helpers for the validation test suite."""

import csv
import json
import pytest
import numpy as np
from pathlib import Path

try:
    from horn_solver.solver import run_simulation_from_step, C0
    from horn_geometry.generator import (
        create_conical_horn,
        create_exponential_horn,
        create_hyperbolic_horn,
    )
    import gmsh

    SOLVER_AVAILABLE = True
except ImportError:
    SOLVER_AVAILABLE = False

REFERENCE_DATA_DIR = Path(__file__).parent / "reference_data"

# Skip all validation tests if dolfinx/solver is not available
pytestmark = pytest.mark.skipif(
    not SOLVER_AVAILABLE, reason="DOLFINx or horn packages not installed"
)


def load_reference(filename: str) -> dict:
    """Load a JSON reference data file."""
    with open(REFERENCE_DATA_DIR / filename) as f:
        return json.load(f)


def assert_spl_within_tolerance(
    computed: np.ndarray,
    reference: np.ndarray | float,
    tolerance_db: float,
    frequencies: np.ndarray,
    label: str = "",
):
    """Compare SPL arrays, print detailed table on failure.

    Parameters
    ----------
    computed : array of computed SPL values (dB)
    reference : array of reference SPL values (dB), or a scalar
    tolerance_db : maximum allowed deviation (dB)
    frequencies : frequency array for diagnostic output
    label : optional label for the comparison
    """
    computed = np.asarray(computed)
    reference = np.atleast_1d(np.asarray(reference, dtype=float))
    if reference.size == 1:
        reference = np.full_like(computed, reference.item())

    diff = np.abs(computed - reference)
    max_diff = np.max(diff)
    mean_diff = np.mean(diff)

    if max_diff > tolerance_db:
        # Build diagnostic table
        lines = [f"\n{'=' * 70}"]
        lines.append(f"SPL VALIDATION FAILURE{f' ({label})' if label else ''}")
        lines.append(f"Tolerance: {tolerance_db:.1f} dB | Max deviation: {max_diff:.2f} dB")
        lines.append(f"{'=' * 70}")
        lines.append(f"{'Freq (Hz)':>10} {'Computed':>10} {'Reference':>10} {'Diff':>8} {'Status':>8}")
        lines.append(f"{'-' * 70}")
        for i in range(len(frequencies)):
            status = "FAIL" if diff[i] > tolerance_db else "ok"
            lines.append(
                f"{frequencies[i]:10.1f} {computed[i]:10.2f} {reference[i]:10.2f} "
                f"{diff[i]:8.2f} {status:>8}"
            )
        lines.append(f"{'=' * 70}")
        detail = "\n".join(lines)
        pytest.fail(
            f"SPL exceeds {tolerance_db} dB tolerance (max diff={max_diff:.2f} dB, "
            f"mean diff={mean_diff:.2f} dB)\n{detail}"
        )


def run_solver_and_get_spl(
    step_file: str | Path,
    freq_range: tuple[float, float],
    num_intervals: int,
    horn_length: float,
    mesh_size: float = 0.01,
    tmp_dir: Path | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Run the solver on a STEP file and return (frequencies, spl) arrays.

    Parameters
    ----------
    step_file : path to the STEP geometry
    freq_range : (min_freq, max_freq) in Hz
    num_intervals : number of frequency points
    horn_length : horn length for boundary tagging
    mesh_size : element size (m), default 5mm
    tmp_dir : directory for output CSV

    Returns
    -------
    (frequencies, spl) : numpy arrays
    """
    import pandas as pd

    if tmp_dir is None:
        import tempfile
        tmp_dir = Path(tempfile.mkdtemp())

    output_file = tmp_dir / "validation_results.csv"

    driver_params = {"length": horn_length}

    run_simulation_from_step(
        step_file=str(step_file),
        freq_range=freq_range,
        num_intervals=num_intervals,
        driver_params=driver_params,
        output_file=str(output_file),
        max_freq_mesh=freq_range[1],
        mesh_size=mesh_size,
    )

    df = pd.read_csv(output_file)
    return df["frequency"].values, df["spl"].values


def _generate_cylinder_step(output_file: Path, radius: float, length: float) -> Path:
    """Generate a cylinder STEP file.

    Tries lofting identical circles first (via create_conical_horn).
    Falls back to gmsh addCylinder if lofting fails.
    """
    try:
        create_conical_horn(
            throat_radius=radius,
            mouth_radius=radius,
            length=length,
            output_file=output_file,
        )
        # Verify the file was created and has content
        if output_file.exists() and output_file.stat().st_size > 0:
            return output_file
    except Exception:
        pass

    # Fallback: use gmsh.model.occ.addCylinder directly
    gmsh.initialize()
    try:
        gmsh.model.add("cylinder")
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, length, radius)
        gmsh.model.occ.synchronize()
        gmsh.write(str(output_file))
    finally:
        gmsh.finalize()

    return output_file


@pytest.fixture
def straight_tube_step(tmp_path):
    """Generate a cylinder STEP file (equal-radii cone) for V1 validation."""
    ref = load_reference("straight_tube_analytical.json")
    geom = ref["geometry"]
    output = tmp_path / "straight_tube.step"
    _generate_cylinder_step(output, geom["throat_radius_m"], geom["length_m"])
    return output, ref


@pytest.fixture
def conical_horn_step(tmp_path):
    """Generate the canonical conical horn STEP for V2 validation."""
    ref = load_reference("conical_horn_webster.json")
    geom = ref["geometry"]
    output = tmp_path / "conical_horn.step"
    create_conical_horn(
        throat_radius=geom["throat_radius_m"],
        mouth_radius=geom["mouth_radius_m"],
        length=geom["length_m"],
        output_file=output,
    )
    return output, ref


# --- Cached solver results (run once per session, share across tests) ---

@pytest.fixture(scope="module")
def straight_tube_results(tmp_path_factory):
    """Solve V1 straight tube once and cache (frequencies, spl, ref)."""
    ref = load_reference("straight_tube_analytical.json")
    geom = ref["geometry"]
    freq_cfg = ref["frequency_range"]
    tmp = tmp_path_factory.mktemp("v1")
    step_file = tmp / "straight_tube.step"
    _generate_cylinder_step(step_file, geom["throat_radius_m"], geom["length_m"])
    frequencies, spl = run_solver_and_get_spl(
        step_file=step_file,
        freq_range=(freq_cfg["min_hz"], freq_cfg["max_hz"]),
        num_intervals=freq_cfg["num_points"],
        horn_length=geom["length_m"],
        tmp_dir=tmp,
    )
    return frequencies, spl, ref


@pytest.fixture(scope="module")
def conical_horn_results(tmp_path_factory):
    """Solve V2 conical horn once and cache (frequencies, spl, ref)."""
    ref = load_reference("conical_horn_webster.json")
    geom = ref["geometry"]
    freq_cfg = ref["frequency_range"]
    tmp = tmp_path_factory.mktemp("v2")
    step_file = tmp / "conical_horn.step"
    create_conical_horn(
        throat_radius=geom["throat_radius_m"],
        mouth_radius=geom["mouth_radius_m"],
        length=geom["length_m"],
        output_file=step_file,
    )
    frequencies, spl = run_solver_and_get_spl(
        step_file=step_file,
        freq_range=(freq_cfg["min_hz"], freq_cfg["max_hz"]),
        num_intervals=freq_cfg["num_points"],
        horn_length=geom["length_m"],
        tmp_dir=tmp,
    )
    return frequencies, spl, ref


# --- CSV reference helpers ---


def load_csv_reference(filename: str) -> tuple[np.ndarray, np.ndarray] | None:
    """Load frequency response reference data from a CSV file.

    Returns (frequencies, spl) arrays, or None if the file has no data rows.
    Header comment lines starting with '#' are skipped.
    """
    csv_path = REFERENCE_DATA_DIR / filename
    if not csv_path.exists():
        return None

    frequencies = []
    spl_values = []

    with open(csv_path) as f:
        reader = csv.DictReader(
            (row for row in f if not row.startswith("#")),
        )
        for row in reader:
            try:
                frequencies.append(float(row["frequency_hz"]))
                spl_values.append(float(row["spl_db"]))
            except (ValueError, KeyError):
                continue

    if not frequencies:
        return None

    return np.array(frequencies), np.array(spl_values)


def has_reference_data(filename: str) -> bool:
    """Check if a reference CSV file has actual data rows."""
    result = load_csv_reference(filename)
    return result is not None


# --- V4/V5 cached solver results ---


@pytest.fixture(scope="module")
def exponential_horn_results(tmp_path_factory):
    """Solve V4 exponential horn once and cache (frequencies, spl, ref)."""
    ref = load_reference("exponential_horn_webster.json")
    geom = ref["geometry"]
    freq_cfg = ref["frequency_range"]
    tmp = tmp_path_factory.mktemp("v4")
    step_file = tmp / "exponential_horn.step"
    create_exponential_horn(
        throat_radius=geom["throat_radius_m"],
        mouth_radius=geom["mouth_radius_m"],
        length=geom["length_m"],
        output_file=step_file,
    )
    frequencies, spl = run_solver_and_get_spl(
        step_file=step_file,
        freq_range=(freq_cfg["min_hz"], freq_cfg["max_hz"]),
        num_intervals=freq_cfg["num_points"],
        horn_length=geom["length_m"],
        mesh_size=0.002,
        tmp_dir=tmp,
    )
    return frequencies, spl, ref


@pytest.fixture(scope="module")
def hyperbolic_horn_results(tmp_path_factory):
    """Solve V5 hyperbolic horn once and cache (frequencies, spl, ref)."""
    ref = load_reference("hyperbolic_horn_webster.json")
    geom = ref["geometry"]
    freq_cfg = ref["frequency_range"]
    tmp = tmp_path_factory.mktemp("v5")
    step_file = tmp / "hyperbolic_horn.step"
    create_hyperbolic_horn(
        throat_radius=geom["throat_radius_m"],
        mouth_radius=geom["mouth_radius_m"],
        length=geom["length_m"],
        output_file=step_file,
    )
    frequencies, spl = run_solver_and_get_spl(
        step_file=step_file,
        freq_range=(freq_cfg["min_hz"], freq_cfg["max_hz"]),
        num_intervals=freq_cfg["num_points"],
        horn_length=geom["length_m"],
        mesh_size=0.005,
        tmp_dir=tmp,
    )
    return frequencies, spl, ref
