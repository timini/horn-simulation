"""Checks that the scientific audit cannot hide invalid or misleading evidence."""
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location("reference_audit", SCRIPTS/"validate_reference_audit.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def test_phase_wrap_does_not_turn_two_degrees_into_358():
    result = audit.compare(np.exp(1j*np.deg2rad([179.])), np.exp(1j*np.deg2rad([-179.])))
    assert result["p95_phase_error_deg"] == pytest.approx(2)


def test_magnitude_pass_does_not_hide_wrong_phase():
    result = audit.compare([1j, 1j], [-1j, -1j])
    assert result["magnitude_agreement"]
    assert result["p95_phase_error_deg"] == pytest.approx(180)
    assert result["p95_scaled_complex_error"] == pytest.approx(1)


@pytest.mark.parametrize("frequency,real", [([120, 3900], [1, 1]),
    ([110, 3800], [1, 1]), ([110, 110, 3900], [1, 1, 1]),
    ([3900, 110], [1, 1]), ([110, 3900], [1, np.nan])])
def test_invalid_reference_cannot_be_extrapolated_or_silently_reordered(frequency, real):
    frame = pd.DataFrame({"frequency": frequency, "real": real, "imag": np.zeros(len(real))})
    with pytest.raises(ValueError, match="entire fixed band"):
        audit.interpolate(frame, np.array([110, 3900]), "real", "imag")


def healthy_frame():
    return pd.DataFrame({"frequency": [110., 3900.], "z_real": [1., 1.], "z_imag": [0., 0.],
        "input_acoustic_power_w": [1., 1.], "mouth_acoustic_power_w": [.5, .5],
        "viscous_wall_power_w": [.3, .3], "thermal_wall_power_w": [.2, .2],
        "relative_residual": [1e-12, 1e-12], "converged_reason": [4, 4], "mesh_cells": [100, 100]})


@pytest.mark.parametrize("column,value", [("z_real", -1), ("relative_residual", 1e-5),
    ("converged_reason", -3), ("input_acoustic_power_w", 2),
    ("viscous_wall_power_w", -.3), ("z_imag", np.inf)])
def test_solver_health_rejects_unreliable_results(column, value):
    frame = healthy_frame()
    assert audit.fem_health(frame)["passed"]
    frame.loc[0, column] = value
    assert not audit.fem_health(frame)["passed"]


def test_all_simulation_configurations_have_explicit_normalization_geometry():
    for case in audit.NUMERICAL_CASE.values():
        assert audit.ALL_CASES[case][0] > 0
    assert audit.ALL_CASES["Cone_U"][3] == "finite_flange"
    assert audit.ALL_CASES["Cone_U"][2] == 0
