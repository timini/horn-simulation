"""Legacy STEP unit detection (no gmsh needed)."""
from horn_core.cad_units import STEP_LENGTH_UNIT, legacy_unit_message


def test_pipeline_uses_metres():
    assert STEP_LENGTH_UNIT == "M"


def test_thousandfold_short_model_is_explained():
    message = legacy_unit_message(0.0, 0.000117906, 0.117906)
    assert message and "millimetres but holds metre values" in message


def test_true_size_model_passes():
    assert legacy_unit_message(0.0, 0.117906, 0.117906) is None
    assert legacy_unit_message(0.0, 0.2, 0.117906) is None
