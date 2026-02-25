import pytest
from pathlib import Path

# Try to import FreeCAD, skip if not available
try:
    import FreeCAD
    import Part
except ImportError:
    FreeCAD = None

from horn_geometry.generator import create_conical_horn

reason = "FreeCAD module not found. Run this test inside the geometry container."
@pytest.mark.skipif(FreeCAD is None, reason=reason)
def test_generator_creates_step_file_smoke_test(tmp_path):
    """
    A simple smoke test to ensure the geometry creation function
    runs without crashing and creates a non-empty output file.
    It does not validate the contents of the geometry.
    """
    output_dir = tmp_path / "smoke_test_output"
    output_dir.mkdir()

    # Run the generator function directly
    output_file = create_conical_horn(
        throat_radius=0.025,
        mouth_radius=0.25,
        length=0.35,
        output_file=output_dir / "horn.step",
    )

    # Assert the file was created
    assert output_file.exists()
    assert output_file.stat().st_size > 0 