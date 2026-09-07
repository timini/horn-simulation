import pytest
from pathlib import Path

from horn_geometry.generator import create_conical_horn

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