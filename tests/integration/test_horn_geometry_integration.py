import pytest
import os
from horn.geometry import horn_generator
from horn.geometry.parameters import HornParameters, FlareProfile

# Mark this test as 'integration' and specifically for the 'freecad' environment
pytestmark = [pytest.mark.integration, pytest.mark.freecad]

def test_create_horn_integration(tmp_path):
    """
    An integration test for the create_horn function.
    
    This test should be run inside the Docker container where FreeCAD is
    installed. It verifies that a real, non-empty STEP file is created.
    
    Args:
        tmp_path: A pytest fixture providing a temporary directory path
                  within the container.
    """
    # Arrange
    params = HornParameters(
        flare_profile=FlareProfile.CONICAL,
        throat_radius=0.01,
        mouth_radius=0.1,
        length=0.2
    )
    output_dir = tmp_path

    # Act
    output_path = horn_generator.create_horn(params, output_dir=str(output_dir))

    # Assert
    # 1. Check that the file exists
    assert os.path.exists(output_path), "Output STEP file was not created."
    
    # 2. Check that the file is not empty
    file_size = os.path.getsize(output_path)
    assert file_size > 0, "Output STEP file is empty."
    
    # 3. Check that the filename is correct
    assert "conical" in output_path
    assert "0.2m" in output_path 