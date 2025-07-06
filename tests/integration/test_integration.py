import pytest
import os
from horn.geometry import horn_generator
from horn.geometry.parameters import HornParameters, FlareProfile
from horn.simulation import meshing

# Mark this test as 'integration' and for the 'gmsh' environment
pytestmark = [pytest.mark.integration, pytest.mark.gmsh]

def test_create_mesh_integration(tmp_path):
    """
    An integration test for the create_mesh function.
    
    This test should be run inside the Docker container where both FreeCAD
    and Gmsh are installed. It verifies that a real STEP file can be
    successfully turned into a real, non-empty MSH file.
    
    Args:
        tmp_path: A pytest fixture providing a temporary directory path
                  within the container.
    """
    # Arrange: First, create a real STEP file to use as input
    horn_params = HornParameters(
        flare_profile=FlareProfile.CONICAL,
        throat_radius=0.01,
        mouth_radius=0.1,
        length=0.2
    )
    step_file_path = horn_generator.create_horn(horn_params, output_dir=str(tmp_path))
    assert os.path.exists(step_file_path)

    # Act: Call the real meshing function
    mesh_file_path = meshing.create_mesh(step_file_path, max_freq=1000, output_dir=str(tmp_path))

    # Assert
    # 1. Check that the mesh file exists
    assert os.path.exists(mesh_file_path), "Output MSH file was not created."
    
    # 2. Check that the file is not empty
    file_size = os.path.getsize(mesh_file_path)
    assert file_size > 0, "Output MSH file is empty."
    
    # 3. Check that the filename is correct
    assert mesh_file_path.endswith(".msh") 