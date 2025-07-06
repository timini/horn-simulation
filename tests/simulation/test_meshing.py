from horn.simulation import meshing

def test_create_mesh_is_callable():
    """
    Test that the create_mesh function exists and is callable.
    """
    assert callable(meshing.create_mesh)

def test_create_mesh_placeholder():
    """
    Test the placeholder implementation of create_mesh.
    It should accept a STEP file path and max frequency, and return a mock mesh file path.
    """
    step_file = "/tmp/dummy_horn.stp"
    max_freq = 1000  # Hz
    
    result = meshing.create_mesh(step_file, max_freq)
    
    assert isinstance(result, str)
    assert result.endswith(".msh")
    assert "1000hz" in result
