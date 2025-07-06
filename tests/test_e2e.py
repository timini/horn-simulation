import pytest
from unittest.mock import patch
from horn import main
from horn.geometry.parameters import HornParameters, FlareProfile

@pytest.mark.e2e
@pytest.mark.db
@patch('horn.main.horn_generator.create_horn')
@patch('horn.main.meshing.create_mesh')
def test_full_pipeline_e2e_mocked(mock_create_mesh, mock_create_horn):
    """
    An end-to-end test for the full pipeline, mocking the containerized
    geometry and meshing stages.

    This test verifies that the main orchestrator can correctly sequence
    all stages of the pipeline.
    """
    # Arrange: Configure the mocks to behave like the real functions
    mock_create_horn.return_value = "/tmp/mock_horn_conical_0.6m.stp"
    mock_create_mesh.return_value = "/tmp/mock_horn_conical_0.6m_1000hz.msh"
    
    # Define the inputs for the pipeline
    driver_id = "test_driver_001"
    horn_params = HornParameters(
        flare_profile=FlareProfile.CONICAL,
        throat_radius=0.03,
        mouth_radius=0.3,
        length=0.6
    )
    freq_range = (20, 1000) # Hz

    # Act: Run the main pipeline function
    final_report = main.run_pipeline(driver_id, horn_params, freq_range)

    # Assert
    # 1. Check that our mocks were called correctly
    mock_create_horn.assert_called_once()
    mock_create_mesh.assert_called_once()

    # 2. Verify the final output
    assert isinstance(final_report, dict)
    assert "metrics" in final_report
    assert "plots" in final_report
    assert "score" in final_report["metrics"]
    
    # 3. Check that the file paths are plausible, based on the mock output
    plot_path = final_report["plots"][0]
    assert "mock_horn_conical" in plot_path
    assert plot_path.endswith("_spl.png")

    # Verify the final output
    assert isinstance(final_report, dict)
    assert "metrics" in final_report
    assert "plots" in final_report
    assert "score" in final_report["metrics"]
    assert final_report["metrics"]["score"] > 0
    assert isinstance(final_report["plots"], list)
    assert len(final_report["plots"]) == 2
    
    # Check that the plot paths are plausible
    plot_path = final_report["plots"][0]
    assert plot_path.startswith("/tmp") # Check it's in a temporary-like directory
    assert plot_path.endswith("_spl.png")
    assert "conical" in plot_path 