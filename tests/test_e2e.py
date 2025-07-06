import pytest
import tempfile
from unittest.mock import patch, call
from horn import main
from horn.geometry.parameters import HornParameters, FlareProfile

@pytest.mark.e2e
@pytest.mark.db
@patch('horn.main.subprocess.run')
@patch('horn.main.driver_db.get_driver_parameters')
def test_full_pipeline_orchestration_e2e(mock_get_driver, mock_subprocess_run):
    """
    An end-to-end test for the pipeline orchestrator.

    This test mocks the database and the Docker calls to verify that the main
    orchestrator correctly sequences the stages and calls the containers
    with the correct arguments and volume mounts.
    """
    # --- Arrange ---
    # 1. Mock the database call
    mock_get_driver.return_value = {
        'manufacturer': 'TestCorp',
        'model_name': 'Test-1'
    }

    # 2. Define pipeline inputs
    driver_id = "test_driver_001"
    horn_params = HornParameters(
        flare_profile=FlareProfile.CONICAL,
        throat_radius=0.01,
        mouth_radius=0.1,
        length=0.2
    )
    freq_range = (100, 1000)

    # --- Act ---
    final_report = main.run_pipeline(driver_id, horn_params, freq_range)

    # --- Assert ---
    # 1. Verify the Docker call for the FreeCAD container
    freecad_call = mock_subprocess_run.call_args_list[0]
    freecad_args = freecad_call.args[0]
    assert "horn-freecad-app" in freecad_args
    assert "/data" in freecad_args[-1] # output-dir
    assert horn_params.flare_profile.value in freecad_args
    assert str(horn_params.length) in freecad_args
    
    # 2. Verify the Docker call for the Solver container (meshing stage)
    solver_call = mock_subprocess_run.call_args_list[1]
    solver_args = solver_call.args[0]
    assert "horn-solver-app" in solver_args
    assert "horn.simulation.meshing_runner" in solver_args
    assert solver_args[-6] == "--step-file"
    assert "horn_conical_0.2m.stp" in solver_args[-5]
    assert solver_args[-2] == "--max-freq"
    
    # 3. Verify the final (placeholder) report
    assert isinstance(final_report, dict)
    assert "metrics" in final_report
    assert "plots" in final_report
    assert final_report["metrics"]["score"] > 0
    assert isinstance(final_report["plots"], list)
    assert len(final_report["plots"]) == 1
    
    # Check that the plot paths are plausible
    plot_path = final_report["plots"][0]
    assert plot_path.startswith(tempfile.gettempdir())
    assert "spl.png" in plot_path 