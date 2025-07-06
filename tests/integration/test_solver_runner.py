import pytest
from pathlib import Path
from unittest.mock import patch
import json
from horn.simulation import solver_runner

@pytest.mark.integration
@pytest.mark.solver
def test_solver_runner(tmp_path: Path):
    """
    An integration test for the solver runner script.

    This test should be run inside the horn-solver-app container.
    It verifies that the runner can correctly parse arguments and
    call the solver function.
    """
    # Arrange:
    # 1. Create a dummy mesh file for the solver to use.
    dummy_mesh_file = tmp_path / "dummy.msh"
    dummy_mesh_file.touch()

    # 2. Mock the actual run_simulation function.
    with patch('horn.simulation.solver_runner.solver.run_simulation') as mock_run_sim:
        mock_run_sim.return_value = str(tmp_path / "results.csv")

        # 3. Prepare the command-line arguments
        driver_params = {"fs": 40}
        driver_params_json = json.dumps(driver_params)
        args = [
            "solver_runner.py",
            "--mesh-file", str(dummy_mesh_file),
            "--output-dir", str(tmp_path),
            "--freq-min", "100",
            "--freq-max", "1000",
            "--driver-params-json", driver_params_json,
        ]

        # 4. Patch sys.argv to simulate command-line execution
        with patch('sys.argv', args):
            # --- Act ---
            solver_runner.main()

            # --- Assert ---
            # Verify that the runner called the solver function with the correct args
            mock_run_sim.assert_called_once()
            call_args, call_kwargs = mock_run_sim.call_args
            assert call_kwargs['mesh_file'] == str(dummy_mesh_file)
            assert call_kwargs['freq_range'] == (100.0, 1000.0)
            assert call_kwargs['driver_params'] == driver_params
            assert call_kwargs['output_dir'] == tmp_path 