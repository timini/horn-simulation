import pytest
from pathlib import Path
from unittest.mock import patch
from horn.simulation import meshing_runner

@pytest.mark.integration
@pytest.mark.solver
def test_meshing_runner(tmp_path: Path):
    """
    An integration test for the meshing runner script.

    This test should be run inside the horn-solver-app container.
    It verifies that the runner can correctly parse arguments and
    call the meshing function to produce a mesh file.
    """
    # Arrange:
    # 1. Create a dummy input file for the mesher to use.
    #    In a real run, this would be a .stp file, but for this test,
    #    we just need a file to exist. We'll mock the meshing call itself.
    dummy_step_file = tmp_path / "dummy.stp"
    dummy_step_file.touch()

    # 2. Mock the actual create_mesh function, since we are only testing
    #    the runner script's logic, not the mesher itself.
    with patch('horn.simulation.meshing_runner.meshing.create_mesh') as mock_create_mesh:
        mock_create_mesh.return_value = str(tmp_path / "dummy.msh")

        # 3. Prepare the command-line arguments
        args = [
            "meshing_runner.py",
            "--step-file", str(dummy_step_file),
            "--output-dir", str(tmp_path),
            "--max-freq", "1000"
        ]

        # 4. Patch sys.argv to simulate command-line execution
        with patch('sys.argv', args):
            # --- Act ---
            meshing_runner.main()

            # --- Assert ---
            # Verify that the runner called the meshing function with the correct args
            mock_create_mesh.assert_called_once_with(
                step_file=str(dummy_step_file),
                output_dir=str(tmp_path),
                max_freq=1000.0
            ) 