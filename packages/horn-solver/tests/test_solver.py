import pytest
from pathlib import Path
from horn_solver.solver import run_simulation_from_step, create_mesh_from_step
import numpy as np

# Define physical group tags for boundaries
INLET_TAG, OUTLET_TAG, WALL_TAG = 2, 3, 4

def test_mesh_boundary_tagging():
    print("\n--- Running test: test_mesh_boundary_tagging ---\n")
    step_file = Path(__file__).parent / "test_box.stp"

    # Call the new, testable meshing function
    domain, facet_tags = create_mesh_from_step(
        step_file=str(step_file),
        mesh_size=1.0, # Coarse mesh for speed
    )

    # Check that the facet tags have been created and contain the correct markers
    assert facet_tags is not None
    
    # Get all unique tags present in the facet_tags
    all_tags = np.unique(facet_tags.values)
    print(f"Found unique facet tags: {all_tags}")

    # Assert that our expected boundary tags are present in the mesh
    assert INLET_TAG in all_tags, "Inlet boundary tag not found in mesh."
    assert OUTLET_TAG in all_tags, "Outlet boundary tag not found in mesh."
    assert WALL_TAG in all_tags, "Wall boundary tag not found in mesh."

    print("--- Test finished ---\n")

def test_e2e_meshing_and_solving(tmp_path):
    print("\n--- Running test: test_e2e_meshing_and_solving ---\n")
    # This test uses a pre-generated STEP file to avoid a dependency on FreeCAD
    step_file = Path(__file__).parent / "test_box.stp"
    output_file = tmp_path / "results.csv"
    print(f"STEP file: {step_file}")
    print(f"Output file: {output_file}")

    # Define dummy inputs
    driver_params = {"Bl": 5.0, "Re": 6.0}
    freq_range = (100.0, 1000.0)
    print("Inputs defined.")

    # Run the simulation from the STEP file.
    print("Calling run_simulation_from_step...")
    result_path = run_simulation_from_step(
        step_file=str(step_file),
        driver_params=driver_params,
        freq_range=freq_range,
        output_file=str(output_file),
        max_freq_mesh=freq_range[1],
        mesh_size=1.0,  # Use a very coarse mesh for speed
    )
    print("run_simulation_from_step finished.")

    # Assert the output file was created
    print("Asserting results...")
    assert result_path.exists()
    assert result_path.name == "results.csv"
    print("--- Test finished ---\n")