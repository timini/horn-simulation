import pytest
import sys
from pathlib import Path

# Try to import FreeCAD, skip if not available
try:
    import FreeCAD
    import Part
except ImportError:
    FreeCAD = None

from horn_geometry.generator import main

reason = "FreeCAD module not found. Run this test inside the geometry container."
@pytest.mark.skipif(FreeCAD is None, reason=reason)
def test_generator_creates_step_file(tmp_path, monkeypatch):
    output_file = tmp_path / "test_horn.stp"
    
    # Mock sys.argv to pass command-line arguments to the main function
    test_args = [
        "script_name",  # sys.argv[0] is the script name
        "--profile", "conical",
        "--throat", "0.025",
        "--mouth", "0.25",
        "--length", "0.35",
        "--output-file", str(output_file),
    ]
    monkeypatch.setattr(sys, "argv", test_args)

    # Run the generator
    main()

    # Assert the file was created
    assert output_file.exists()
    assert output_file.stat().st_size > 0 