import pytest
from unittest.mock import patch, MagicMock
from horn.geometry import horn_generator
from horn.geometry.parameters import HornParameters, FlareProfile

# Mock the FreeCAD modules at the module level where they are imported
# This ensures that any call to 'horn_generator.FreeCAD' or 'horn_generator.Part'
# will be intercepted by our mock objects.
@patch('horn.geometry.horn_generator.FreeCAD', MagicMock())
@patch('horn.geometry.horn_generator.Part', MagicMock())
def test_create_horn_with_mocked_freecad(tmp_path):
    """
    Test the create_horn function's logic by mocking the FreeCAD API.

    This test verifies that the function calls the correct FreeCAD methods
    with the correct parameters, without actually running FreeCAD.
    
    Args:
        tmp_path: A pytest fixture providing a temporary directory path.
    """
    # Arrange
    params = HornParameters(
        flare_profile=FlareProfile.EXPONENTIAL,
        throat_radius=0.025,
        mouth_radius=0.25,
        length=0.5
    )
    output_dir = tmp_path / "test_output"
    
    # Act
    # We call the real function. The patches will intercept the FreeCAD calls.
    output_path = horn_generator.create_horn(params, output_dir=str(output_dir))

    # Assert
    # Check that the main FreeCAD functions were called
    horn_generator.FreeCAD.newDocument.assert_called_once_with("horn")
    doc_mock = horn_generator.FreeCAD.newDocument.return_value
    doc_mock.addObject.assert_called_with("Part::Feature", "Horn")
    
    # Check that the Part geometry creation function was called with correct args
    horn_generator.Part.makeCylinder.assert_called_once_with(
        params.throat_radius, params.length
    )
    
    # Check that the export function was called
    part_obj_mock = doc_mock.addObject.return_value
    horn_generator.Part.export.assert_called_once_with(
        [part_obj_mock], output_path
    )
    
    # Check that the output path is correctly formed
    expected_file = output_dir / "horn_exponential_0.5m.stp"
    assert output_path == str(expected_file)
