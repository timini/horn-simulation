import pytest
import numpy as np
from pathlib import Path

# Try to import FreeCAD, skip if not available
try:
    import FreeCAD
    import Part
except ImportError:
    FreeCAD = None

from horn_geometry.generator import create_conical_horn

reason = "FreeCAD module not found. Run this test inside the geometry container."
@pytest.mark.skipif(FreeCAD is None, reason=reason)
class TestConicalHornValidation:
    """
    These tests validate the geometric properties of the generated
    conical horn STEP file. They require FreeCAD to be installed.
    """
    
    @pytest.fixture(scope="class")
    def horn_parameters(self):
        return {
            "throat_radius": 0.025,
            "mouth_radius": 0.1,
            "length": 0.2,
        }

    @pytest.fixture(scope="class")
    def generated_shape(self, tmp_path_factory, horn_parameters):
        """
        Generates the conical horn once and returns the FreeCAD shape object.
        """
        output_dir = tmp_path_factory.mktemp("horn_test")
        step_file = create_conical_horn(
            output_dir=output_dir,
            **horn_parameters
        )
        # Load the shape from the generated STEP file
        shape = Part.Shape()
        shape.read(str(step_file))
        return shape

    def test_horn_volume(self, generated_shape, horn_parameters):
        """
        Validates the volume of the generated horn against the analytical
        formula for a conical frustum.
        V = (1/3) * pi * h * (R^2 + R*r + r^2)
        """
        r = horn_parameters["throat_radius"]
        R = horn_parameters["mouth_radius"]
        h = horn_parameters["length"]
        
        expected_volume = (1/3) * np.pi * h * (R**2 + R*r + r**2)
        
        # Check if the generated volume is within a small tolerance
        assert np.isclose(generated_shape.Volume, expected_volume, rtol=1e-3)

    def test_horn_faces(self, generated_shape, horn_parameters):
        """
        Validates the circular faces at the throat and mouth.
        Checks for their existence, area, and center position.
        """
        r_throat = horn_parameters["throat_radius"]
        r_mouth = horn_parameters["mouth_radius"]
        length = horn_parameters["length"]

        throat_area_expected = np.pi * r_throat**2
        mouth_area_expected = np.pi * r_mouth**2

        # Find the circular faces
        circular_faces = [f for f in generated_shape.Faces if isinstance(f.Surface, Part.Cylinder) or f.Curve.isclosed()]
        
        # We expect to find two planar, circular faces
        planar_circles = []
        for face in generated_shape.Faces:
            if face.Surface.TypeId == "Part::GeomPlane" and len(face.Wires) == 1 and face.Wires[0].Curve.TypeId == "Part::GeomCircle":
                planar_circles.append(face)

        assert len(planar_circles) == 2, "Expected to find exactly two planar, circular faces"

        # Identify throat and mouth by area
        face1, face2 = planar_circles
        if face1.Area > face2.Area:
            mouth_face, throat_face = face1, face2
        else:
            mouth_face, throat_face = face2, face1

        # Validate areas
        assert np.isclose(throat_face.Area, throat_area_expected, rtol=1e-3)
        assert np.isclose(mouth_face.Area, mouth_area_expected, rtol=1e-3)

        # Validate positions (assuming horn is aligned along the Y-axis)
        throat_center = throat_face.CenterOfMass
        mouth_center = mouth_face.CenterOfMass
        
        assert np.isclose(throat_center.y, 0)
        assert np.isclose(mouth_center.y, length)
        # X and Z should be at the origin
        assert np.isclose(throat_center.x, 0) and np.isclose(throat_center.z, 0)
        assert np.isclose(mouth_center.x, 0) and np.isclose(mouth_center.z, 0) 