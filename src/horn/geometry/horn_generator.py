import os
from horn.geometry.parameters import HornParameters

# The 'FreeCAD' and 'Part' libraries are expected to be available in the
# environment where this script is run (e.g., inside our Docker container).
# We add a try-except block to allow the code to be imported in other
# environments (like for unit testing with mocks) without raising an error.
try:
    import FreeCAD
    import Part
except ImportError:
    # If the import fails, we define dummy variables. This allows 'patch'
    # to find these names on the module and replace them with mocks for
    # the unit tests, while still raising an error if the code is run
    # for real without FreeCAD installed.
    FreeCAD = None
    Part = None
    print("Warning: FreeCAD libraries not found. This is expected during unit testing outside of the Docker container.")


def create_horn(params: HornParameters, output_dir: str = "/tmp") -> str:
    """
    Generates a 3D model of a horn and exports it as a STEP file.

    This implementation uses the FreeCAD scripting API to create the geometry.

    Args:
        params: A dataclass object containing the horn's geometric parameters.
        output_dir: The directory where the output file will be saved.

    Returns:
        The path to the generated STEP file.
    """
    if not FreeCAD or not Part:
        raise ImportError("FreeCAD libraries are required for geometry generation.")

    print(f"Generating horn with parameters: {params}")

    # Create a new document
    doc = FreeCAD.newDocument("horn")

    # For now, we'll create a simple cylinder as a proof-of-concept
    # to ensure the FreeCAD API is working correctly.
    # The real lofting logic will be added next.
    radius = params.throat_radius
    length = params.length
    cylinder = Part.makeCylinder(radius, length)

    # Add the part to the document
    part_obj = doc.addObject("Part::Feature", "Horn")
    part_obj.Shape = cylinder

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Define output path
    file_name = f"horn_{params.flare_profile.value}_{params.length}m.stp"
    output_path = os.path.join(output_dir, file_name)

    # Export the document to a STEP file
    Part.export([part_obj], output_path)
    
    print(f"Successfully generated STEP file at: {output_path}")
    
    return output_path
