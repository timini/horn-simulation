import os

try:
    import gmsh
except ImportError:
    gmsh = None
    print("Warning: gmsh library not found. This is expected during unit testing outside of the Docker container.")

# Speed of sound in air in meters/second
C0 = 343.0

def create_mesh(step_file: str, max_freq: float, output_dir: str = "/tmp") -> str:
    """
    Generates a 3D mesh from a STEP file using Gmsh.

    Args:
        step_file: The path to the input STEP file.
        max_freq: The maximum frequency (in Hz) to be simulated, used to
                  determine the required mesh density.
        output_dir: The directory where the output mesh file will be saved.

    Returns:
        The path to the generated MSH file.
    """
    if not gmsh:
        raise ImportError("Gmsh library is required for meshing.")

    gmsh.initialize()

    try:
        print(f"Meshing STEP file: {step_file}")
        gmsh.model.add(f"horn_model")
        gmsh.model.occ.importShapes(step_file)

        # Synchronize the CAD model with the Gmsh model
        gmsh.model.occ.synchronize()

        # For this basic implementation, we'll set a coarse global mesh size.
        # A real implementation would use max_freq to calculate a wavelength-
        # dependent element size. e.g. c0 / max_freq / points_per_wavelength
        gmsh.option.setNumber("Mesh.MeshSizeMax", 0.1)
        gmsh.option.setNumber("Mesh.MeshSizeMin", 0.01)

        # Generate the 3D volume mesh
        print("Generating 3D mesh...")
        gmsh.model.mesh.generate(3)

        # Define the output path
        base_name = os.path.splitext(os.path.basename(step_file))[0]
        output_path = os.path.join(output_dir, f"{base_name}.msh")

        # Write the mesh file
        # We use version 2.2 of the MSH format for compatibility with older solvers.
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.write(output_path)
        print(f"Successfully generated mesh file: {output_path}")

    finally:
        # Finalize Gmsh to clean up
        gmsh.finalize()

    return output_path
