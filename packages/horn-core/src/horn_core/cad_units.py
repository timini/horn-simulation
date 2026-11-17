"""STEP length units shared by every gmsh/OCC reader and writer in the pipeline.

The pipeline models horns in metres. gmsh's OCC kernel defaults to
millimetres when it reads or writes STEP, so a model built in metres and
written with the default is labelled millimetres while holding metre
values, and opens 1000 times too small in CAD.

Setting ``Geometry.OCCTargetUnit`` to metres makes OCC convert at the file
boundary. On export it writes millimetre coordinates under a millimetre
unit declaration (the file is true size). On import it converts whatever
unit the file declares to metres, so CAD exports in millimetres, metres
or inches all arrive at the right size.

This module has no gmsh dependency; callers pass their gmsh module.
"""

STEP_LENGTH_UNIT = "M"


def set_step_units(gmsh_module) -> None:
    """Make OCC STEP import and export use metres. Call after gmsh.initialize()."""
    gmsh_module.option.setString("Geometry.OCCTargetUnit", STEP_LENGTH_UNIT)


def legacy_unit_message(z_min: float, z_max: float, horn_length: float, rel_tol: float = 0.02):
    """Explain an imported model whose axial extent is 1000x shorter than expected.

    Before STEP units were honoured, the pipeline wrote metre coordinates
    under a millimetre declaration. Read with units honoured, such a file
    comes out 1000 times too small. Returns a message for that case, or
    None when the extent is not a factor of 1000 short. The tolerance allows
    for the small padding gmsh adds to bounding boxes.
    """
    extent = z_max - z_min
    if horn_length > 0 and abs(extent * 1000.0 - horn_length) <= rel_tol * horn_length:
        return (f"STEP model is {extent:.6g} m long, 1/1000 of the requested length "
                f"{horn_length:.6g} m. The file declares millimetres but holds metre values "
                "(written by an older pipeline version). Re-export it with true units, or "
                "change its SI_UNIT(.MILLI.,.METRE.) declaration to SI_UNIT($,.METRE.).")
    return None


def import_step_metres(gmsh_module, step_file: str, horn_length: float) -> None:
    """Import a STEP file into the current gmsh model in metres, honouring its units.

    Raises ValueError (after finalising gmsh) for a legacy file whose
    declared unit disagrees with its values by a factor of 1000.
    """
    set_step_units(gmsh_module)
    gmsh_module.model.occ.importShapes(str(step_file))
    gmsh_module.model.occ.synchronize()
    box = gmsh_module.model.getBoundingBox(-1, -1)
    message = legacy_unit_message(box[2], box[5], horn_length)
    if message:
        gmsh_module.finalize()
        raise ValueError(message)
