"""STEP files must carry their true size, and imports must honour declared units."""
import re

import gmsh
import numpy as np
import pytest

from horn_core.cad_units import import_step_metres
from horn_geometry.generator import create_horn

R_T, R_M, L = 0.032498, 0.08870898640584517, 0.11790625  # worked 800-1600 Hz example
_UNIT_SCALE = {"$": 1.0, ".MILLI.": 1e-3, ".CENTI.": 1e-2, ".KILO.": 1e3}


def _declared_metres_per_unit(text):
    prefixes = set(re.findall(r"LENGTH_UNIT\(\)\s*NAMED_UNIT\(\*\)\s*SI_UNIT\(([^,]+),\.METRE\.\)", text))
    assert len(prefixes) == 1, prefixes
    return _UNIT_SCALE[prefixes.pop()]


def _max_axial_coordinate(text):
    points = re.findall(r"CARTESIAN_POINT\('[^']*',\(([^)]*)\)\)", text)
    return max(float(c[2]) for c in (p.split(",") for p in points) if len(c) == 3)


def _ports(step, unit):
    """Throat/mouth areas and length read by a units-honouring OCC import."""
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setString("Geometry.OCCTargetUnit", unit)
        gmsh.model.occ.importShapes(str(step))
        gmsh.model.occ.synchronize()
        return _summarise()
    finally:
        gmsh.finalize()


def _summarise():
    discs = sorted((gmsh.model.occ.getCenterOfMass(d, t)[2], gmsh.model.occ.getMass(d, t))
                   for d, t in gmsh.model.getEntities(2) if gmsh.model.getType(d, t) == "Plane")
    return {"length": discs[-1][0] - discs[0][0], "throat_area": discs[0][1], "mouth_area": discs[-1][1]}


@pytest.fixture(scope="module")
def generated(tmp_path_factory):
    step = tmp_path_factory.mktemp("units") / "horn.step"
    return create_horn("hyperbolic", R_T, R_M, L, step, num_sections=20)


def test_declared_unit_matches_coordinates(generated):
    text = generated.read_text()
    metres = _max_axial_coordinate(text) * _declared_metres_per_unit(text)
    assert metres == pytest.approx(L, rel=1e-6)


def test_reimport_in_millimetres_gives_worked_example_dimensions(generated):
    ports = _ports(generated, "MM")
    assert ports["length"] == pytest.approx(117.906, abs=1e-3)
    assert 2 * np.sqrt(ports["throat_area"] / np.pi) == pytest.approx(64.996, abs=1e-3)
    assert 2 * np.sqrt(ports["mouth_area"] / np.pi) == pytest.approx(177.418, abs=1e-3)


def test_round_trip_in_metres(generated):
    ports = _ports(generated, "M")
    assert ports["length"] == pytest.approx(L, rel=1e-6)
    assert ports["throat_area"] == pytest.approx(np.pi * R_T**2, rel=1e-6)
    assert ports["mouth_area"] == pytest.approx(np.pi * R_M**2, rel=1e-6)


def _write_cad_export(path, unit, scale):
    """A CAD-style export: a cone built in `unit` and written in that unit."""
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setString("Geometry.OCCTargetUnit", unit)
        gmsh.model.occ.addCone(0, 0, 0, 0, 0, 0.08 * scale, 0.02 * scale, 0.05 * scale)
        gmsh.model.occ.synchronize()
        gmsh.write(str(path))
    finally:
        gmsh.finalize()


@pytest.mark.parametrize("unit,scale", [("MM", 1000.0), ("M", 1.0)])
def test_solver_import_path_reads_cad_exports_at_true_size(tmp_path, unit, scale):
    step = tmp_path / f"cad_{unit}.step"
    _write_cad_export(step, unit, scale)
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        import_step_metres(gmsh, step, 0.08)
        ports = _summarise()
    finally:
        gmsh.finalize()
    assert ports["length"] == pytest.approx(0.08, rel=1e-6)
    assert ports["throat_area"] == pytest.approx(np.pi * 0.02**2, rel=1e-6)
    assert ports["mouth_area"] == pytest.approx(np.pi * 0.05**2, rel=1e-6)


def test_legacy_metre_values_under_mm_declaration_are_rejected(tmp_path):
    step = tmp_path / "legacy.step"
    # The old pipeline: metre-valued model written with OCC's default millimetre label.
    _write_cad_export(step, "MM", 1.0)
    gmsh.initialize()
    with pytest.raises(ValueError, match="declares millimetres but holds metre values"):
        import_step_metres(gmsh, step, 0.08)
    assert not gmsh.isInitialized()
