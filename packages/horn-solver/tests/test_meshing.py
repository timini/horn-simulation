"""Mesh qualification rules. Needs gmsh but not dolfinx, so it runs locally."""
import json

import numpy as np
import pytest

gmsh = pytest.importorskip("gmsh")
from horn_core.cad_units import import_step_metres
from horn_geometry.generator import create_horn
from horn_solver import meshing

C = 343.0


def test_element_size_depends_on_degree():
    assert meshing.wavelength_mesh_size(C, 8000., 1) == pytest.approx(C / 80000.)
    assert meshing.wavelength_mesh_size(C, 8000., 2) == pytest.approx(C / 48000.)
    # The user size still wins when it is finer than the wavelength rule.
    assert meshing.target_mesh_size(0.002, C, 8000., 1) == 0.002
    with pytest.raises(ValueError):
        meshing.wavelength_mesh_size(C, 8000., 3)


def test_node_limit_may_coarsen_p1_only_to_lambda_over_six():
    assert meshing.coarsest_mesh_size(0.01, C, 8000., 1) == pytest.approx(C / 48000.)
    assert meshing.coarsest_mesh_size(0.01, C, 8000., 2) == pytest.approx(C / 48000.)
    # A user size coarser than lambda/10 but finer than lambda/6 is never exceeded.
    assert meshing.coarsest_mesh_size(0.005, C, 8000., 1) == pytest.approx(0.005)


def test_documented_dispersion_errors():
    assert meshing.dispersion_phase_error(1, 6) == pytest.approx(0.0457, abs=5e-4)
    assert meshing.dispersion_phase_error(1, 10) == pytest.approx(0.0164, abs=5e-4)
    assert meshing.dispersion_phase_error(2, 6) == pytest.approx(0.00084, abs=5e-5)


def test_node_limit_reads_environment():
    assert meshing.max_mesh_nodes({}) == meshing.DEFAULT_MAX_MESH_NODES
    assert meshing.max_mesh_nodes({"HORN_MAX_MESH_NODES": "1234"}) == 1234


@pytest.fixture(scope="module")
def small_throat_horn(tmp_path_factory):
    step = tmp_path_factory.mktemp("mesh") / "horn.step"
    create_horn("conical", 0.0127, 0.05, 0.08, step, num_sections=20)
    return step


class Model:
    def __init__(self, step, length):
        self.step, self.length = step, length

    def __enter__(self):
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        import_step_metres(gmsh, self.step, self.length)
        surfaces = [s for _, s in gmsh.model.occ.getEntities(2)]
        self.inlet = [s for s in surfaces
                      if np.allclose([gmsh.model.getBoundingBox(2, s)[i] for i in (2, 5)], 0., atol=1e-6)]
        self.area = sum(gmsh.model.occ.getMass(2, s) for s in self.inlet)
        return self

    def __exit__(self, *exc):
        gmsh.finalize()


def longest_inlet_edge():
    tags, coords, _ = gmsh.model.mesh.getNodes()
    xyz = dict(zip(tags.astype(int), coords.reshape(-1, 3)))
    longest = 0.
    for surface in [s for _, s in gmsh.model.getEntities(2)
                    if np.allclose([gmsh.model.getBoundingBox(2, s)[i] for i in (2, 5)], 0., atol=1e-6)]:
        _, _, nodes = gmsh.model.mesh.getElements(2, surface)
        for a, b, c in np.asarray(nodes[0], dtype=int).reshape(-1, 3):
            for p, q in ((a, b), (b, c), (c, a)):
                longest = max(longest, float(np.linalg.norm(xyz[p] - xyz[q])))
    return longest


def test_throat_is_refined_to_a_fifth_radius_and_area_guarded(small_throat_horn, capsys):
    with Model(small_throat_horn, 0.08) as model:
        info = meshing.generate_qualified_mesh(gmsh, model.inlet, model.area, user_size=0.01,
                                               c=C, f_max=1000., degree=1, max_nodes=10**6)
        edge = longest_inlet_edge()
    # The first attempt is the documented rule, not the retry.
    assert info["throat_mesh_size_m"] == pytest.approx(0.0127 / 5)
    assert "remeshing with the throat cap halved" not in capsys.readouterr().out
    assert info["mesh_size_m"] == pytest.approx(0.01)
    # gmsh edges can exceed the target size somewhat; a uniform 10 mm mesh
    # would put edges near 10 mm across a 25 mm throat.
    assert edge < 1.5 * 0.0127 / 5
    # Margin inside the 1 % guard: r_t/4 gives about 0.990, a 0.03 % margin.
    assert info["throat_area_ratio"] >= 0.993
    assert abs(info["throat_area_ratio"] - 1) <= meshing.THROAT_AREA_TOLERANCE


def test_throat_area_retry_is_logged(small_throat_horn, monkeypatch, capsys):
    monkeypatch.setattr(meshing, "THROAT_AREA_TOLERANCE", 0.0)
    with Model(small_throat_horn, 0.08) as model:
        with pytest.raises(meshing.ThroatAreaGuardError):
            meshing.generate_qualified_mesh(gmsh, model.inlet, model.area, user_size=0.01,
                                            c=C, f_max=1000., degree=1, max_nodes=10**6)
    assert "remeshing with the throat cap halved" in capsys.readouterr().out


def test_coarse_throat_fails_the_one_percent_area_guard(small_throat_horn, monkeypatch):
    monkeypatch.setattr(meshing, "THROAT_SIZE_FRACTION", 4.0)
    with Model(small_throat_horn, 0.08) as model:
        with pytest.raises(meshing.ThroatAreaGuardError) as caught:
            meshing.generate_qualified_mesh(gmsh, model.inlet, model.area, user_size=0.03,
                                            c=C, f_max=1000., degree=1, max_nodes=10**6)
    assert caught.value.reason == "throat_area_guard_failed"


def test_node_limit_coarsens_p1_then_rejects(small_throat_horn):
    with Model(small_throat_horn, 0.08) as model:
        fine = meshing.generate_qualified_mesh(gmsh, model.inlet, model.area, user_size=0.01,
                                               c=C, f_max=8000., degree=1, max_nodes=10**6)
        assert not fine["mesh_coarsened_for_node_limit"]
        assert fine["elements_per_wavelength"] == pytest.approx(10.)
        limit = int(fine["mesh_nodes"] * 0.6)
        coarse = meshing.generate_qualified_mesh(gmsh, model.inlet, model.area, user_size=0.01,
                                                 c=C, f_max=8000., degree=1, max_nodes=limit)
        assert coarse["mesh_coarsened_for_node_limit"]
        assert coarse["mesh_nodes"] <= limit
        assert 6. - 1e-9 <= coarse["elements_per_wavelength"] < 10.
        with pytest.raises(meshing.MeshTooLargeError) as caught:
            meshing.generate_qualified_mesh(gmsh, model.inlet, model.area, user_size=0.01,
                                            c=C, f_max=8000., degree=1, max_nodes=200)
    assert caught.value.reason == "mesh_node_limit_exceeded"


def test_p2_is_never_coarsened_below_lambda_over_six(small_throat_horn):
    with Model(small_throat_horn, 0.08) as model:
        info = meshing.generate_qualified_mesh(gmsh, model.inlet, model.area, user_size=0.01,
                                               c=C, f_max=8000., degree=2, max_nodes=10**6)
        assert info["elements_per_wavelength"] == pytest.approx(6.)
        with pytest.raises(meshing.MeshTooLargeError):
            meshing.generate_qualified_mesh(gmsh, model.inlet, model.area, user_size=0.01,
                                            c=C, f_max=8000., degree=2, max_nodes=info["mesh_nodes"] - 1)


def test_rejection_record_names_the_reason(tmp_path):
    path = tmp_path / "rejected.json"
    meshing.write_rejection(path, meshing.MeshTooLargeError("big"), max_freq_hz=8000.)
    record = json.loads(path.read_text())
    assert record["rejection_reasons"] == ["mesh_node_limit_exceeded"]
    assert record["max_freq_hz"] == 8000.
    assert meshing.rejection_reason(RuntimeError("solver fault")) is None
    with pytest.raises(ValueError):
        meshing.write_rejection(path, RuntimeError("solver fault"))
