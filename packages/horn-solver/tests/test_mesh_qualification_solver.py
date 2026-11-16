"""Mesh qualification and directivity through the production solver (needs dolfinx)."""
import json
import sys

import numpy as np
import pandas as pd
import pytest
from horn_geometry.generator import create_horn
from horn_solver import solver


def test_response_records_the_mesh_actually_used(tmp_path):
    step = tmp_path / "horn.step"
    create_horn("conical", .0127, .04, .08, step, num_sections=8)
    output = tmp_path / "out.csv"
    solver.run_simulation_from_step(str(step), (2000., 4000.), 2, {"length": .08}, str(output), 4000.,
                                    mesh_size=.01, radiation_model="modal_baffled")
    frame = pd.read_csv(output)
    assert (frame.mesh_elements_per_wavelength >= 10 - 1e-9).all()
    assert (frame.mesh_throat_size_m <= .0127 / 4 + 1e-12).all()
    assert (abs(frame.mesh_throat_area_ratio - 1) <= .01).all()
    assert not frame.mesh_coarsened_for_node_limit.any()


def test_oversized_mesh_is_written_as_a_rejection_not_a_failure(tmp_path, monkeypatch):
    step = tmp_path / "horn.step"
    create_horn("conical", .0127, .04, .08, step, num_sections=8)
    output, rejection = tmp_path / "out.csv", tmp_path / "out.rejected.json"
    monkeypatch.setattr(sys, "argv", ["horn-solver", "--step-file", str(step), "--output-file", str(output),
                                      "--rejection-file", str(rejection), "--min-freq", "2000",
                                      "--max-freq", "4000", "--num-intervals", "2", "--length", ".08",
                                      "--radiation-model", "modal_baffled", "--max-mesh-nodes", "50"])
    solver.main()
    assert not output.exists()
    assert json.loads(rejection.read_text())["rejection_reasons"] == ["mesh_node_limit_exceeded"]
    # Without a rejection file the same request still fails loudly.
    monkeypatch.setattr(sys, "argv", [a for a in sys.argv if a not in ("--rejection-file", str(rejection))])
    with pytest.raises(Exception, match="nodes"):
        solver.main()


def test_cli_lists_only_supported_radiation_models(capsys, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["horn-solver", "--help"])
    with pytest.raises(SystemExit):
        solver.main()
    text = capsys.readouterr().out
    assert "bem" not in text and "directivity" not in text
    assert "modal_baffled" in text


@pytest.mark.validation
def test_os_beamwidth_varies_less_than_conical_with_the_same_mouth(tmp_path):
    """WP10 acceptance: an OS waveguide holds its beamwidth better than a cone."""
    from horn_analysis.directivity import directivity_table
    from horn_analysis.directivity_plot import compute_beamwidth
    spread = {}
    for profile in ("os", "conical"):
        step = tmp_path / f"{profile}.step"
        create_horn(profile, .0127, .1, .12, step)
        output = tmp_path / f"{profile}.csv"
        solver.run_simulation_from_step(str(step), (2000., 8000.), 9, {"length": .12}, str(output), 8000.,
                                        mesh_size=.01, radiation_model="modal_baffled")
        f, width = compute_beamwidth(directivity_table(pd.read_csv(output)))
        spread[profile] = float(np.std(np.log(width)))
    assert spread["os"] < spread["conical"]
