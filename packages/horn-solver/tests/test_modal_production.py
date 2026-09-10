"""The optional modal path must conserve power and preserve coupled outputs."""
import numpy as np
import pandas as pd
import pytest
from horn_geometry.generator import create_horn
from horn_solver.solver import run_simulation_from_step


def test_modal_pressure_and_velocity_transfers_share_impedance_and_power(tmp_path, monkeypatch):
    from scipy.sparse.linalg import splu
    import horn_solver.modal_boundary as boundary
    solve = boundary.solve_sparse_system

    def compare_backends(system, forcing):
        result = solve(system, forcing)
        reference = splu(system).solve(forcing)
        np.testing.assert_allclose(result, reference, rtol=1e-9, atol=1e-10)
        return result

    monkeypatch.setattr(boundary, 'solve_sparse_system', compare_backends)
    step=tmp_path/'tube.step'
    create_horn('conical',.02,.02,.08,step,num_sections=4)
    frames={}
    for mode in ('dirichlet','velocity'):
        output=tmp_path/f'{mode}.csv'
        run_simulation_from_step(str(step),(400.,800.),3,{'length':.08},str(output),800.,
            mesh_size=.008,radiation_model='modal_baffled',bc_mode=mode)
        frame=pd.read_csv(output);frames[mode]=frame
        assert (frame.modal_mode_count==16).all()
        assert (frame.modal_interface_relative_error<1e-8).all()
        assert (frame.relative_residual<1e-8).all()
        np.testing.assert_allclose(frame.input_acoustic_power_w,frame.mouth_acoustic_power_w,rtol=1e-8)
        assert (frame.mouth_acoustic_power_w>0).all()
    a=frames['dirichlet'];b=frames['velocity']
    np.testing.assert_allclose(a.z_real+1j*a.z_imag,b.z_real+1j*b.z_imag,rtol=.002)


@pytest.mark.parametrize('shape',['annulus','offset','square','internal_cut'])
def test_modal_cad_rejects_unsupported_apertures_and_interiors(tmp_path,shape):
    import gmsh
    gmsh.initialize()
    try:
        if shape=='square':gmsh.model.occ.addBox(-.02,-.02,0,.04,.04,.08)
        else:
            volume=gmsh.model.occ.addCylinder(.01 if shape=='offset' else 0,0,0,0,0,.08,.02)
            if shape=='annulus':
                hole=gmsh.model.occ.addCylinder(0,0,0,0,0,.08,.005)
                gmsh.model.occ.cut([(3,volume)],[(3,hole)])
            elif shape=='internal_cut':
                hole=gmsh.model.occ.addBox(.015,-.005,.02,.02,.01,.02)
                gmsh.model.occ.cut([(3,volume)],[(3,hole)])
        gmsh.model.occ.synchronize();step=tmp_path/f'{shape}.step';gmsh.write(str(step))
    finally:gmsh.finalize()
    with pytest.raises(ValueError,match='Modal aperture'):
        run_simulation_from_step(str(step),(400.,800.),2,{'length':.08},str(tmp_path/'unused.csv'),800.,mesh_size=.008,radiation_model='modal_baffled')
