"""Production thin-layer FEM against circular Kirchhoff/TMM and energy balance."""
import numpy as np
import pandas as pd
from horn_core.duct import AirProperties
from horn_core.webster import compute_horn_transfer_tmm
from horn_solver.solver import run_simulation_from_step
from .conftest import _generate_cylinder_step


def test_lossy_pipe_convergence_energy_and_complex_impedance(tmp_path):
    air=AirProperties(c=346.28592,rho=1.184490,gamma=1.402,viscosity=1.83183e-5,
                      conductivity=.02613336,heat_capacity=1004.16)
    step=tmp_path/'cylinder.step';_generate_cylinder_step(step,.007,.18)
    responses=[]
    for h in (.006,.004,.003):
        output=tmp_path/f'pipe-{h}.csv'
        run_simulation_from_step(str(step),(110,3900),31,{'length':.18},str(output),3900,
                                 mesh_size=h,element_degree=2,loss_model='boundary_layer',
                                 minimum_wall_scale=.007,radiation_model='finite_flange',
                                 flange_width=.002,air=air)
        d=pd.read_csv(output);responses.append(d)
        assert (d.z_real > 0).all()
        loss=d.viscous_wall_power_w+d.thermal_wall_power_w
        assert (loss>0).all()
        np.testing.assert_allclose(d.input_acoustic_power_w,d.mouth_acoustic_power_w+loss,rtol=1e-7,atol=1e-15)
    last=responses[-1]
    t=compute_horn_transfer_tmm(last.frequency.to_numpy(),lambda _: .007,.18,.007,.007,
                              loss_model='boundary_layer',radiation_model='finite_flange',flange_width=.002,air=air)
    actual=last.z_real+1j*last.z_imag;expected=t['z_real']+1j*t['z_imag']
    np.testing.assert_allclose(actual,expected,rtol=.05,atol=1.)
    assert max(abs(responses[-1].spl-responses[-2].spl)) < .5
    assert max(abs(responses[-1].spl-responses[-2].spl)) < max(abs(responses[-1].spl-responses[0].spl))
