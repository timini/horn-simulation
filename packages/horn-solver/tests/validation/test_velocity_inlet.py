"""Independent travelling-wave checks for the prescribed-velocity inlet."""
import numpy as np
import pandas as pd
import pytest
from horn_solver.solver import run_simulation, run_simulation_from_step
from .conftest import _generate_cylinder_step


pytestmark = pytest.mark.validation

def test_complex_velocity_travelling_wave_and_power(tmp_path):
    """Check phase, normalization and RMS energy against an exact tube wave."""
    radius, length = .01, .1
    step = _generate_cylinder_step(tmp_path / 'tube.step', radius, length)
    velocity = .003 + .004j
    output = run_simulation_from_step(
        str(step), (200., 800.), 5, {'length': length}, str(tmp_path/'response.csv'),
        800., mesh_size=.004, element_degree=2, bc_mode='velocity',
        inlet_velocity_rms=velocity, radiation_model='plane_wave')
    frame = pd.read_csv(output)
    inlet_p = frame.inlet_p_real.to_numpy() + 1j*frame.inlet_p_imag.to_numpy()
    inlet_u = frame.inlet_u_real.to_numpy() + 1j*frame.inlet_u_imag.to_numpy()
    mouth_p = frame.mouth_p_real.to_numpy() + 1j*frame.mouth_p_imag.to_numpy()
    expected = frame.air_rho_kg_m3.to_numpy()*frame.air_c_m_s.to_numpy()*velocity
    phase = np.exp(-2j*np.pi*frame.frequency.to_numpy()*length/frame.air_c_m_s.to_numpy())
    np.testing.assert_allclose(inlet_p, expected, rtol=2e-3)
    np.testing.assert_allclose(mouth_p, expected*phase, rtol=2e-3)
    np.testing.assert_allclose(inlet_u, velocity*frame.mesh_inlet_area_m2, rtol=1e-12)
    expected_power = np.real(expected*np.conj(velocity))*frame.mesh_inlet_area_m2
    np.testing.assert_allclose(frame.input_acoustic_power_w, expected_power, rtol=2e-3)
    np.testing.assert_allclose(frame.input_acoustic_power_w, frame.mouth_acoustic_power_w, rtol=1e-8)
    assert set(frame.bc_mode) == {'velocity'}


@pytest.mark.parametrize('velocity', [0., float('nan'), float('inf'), complex(0,float('inf'))])
def test_invalid_velocity_is_rejected_before_assembly(tmp_path, velocity):
    with pytest.raises(ValueError, match='finite and nonzero'):
        run_simulation(None, None, (200,800), 5, {}, str(tmp_path/'unused.csv'),
                       bc_mode='velocity', inlet_velocity_rms=velocity)
