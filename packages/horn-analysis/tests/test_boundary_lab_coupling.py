"""Regression against complex motor outputs from a pinned independent solver.

No Julia dependency in normal CI. The recorded loads and expected motor values
come from independently coupled Boundary Lab solves, not this implementation.
"""
import json
from pathlib import Path
import numpy as np
import pytest
from horn_core.parameters import DriverParameters
from horn_analysis.transfer_function import compute_driver_operating_point

EVIDENCE = json.loads((Path(__file__).resolve().parents[3]/'data/validation/boundary_lab_reference.json').read_text())


@pytest.mark.parametrize('case', EVIDENCE['cases'], ids=lambda row: row['case'])
def test_motor_matches_independent_coupled_fem(case):
    assert EVIDENCE['passed'] and case['passed']
    protocol = EVIDENCE['definition']
    assert protocol['reference_kind'] == 'independent_numerical_not_physical'
    frequencies = np.geomspace(protocol['frequency_min_hz'], protocol['frequency_max_hz'], protocol['frequency_count'])
    assert case['frequencies'] == len(frequencies)
    fields = {key: np.array(value)[:,0]+1j*np.array(value)[:,1]
              for key,value in case['reference_values'].items()}
    for values in fields.values():
        assert values.shape == frequencies.shape and np.isfinite(values).all()
    d = protocol['driver']
    area = case['mesh_area_m2']
    driver = DriverParameters(driver_id='boundary-lab-reference', manufacturer='Synthetic',
        model_name='Independent characterized piston',
        fs_hz=1/(2*np.pi*np.sqrt(d['mmd_kg']*d['cms_m_per_n'])),
        re_ohm=d['re_ohm'], le_h=d['le_h'], bl_tm=d['bl_n_per_a'], sd_m2=area,
        mms_kg=d['mmd_kg'], mmd_kg=d['mmd_kg'], rear_load_mass_kg=0.,
        cms_m_per_n=d['cms_m_per_n'], rms_kg_per_s=d['rms_n_s_per_m'])
    actual = compute_driver_operating_point(driver, frequencies, fields['z'].real, fields['z'].imag,
                                             area, protocol['voltage_rms'])
    for name, reference in [('current_rms','current'),('velocity_rms','velocity'),('throat_pressure','inlet_pressure')]:
        np.testing.assert_allclose(actual[name],fields[reference],rtol=protocol['relative_complex_error_limit'],atol=0)
    np.testing.assert_allclose(actual['electrical_impedance'], protocol['voltage_rms']/fields['current'],
                               rtol=protocol['relative_complex_error_limit'],atol=0)
