"""Regression gates for candidate resolution, including unsampled peaks."""
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

PATH=Path(__file__).resolve().parents[1]/'validate_candidate_resolution.py'
spec=importlib.util.spec_from_file_location('candidate_resolution',PATH)
v=importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)


def test_extra_fine_grid_peak_cannot_disappear_through_downsampling():
    coarse=dict(frequency=np.array([1.,3.]),spl=np.zeros(2),z=np.ones(2,dtype=complex))
    fine=dict(frequency=np.array([1.,2.,3.]),spl=np.array([0.,4.,0.]),z=np.array([1.,2j,1.]))
    change=v.curve_change(coarse,fine)
    assert change['spl_db']==4.
    assert change['phase_deg']==90.
    assert change['impedance_db']==pytest.approx(20*np.log10(2))
    assert v.ripple(fine,[1.5,2.5])==2.


def test_partial_band_and_zero_impedance_cannot_pass():
    full=dict(frequency=np.array([1.,3.]),spl=np.zeros(2),z=np.ones(2,dtype=complex))
    partial=dict(full,frequency=np.array([2.,3.]))
    with pytest.raises(ValueError,match='same band'):
        v.curve_change(partial,full)
    with pytest.raises(ValueError,match='not covered'):
        v.ripple(partial,[1.,3.])
    with pytest.raises(ValueError,match='Zero impedance'):
        v.curve_change(full,dict(full,z=np.zeros(2)))


def frame():
    return pd.DataFrame(dict(frequency=[1.,2.],spl=[0.,0.],z_real=[1.,1.],
        schema_version=2,bc_mode='dirichlet',phasor_convention='exp(+iwt)_rms',
        radiation_model='flanged_piston',loss_model='lossless',element_degree=1,
        input_acoustic_power_w=1.,mouth_acoustic_power_w=1.,viscous_wall_power_w=0.,
        thermal_wall_power_w=0.,relative_residual=1e-12,converged_reason=4,mesh_cells=100))


@pytest.mark.parametrize('column,value', [('frequency',1.),('spl',np.nan),
    ('relative_residual',-1.),('relative_residual',1e-3),('converged_reason',0),
    ('mouth_acoustic_power_w',.8),('bc_mode','velocity'),('element_degree',2)])
def test_failed_health_or_changed_contract_cannot_pass(column,value):
    good=frame()
    assert v.check_frame(good,[1.,2.],2)['mesh_cells']==100
    good[column]=value
    with pytest.raises(ValueError):
        v.check_frame(good,[1.,2.],2)


def test_preparation_freezes_candidate_driver_source_and_limits(tmp_path):
    candidate=dict(driver_id='test',loss_model='lossless',radiation_model='flanged_piston',
        element_degree=1,throat_radius=.01,mouth_radius=.03,length=.1,
        drive_voltage_rms=2.83,observation_distance_m=1.,profile='conical')
    ranking=tmp_path/'input-ranking.json'
    driver=tmp_path/'input-driver.json'
    ranking.write_text(json.dumps([candidate]))
    driver.write_text(json.dumps(dict(driver_id='wrong')))
    out=tmp_path/'study'
    with pytest.raises(ValueError,match='does not match'):
        v.prepare(out,ranking,driver,800.,1600.,0)
    assert not out.exists()
    driver.write_text(json.dumps(dict(driver_id='test')))
    v.prepare(out,ranking,driver,800.,1600.,0)
    assert v.verify(out)['candidate']==candidate
    driver_copy=out/'driver.json'
    driver_copy.write_text('{}')
    with pytest.raises(ValueError,match='Frozen input'):
        v.verify(out)
