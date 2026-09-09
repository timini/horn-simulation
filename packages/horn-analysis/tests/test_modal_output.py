"""Nonuniform mouth flow must survive report-side driver coupling."""
from types import SimpleNamespace
import numpy as np
import pandas as pd
import pytest
from horn_analysis import evaluation
from horn_core.acoustics import baffled_piston_on_axis,pressure_level


def modal_frame():
    frame=pd.DataFrame(dict(frequency=[400.,800.],spl=80.,z_real=400.,z_imag=20.,
        schema_version=2,phasor_convention='exp(+iwt)_rms',bc_mode='dirichlet',
        inlet_area_m2=.001,mouth_area_m2=.01,mouth_u_real=.0001,mouth_u_imag=0.,
        radiation_model='modal_baffled',modal_mode_count=16))
    for n in range(16):
        frame[f'modal_v_{n}_real']=.01 if n==0 else 0.
        frame[f'modal_v_{n}_imag']=0.
    return frame


def test_modal_uniform_limit_and_voltage_scaling(monkeypatch):
    monkeypatch.setattr(evaluation,'compute_driver_response',lambda driver,f,zr,zi,area,voltage:np.full(len(f),voltage,dtype=complex))
    frame=modal_frame();target=SimpleNamespace(voltage_rms=2.,observation_distance_m=1.)
    actual,area,metric=evaluation.coupled_output(frame,None,target)
    expected=pressure_level(baffled_piston_on_axis(frame.frequency,np.full(2,.0002),.01,1.))
    np.testing.assert_allclose(actual,expected,atol=1e-10)
    assert area==.001 and metric=='modal_baffled_on_axis'
    target.voltage_rms=4.
    doubled,_,_=evaluation.coupled_output(frame,None,target)
    np.testing.assert_allclose(doubled-actual,20*np.log10(2),atol=1e-12)


def test_zero_total_flow_higher_mode_still_radiates(monkeypatch):
    monkeypatch.setattr(evaluation,'compute_driver_response',lambda driver,f,zr,zi,area,voltage:np.ones(len(f),complex))
    frame=modal_frame();frame['mouth_u_real']=0.;frame['modal_v_0_real']=0.;frame['modal_v_1_real']=.01
    actual,_,_=evaluation.coupled_output(frame,None,SimpleNamespace(voltage_rms=1.,observation_distance_m=.1))
    assert np.isfinite(actual).all() and (actual>0).all()


@pytest.mark.parametrize('failure',['missing','count','mixed','flow','nonfinite'])
def test_incomplete_or_inconsistent_modal_contract_is_rejected(monkeypatch,failure):
    monkeypatch.setattr(evaluation,'compute_driver_response',lambda driver,f,zr,zi,area,voltage:np.ones(len(f),complex))
    frame=modal_frame()
    if failure=='missing':frame=frame.drop(columns=['modal_v_15_imag'])
    elif failure=='count':frame['modal_mode_count']=8
    elif failure=='mixed':frame.loc[0,'radiation_model']='flanged_piston'
    elif failure=='flow':frame['mouth_u_real']=.0002
    else:frame['modal_v_2_real']=np.nan
    with pytest.raises(ValueError,match='modal aperture'):
        evaluation.coupled_output(frame,None,SimpleNamespace(voltage_rms=1.,observation_distance_m=1.))
