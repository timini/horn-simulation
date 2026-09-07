import numpy as np
import pandas as pd
import pytest
from horn_core.acoustics import baffled_piston_on_axis, inlet_area_from_frame, pressure_level, C0, RHO0
from horn_core.webster import compute_horn_transfer_tmm
from horn_analysis.evaluation import evaluate_response
from horn_analysis.scoring import TargetSpec
from horn_analysis.merge import merge_bands


def test_uniform_tube_has_analytic_complex_transfer_and_positive_load():
    f = np.geomspace(100, 2000, 30)
    r, length = .025, .2
    t = compute_horn_transfer_tmm(f, lambda z: r, length, r, r, radiation_model="plane_wave")
    expected = np.exp(-1j*2*np.pi*f*length/C0)
    np.testing.assert_allclose(t["mouth_pressure_transfer"], expected, atol=1e-12)
    np.testing.assert_allclose(t["z_real"], RHO0*C0, rtol=1e-12)
    np.testing.assert_allclose(t["z_imag"], 0, atol=1e-10)


def test_radiated_pressure_matches_far_field_small_aperture_limit():
    f = np.geomspace(100, 1000, 10)
    u = np.full(10, 1e-5)
    p = baffled_piston_on_axis(f, u, 1e-5, 5)
    reference = 1j*RHO0*f*u/5*np.exp(-1j*2*np.pi*f*5/C0)
    np.testing.assert_allclose(p, reference, rtol=1e-5)
    np.testing.assert_allclose(pressure_level(2*p)-pressure_level(p), 6.020599913, atol=1e-8)


def test_annular_mesh_area_overrides_nominal_outer_radius():
    annulus = np.pi*(.038**2-.035**2)
    frame = pd.DataFrame({"inlet_area_m2": [annulus]*3})
    assert inlet_area_from_frame(frame, .038) == pytest.approx(annulus)
    frame.loc[1, "inlet_area_m2"] *= 2
    with pytest.raises(ValueError, match="Inconsistent"):
        inlet_area_from_frame(frame, .038)


def test_internal_notch_cannot_be_hidden_by_outer_cutoffs():
    target = TargetSpec(100, 1000)
    result = evaluate_response(np.array([100, 250, 500, 750, 1000]), np.array([90,90,60,90,90]), target)
    assert not result["model_feasible"]
    assert result["composite_score"] == 0
    assert result["bandwidth_coverage"] < 1
    assert result["passband_ripple_db"] == 30


def test_outside_target_peak_does_not_determine_ripple():
    result = evaluate_response(np.array([50,100,200,1000,2000]), np.array([130,90,90,90,140]), TargetSpec(100,1000))
    assert result["model_feasible"]
    assert result["passband_ripple_db"] == 0


def test_incomplete_band_is_an_error():
    with pytest.raises(ValueError, match="bracketed"):
        evaluate_response(np.array([200,500,1000]), np.array([90,90,90]), TargetSpec(100,1000))


def write_band(tmp_path, index, low, high):
    f = np.geomspace(low,high,3)
    frame = pd.DataFrame({"frequency":f,"spl":94.,"z_real":420.,"z_imag":0.,
        "schema_version":2,"inlet_area_m2":.01,"mouth_area_m2":.02,
        "mouth_u_real":.0001,"mouth_u_imag":0.,"mouth_p_real":1.,"mouth_p_imag":0.,"radiation_model":"plane_wave",
        "phasor_convention":"exp(+iwt)_rms","bc_mode":"dirichlet"})
    p = tmp_path/f"results_candidate_{index}.csv"
    frame.to_csv(p,index=False)
    return p


def test_merger_rejects_missing_band_and_deduplicates_valid_endpoints(tmp_path):
    p0 = write_band(tmp_path,0,100,550)
    args = dict(num_bands=2,min_freq=100,max_freq=1000,points_per_band=3,output=tmp_path/'out.csv')
    with pytest.raises(ValueError, match="Incomplete simulation"):
        merge_bands([p0],**args)
    p1 = write_band(tmp_path,1,550,1000)
    result = merge_bands([p1,p0],**args)
    assert len(result) == 5
    assert np.all(np.diff(result.frequency)>0)
    df = pd.read_csv(p1);df.loc[0,'spl']=96;df.to_csv(p1,index=False)
    with pytest.raises(ValueError, match="mismatch"):
        merge_bands([p0,p1],**args)


@pytest.mark.parametrize('low,high', [(0,100), (500,100), (100,np.inf), (np.nan,1000)])
def test_bad_targets_fail_before_simulation(low,high):
    with pytest.raises(ValueError):
        TargetSpec(low,high)


@pytest.mark.parametrize('prefix,value', [('z',-420+0j),('mouth_p',1j),('mouth_u',-.0001+0j)])
def test_overlap_rejects_complex_disagreement_even_when_spl_matches(tmp_path,prefix,value):
    p0=write_band(tmp_path,0,100,550);p1=write_band(tmp_path,1,550,1000)
    frame=pd.read_csv(p1)
    frame.loc[0,prefix+'_real']=value.real;frame.loc[0,prefix+'_imag']=value.imag
    frame.to_csv(p1,index=False)
    with pytest.raises(ValueError,match='complex '+prefix+' mismatch'):
        merge_bands([p0,p1],num_bands=2,min_freq=100,max_freq=1000,points_per_band=3,output=tmp_path/'out.csv')
    assert not (tmp_path/'out.csv').exists()


def test_overlap_allows_small_complex_mesh_error_near_zero(tmp_path):
    p0=write_band(tmp_path,0,100,550);p1=write_band(tmp_path,1,550,1000)
    frame=pd.read_csv(p1);frame.loc[0,'z_imag']=.01;frame.loc[0,'mouth_u_imag']=1e-9;frame.to_csv(p1,index=False)
    assert len(merge_bands([p0,p1],num_bands=2,min_freq=100,max_freq=1000,points_per_band=3,output=tmp_path/'out.csv'))==5


def test_single_coupling_rejects_already_driven_neumann_results(tmp_path,monkeypatch):
    import horn_analysis.couple_single as single
    from horn_core.parameters import DriverParameters
    d=DriverParameters('test','Test','Motor',200.,6.,5.,.002,.004,.0001,qms=5.,qes=.4)
    monkeypatch.setattr(single,'load_driver',lambda *_:d)
    p=write_band(tmp_path,0,100,1000);frame=pd.read_csv(p);frame['bc_mode']='neumann';frame.to_csv(p,index=False)
    with pytest.raises(ValueError,match='Dirichlet transfer'):
        single.couple(str(p),'unused','test',.025,output_csv=str(tmp_path/'coupled.csv'))
    assert not (tmp_path/'coupled.csv').exists()
