import numpy as np
from horn_core.parameters import DriverParameters
from horn_analysis.transfer_function import compute_driver_response, compute_driver_operating_point
from horn_analysis.evaluation import evaluate_response
from horn_analysis.scoring import TargetSpec


def test_driver_pressure_matches_independent_motor_force_balance():
    d=DriverParameters('fixture','Test','Motor',200.,6.,5.,.002,.004,.0001,qms=5.,qes=.4)
    f=np.array([100.,500.,2000.]);area=.001
    z=np.array([100+20j,420+100j,600-20j])
    pressure=compute_driver_response(d,f,z.real,z.imag,area,v_g=2.)
    reference=[]
    for frequency,load in zip(f,z):
        omega=2*np.pi*frequency
        mechanical=d.rms_kg_per_s+1j*omega*d.mms_kg+1/(1j*omega*d.cms_m_per_n)
        electrical=d.re_ohm+1j*omega*d.le_h
        # Voltage = Ze I + BL v; force balance BL I = (Zm + Zload) v.
        current,velocity=np.linalg.solve([[electrical,d.bl_tm],[-d.bl_tm,mechanical+load*d.sd_m2**2/area]],[2.,0.])
        reference.append(load*d.sd_m2*velocity/area)
        electrical_power=np.real(2.*np.conj(current))
        losses=d.re_ohm*abs(current)**2 + d.rms_kg_per_s*abs(velocity)**2
        acoustic_power=np.real(load)*abs(velocity*d.sd_m2)**2/area
        np.testing.assert_allclose(electrical_power,losses+acoustic_power,rtol=1e-12)
    np.testing.assert_allclose(pressure,reference,rtol=1e-12)
    np.testing.assert_allclose(compute_driver_response(d,f,z.real,z.imag,area,v_g=4.),2*pressure)


def test_known_band_and_compression_are_hard_constraints():
    d=DriverParameters('fixture','Test','Motor',200.,6.,5.,.002,.004,.0001,qms=5.,qes=.4,usable_f_low_hz=300,usable_f_high_hz=1500)
    target=TargetSpec(500,2000)
    result=evaluate_response(np.array([500,1000,2000]),np.array([90,90,90]),target,d,.0001)
    assert set(result['rejection_reasons']) == {'compression_ratio_exceeded','outside_driver_usable_band'}
    assert not result['model_feasible'] and result['composite_score']==0
    missing=evaluate_response(np.array([500,1000,1500]),np.array([90,90,90]),TargetSpec(500,1500),d,.002)
    assert missing['eligibility_status']=='insufficient_evidence'
    assert 'driver_interface_unverified' in missing['evidence_gaps']


def test_voltage_operating_limits_and_energy_balance():
    d=DriverParameters('fixture','Test','Motor',100.,6.,5.,.002,.004,.0001,
                       qms=5.,qes=.4,xmax_m=.001,power_w=5.)
    f=np.array([100.,200.,400.]);z=np.full(3,420.)
    quiet=compute_driver_operating_point(d,f,z,z*0,.002,1.)
    loud=compute_driver_operating_point(d,f,z,z*0,.002,100.)
    np.testing.assert_allclose(loud['displacement_peak_m'],100*quiet['displacement_peak_m'])
    np.testing.assert_allclose(loud['input_power_w'],10000*quiet['input_power_w'])
    np.testing.assert_allclose(loud['input_power_w'],loud['copper_power_w']+loud['mechanical_loss_w']+loud['horn_power_w'],rtol=1e-12)
    result=evaluate_response(f,np.full(3,95.),TargetSpec(100,400),d,.002,operating_point=loud)
    assert {'driver_excursion_exceeded','driver_nominal_power_exceeded'} <= set(result['rejection_reasons'])
    assert result['composite_score']==0
    quiet_result=evaluate_response(f,np.full(3,95.),TargetSpec(100,400),d,.002,operating_point=quiet)
    assert quiet_result['model_feasible']
    assert quiet_result['eligibility_status']=='insufficient_evidence'


def test_zero_mechanical_impedance_has_finite_motor_solution():
    d=DriverParameters('fixture','Test','Undamped motor',100.,6.,5.,.002,.004,0.,rms_kg_per_s=0.)
    f=np.array([100.]);z=np.zeros(1)
    point=compute_driver_operating_point(d,f,z,z,.002,2.)
    np.testing.assert_allclose(point['velocity_rms'],[2./5.],atol=1e-12)
    np.testing.assert_allclose(point['current_rms'],[0.],atol=1e-12)
    np.testing.assert_allclose(point['throat_pressure'],[0.])
