"""Physical limiting cases for the isolated T90A reduced-order study."""
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

path = Path(__file__).resolve().parents[1] / 'study_t90a_annular.py'
spec = importlib.util.spec_from_file_location('t90a_study', path)
study = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = study
spec.loader.exec_module(study)


def test_annular_open_area_retains_tweeter_clearance():
    g = study.Geometry('conical', .047, .18, .25)
    actual = np.pi*(g.outer_throat_radius**2-(.064/2)**2)
    assert actual == pytest.approx(np.pi*(.047/2)**2)
    assert 2*g.outer_throat_radius > .064 > .060
    assert g.dimensions()['radial_gap_mm'] == pytest.approx(7.7020150622)


def test_zero_body_and_chamber_recovers_existing_driver_horn_chain():
    g=study.Geometry('os',.076,.24480268,.3,0.,0.)
    f=np.geomspace(250,10000,71)
    reference=study.compute_horn_transfer_tmm(f,g.radius,g.length_m,
                 g.outer_throat_radius,g.mouth_diameter_m/2,n_segments=300)
    point=study.compute_driver_operating_point(study.BASE_DRIVER,f,
                 reference['z_real'],reference['z_imag'],g.throat_area,2.83)
    uout=reference['mouth_volume_velocity_transfer']*point['throat_pressure']
    pressure=study.baffled_piston_on_axis(f,uout,reference['mouth_area_m2'],1.)
    actual=study.response(g,f,front_cc=0,rear_litres=None)
    np.testing.assert_allclose(actual['pressure'],pressure,rtol=1e-10,atol=1e-12)


def test_front_compliance_conserves_real_power_and_affects_upper_band():
    g=study.Geometry('conical',.060,.18,.25)
    f=np.geomspace(320,10000,101)
    small=study.response(g,f,front_cc=10)
    large=study.response(g,f,front_cc=60)
    for result in (small,large):
        np.testing.assert_allclose(result['radiation_power'],result['horn_power'],rtol=1e-10,atol=1e-14)
        assert np.all(result['radiation_power']>=0)
    assert abs(large['spl'][-1]-small['spl'][-1])>3


def test_centre_body_cannot_block_or_extend_through_mouth():
    with pytest.raises(ValueError):
        study.transfer(study.Geometry('conical',.047,.06,.25),np.array([1000.]))
    with pytest.raises(ValueError):
        study.transfer(study.Geometry('conical',.047,.18,.1),np.array([1000.]))


def test_full_numerical_checks():
    g=study.Geometry('conical',.060,.18,.25)
    f=np.unique(np.r_[np.geomspace(226.274,12000,601),320,7000])
    result=study.checks(g,f)
    assert result['segments_300_to_600_max_spl_db']<.25
