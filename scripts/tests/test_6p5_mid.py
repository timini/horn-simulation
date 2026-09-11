import sys
from pathlib import Path
import json
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import study_6p5_mid as s


def test_shared_catalogue_exact_size_and_provenance():
    records,drivers=s.load_candidates()
    assert len(records)>20
    assert 'beyma-6mcf200nd' in drivers and 'bc-6mdn44-8' in drivers
    assert not any(r['manufacturer'] in ('PRV','SB-Audience') for r in records)
    assert all(r['catalogue_status']=='manufacturer_verified' for r in records)


def test_sealed_back_does_not_receive_second_rear_compliance():
    _,drivers=s.load_candidates();d=drivers['beyma-6mcf200nd']
    g=s.core.Geometry('conical',.052,.210,.2);f=np.array([500.,1000.,6500.])
    s.SEALED_IDS.add(d.driver_id)
    try:
        result=s.response(g,f,driver=d,rear_litres=.5)
        expected=s.core.response(g,f,driver=d,rear_litres=None)
        assert np.allclose(result['spl'],expected['spl'],rtol=0,atol=1e-10)
    finally:s.SEALED_IDS.clear()


def test_band_metrics_include_actual_6500_endpoint():
    f=np.array([300,500,1000,6000,6500,10000.])
    result=dict(spl=np.array([0,100,101,102,90,0.]),x_mm=np.ones(6),throat_speed=np.ones(6))
    m=s.metrics(f,result)
    assert m['spl_6500_db']==90 and m['ripple_db']==12


def test_annular_clearance_accounts_for_tweeter_obstruction():
    g=s.core.Geometry('conical',.042,.190,.125)
    r=s.load_candidates()[1]['celestion-cf0617m']
    dims=g.dimensions(r)
    assert dims['outer_throat_diameter_mm']==pytest.approx(np.hypot(42,64))
    assert dims['compression_ratio']==pytest.approx(r.sd_m2/(np.pi*.042**2/4))
