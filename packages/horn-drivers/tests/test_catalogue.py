import copy
import json
from pathlib import Path
import pytest
from horn_drivers.catalogue import audit_records, redcatt_record
from horn_drivers.loader import _driver_from_dict, load_drivers_raw
from horn_drivers.validator import validate_driver
from horn_drivers import scraper

ROOT=Path(__file__).resolve().parents[3]


def record():
    return json.loads((ROOT/'data/drivers/Celestion/celestion-cf0617m.json').read_text())


def test_repeated_source_quarantine_and_recovered_identity():
    rows=[]
    for i in range(12):
        r=record();r['driver_id']=f'wrong-{i}';r.pop('catalogue_status');rows.append(r)
    groups=audit_records(rows)
    assert groups[0]['count']==12
    assert all(r['catalogue_status']=='quarantined' for r in rows)
    with pytest.raises(ValueError,match='quarantined'):_driver_from_dict(rows[0])
    recovered=record();rows[0].update(recovered)
    audit_records(rows)
    assert _driver_from_dict(rows[0]).fs_hz==116.6


def test_small_family_not_automatically_quarantined():
    a=record();a.pop('catalogue_status');b=copy.deepcopy(a);b['driver_id']='family-variant'
    assert audit_records([a,b])==[]
    assert a['catalogue_status']=='legacy_unverified'


def test_catalogue_preserves_quarantines_and_corrects_6mdn44():
    rows=load_drivers_raw(str(ROOT/'data/drivers'))
    quarantined=[r for r in rows if r.get('catalogue_status')=='quarantined']
    assert len(quarantined)>=375
    for r in quarantined:
        with pytest.raises(ValueError,match='quarantined'):_driver_from_dict(r)
    bc=next(r for r in rows if r['driver_id']=='bc-6mdn44-8')
    d=_driver_from_dict(bc)
    assert d.nominal_diameter=='6.5in' and d.fs_hz==140 and d.sd_m2==pytest.approx(.0132)
    for r in rows:
        if r.get('catalogue_status')=='manufacturer_verified':
            assert r['parameter_source'].startswith('https://')
            assert not validate_driver(r).errors, (r['driver_id'], validate_driver(r).errors)


def test_unverified_size_not_advertised_by_loader():
    r=record();r['nominal_diameter_verified']=False
    assert _driver_from_dict(r).nominal_diameter is None


@pytest.mark.parametrize('bad',[float('nan'),float('inf'),-1,'6.2',True])
def test_invalid_numeric_validation_does_not_crash(bad):
    r=record();r['parameters']['re_ohm']=bad
    assert validate_driver(r).errors


def test_zero_inductance_valid_but_missing_is_not():
    r=record();r['parameters']['le_h']=0
    assert not validate_driver(r).errors
    r['parameters'].pop('le_h')
    assert validate_driver(r).errors
    with pytest.raises(ValueError,match='inductance'):_driver_from_dict(r)


def test_optional_scraper_values_cannot_poison_catalogue():
    base=dict(fs=100,re=6,bl=10,sd=140,mmd=12,le=0)
    assert scraper._parse_data_woofer(json.dumps(base))['le_h']==0
    for field in ['le','qts','z','xmax','pmax']:
        assert scraper._parse_data_woofer(json.dumps({**base,field:'NaN'})) is None


def test_recommended_product_is_not_primary_driver():
    from types import SimpleNamespace
    base=dict(fs=100,re=6,bl=10,sd=140,mmd=12,le=.3)
    html="<meta property='og:url' content='/Test/Model'>"
    html+="<article class='woofer_card' data-woofer='"+json.dumps({**base,'fs':50})+"'></article>"
    html+="<div data-graph-size='normal' data-woofer='"+json.dumps(base)+"'></div>"
    class Session:
        def get(self,*args,**kwargs):return SimpleNamespace(status_code=200,text=html,headers={})
    assert scraper.scrape_driver_page(scraper.BASE_URL+'/Test/Model',Session())['fs_hz']==100


def test_redcatt_units_identity_nulls_and_dual_ratings():
    p=json.loads((Path(__file__).parent/'fixtures/redcatt-6npm.json').read_text())
    r=redcatt_record(p,'2026-09-11')
    d=_driver_from_dict(r)
    assert d.sd_m2==pytest.approx(.01431) and d.le_h==pytest.approx(.00016)
    assert d.mms_kg==pytest.approx(.0111) and d.power_w==180
    assert r['ordering_code']=='DR-6.5-006-8R-B2-C0001'
    assert r['usable_f_high_hz'] is None
    assert r.get('overall_diameter_m') is None  # 162 < 172 bolt circle, ears unmeasured
    p['power_handling']['aes_watts']='600/80'
    with pytest.raises(ValueError):redcatt_record(p,'2026-09-11')
    p['category']='Coax'
    with pytest.raises(ValueError):redcatt_record(p,'2026-09-11')


@pytest.mark.parametrize('value,expected',[('1.23 mH',1.23),('1,23',1.23),('1.2e-3 H',.0012),('1 to 2',None),('NaN',None)])
def test_single_numeric_scraper_value(value,expected):
    assert scraper._parse_float(value)==expected
