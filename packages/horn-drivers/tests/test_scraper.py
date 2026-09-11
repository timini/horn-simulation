"""Offline regressions for failed sources, resume and non-destructive writes."""
import json
import os
import stat
from types import SimpleNamespace
import pytest
from horn_drivers import scraper as s


class Session:
    def __init__(self, *responses):
        self.responses = list(responses)
        self.calls = []

    def get(self, url, **kwargs):
        self.calls.append(url)
        response = self.responses.pop(0)
        if isinstance(response, Exception):
            raise response
        return response


def response(status=200, text="", headers=None):
    return SimpleNamespace(status_code=status, text=text, headers=headers or {})


@pytest.fixture
def clock(monkeypatch):
    now = [0.]
    monkeypatch.setattr(s.time, "monotonic", lambda: now[0])
    monkeypatch.setattr(s.time, "sleep", lambda wait: now.__setitem__(0, now[0]+wait))
    monkeypatch.setattr(s, "_polite_sleep", lambda *args: None)
    return now


@pytest.mark.parametrize("status", [429, 503, 500, 520])
def test_failed_requests_have_finite_attempts(clock, status):
    session = Session(*[response(status)]*3)
    assert s.request(session, "url", delay=0, max_attempts=3) is None
    assert len(session.calls) == 3


def test_origin_patience_stops_at_deadline(clock):
    session = Session(*[response(522)]*10)
    with pytest.raises(s.OriginUnavailable):
        s.request(session, "url", delay=0, patience_s=3)
    assert clock[0] == 3
    assert len(session.calls) == 2


def test_retry_after_is_honored_or_run_stops(clock):
    session = Session(response(429, headers={"Retry-After":"5"}), response())
    assert s.request(session, "url", delay=0).status_code == 200
    assert clock[0] == 5
    session = Session(response(429, headers={"Retry-After":"3600"}))
    assert s.request(session, "url", delay=0, max_backoff=300) is None
    assert len(session.calls) == 1


def test_session_only_advertises_installed_compression():
    session = s.build_session()
    assert session.headers['Accept-Encoding'] == 'gzip, deflate'
    session.close()


@pytest.mark.parametrize("probe_status", [200, 403, 429, 503])
def test_origin_probe_requires_rejection_of_unknown_routes(clock, probe_status):
    assert not s.origin_is_healthy(Session(response(), response(probe_status)))


def test_origin_probe_accepts_normal_missing_route(clock):
    assert s.origin_is_healthy(Session(response(), response(404)))


RAW = {"fs":100, "re":5, "bl":8, "sd":100, "mmd":10, "pmax":200}
PARAMS = {"fs_hz":100., "re_ohm":5., "bl_tm":8., "sd_m2":.01, "mms_kg":.012, "mmd_kg":.01, "le_h":.0003}
ENTRY = {"manufacturer":"Test", "name":"Model", "url":s.BASE_URL+'/Test/Model'}


def html_page(identity='/Test/Model', with_mms=True):
    meta = f'<meta property="og:url" content="{s.BASE_URL}{identity}">' if identity else ''
    mms = '<span data-highlight="mms"></span><span data-highlight="mms"><b>12</b></span>' if with_mms else ''
    return meta + '<div data-woofer=\'' + json.dumps(RAW) + '\'></div>' + mms


def test_mass_and_power_are_not_invented_from_different_quantities(clock):
    parsed = s._parse_data_woofer(json.dumps(RAW))
    assert parsed['mmd_kg'] == .01
    assert 'mms_kg' not in parsed and 'power_w' not in parsed
    parsed = s.scrape_driver_page(ENTRY['url'], Session(response(text=html_page())))
    assert parsed['mms_kg'] == .012 and parsed['mmd_kg'] == .01
    assert parsed['peak_power_w'] == 200 and 'power_w' not in parsed


@pytest.mark.parametrize("identity", [None, '/Wrong/Driver'])
def test_unidentified_or_wrong_driver_page_is_rejected(clock, identity):
    with pytest.raises(s.OriginWedged):
        s.scrape_driver_page(ENTRY['url'], Session(response(text=html_page(identity))))


@pytest.mark.parametrize("raw", [[], {**RAW, 'mmd':float('inf')}, {**RAW, 'fs':'bad'}])
def test_bad_numeric_source_is_rejected(raw):
    assert s._parse_data_woofer(json.dumps(raw)) is None


def test_failed_discovery_does_not_return_successful_partial_list(clock):
    with pytest.raises(s.ScrapeError):
        s.discover_manufacturers(Session(response(404)))
    first = '<script class="count">41</script><a href="/Test/Model">Model</a>'
    with pytest.raises(s.ScrapeError, match='pagination'):
        s.discover_drivers(Session(response(text=first), response(404)), 'Test', delay=0)


def batch(monkeypatch):
    monkeypatch.setattr(s, 'build_session', lambda: Session())
    monkeypatch.setattr(s, 'discover_drivers', lambda *a, **k: [ENTRY])
    monkeypatch.setattr(s, 'scrape_driver_page', lambda *a, **k: dict(PARAMS))


def test_resume_writes_once_then_skips_valid_record(monkeypatch, tmp_path):
    batch(monkeypatch)
    kwargs = dict(db_dir=tmp_path, manufacturer_filter=['Test'], required_fields=('mmd_kg',))
    assert s.scrape_all(**kwargs) == 1
    path = tmp_path/'Test/test-model.json'
    before = path.read_bytes()
    monkeypatch.setattr(s, 'scrape_driver_page', lambda *a, **k: pytest.fail('Resume must not refetch a current record'))
    assert s.scrape_all(**kwargs) == 0
    assert path.read_bytes() == before
    record = json.loads(before); record['parameters']['mmd_kg'] = None
    path.write_text(json.dumps(record))
    assert not s._driver_is_current(tmp_path, 'Test', 'test-model', ('mmd_kg',))


def test_failed_refresh_preserves_existing_record(monkeypatch, tmp_path):
    batch(monkeypatch)
    kwargs = dict(db_dir=tmp_path, manufacturer_filter=['Test'])
    s.scrape_all(**kwargs)
    path = tmp_path/'Test/test-model.json'; before = path.read_bytes()
    monkeypatch.setattr(s, 'scrape_driver_page', lambda *a, **k: None)
    state = tmp_path/'state.json'
    with pytest.raises(s.ScrapeError, match='1 drivers failed'):
        s.scrape_all(**kwargs, refresh=True, state_path=state)
    assert path.read_bytes() == before
    assert json.loads(state.read_text())['Test']['complete'] is False


def test_interrupted_atomic_write_preserves_existing_file(monkeypatch, tmp_path):
    path = tmp_path/'record.json';path.write_text('old record')
    def fail(*args): raise OSError('interrupted')
    monkeypatch.setattr(s.os, 'replace', fail)
    with pytest.raises(OSError): s._atomic_json(path, {'new':'record'})
    assert path.read_text() == 'old record'
    assert list(tmp_path.iterdir()) == [path]


def test_cli_failed_discovery_and_noop_resume_have_different_exit_codes(monkeypatch, tmp_path):
    monkeypatch.setattr('sys.argv', ['horn-scrape-drivers', '--db', str(tmp_path)])
    def fail(**kwargs): raise s.ScrapeError('discovery failed')
    monkeypatch.setattr(s, 'scrape_all', fail)
    assert s.main() == 1
    monkeypatch.setattr(s, 'scrape_all', lambda **kwargs: 0)
    assert s.main() == 0


def test_database_paths_cannot_escape_output_directory(tmp_path):
    with pytest.raises(ValueError):
        s._save_driver(tmp_path, {'manufacturer':'..','driver_id':'bad'})
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('sd', [.00126, .0045, .005])
def test_refresh_preserves_compression_category_without_cone_diameter(monkeypatch, tmp_path, sd):
    batch(monkeypatch)
    s._save_driver(tmp_path, dict(driver_id='test-model', manufacturer='Test',
        driver_type='compression', parameters={**PARAMS, 'sd_m2': sd}))
    monkeypatch.setattr(s, 'scrape_driver_page', lambda *a, **k: {**PARAMS, 'sd_m2': sd})
    assert s.scrape_all(db_dir=tmp_path, manufacturer_filter=['Test'], refresh=True) == 1
    record = json.loads((tmp_path/'Test/test-model.json').read_text())
    assert record['driver_type'] == 'compression'
    assert 'nominal_diameter' not in record


def test_new_ambiguous_driver_has_no_invented_category_or_diameter(monkeypatch, tmp_path):
    batch(monkeypatch)
    monkeypatch.setattr(s, 'scrape_driver_page', lambda *a, **k: {**PARAMS, 'sd_m2': .00126})
    s.scrape_all(db_dir=tmp_path, manufacturer_filter=['Test'])
    record = json.loads((tmp_path/'Test/test-model.json').read_text())
    assert record['driver_type'] == 'unknown'
    assert 'nominal_diameter' not in record


def test_resume_migrates_provisional_mass_and_missing_provenance(monkeypatch, tmp_path):
    batch(monkeypatch)
    s._save_driver(tmp_path, dict(driver_id='test-model', manufacturer='Test',
        parameters={**PARAMS, 'mms_kg': PARAMS['mmd_kg']}))
    assert not s._driver_is_current(tmp_path, 'Test', 'test-model', ('mmd_kg',))
    assert s.scrape_all(db_dir=tmp_path, manufacturer_filter=['Test']) == 1
    path = tmp_path/'Test/test-model.json'
    record = json.loads(path.read_text())
    assert record['parameters']['mms_kg'] == .012
    assert record['parameter_source'] == ENTRY['url']
    assert s._driver_is_current(tmp_path, 'Test', 'test-model')
    del record['parameter_source']
    path.write_text(json.dumps(record))
    assert not s._driver_is_current(tmp_path, 'Test', 'test-model')


def test_atomic_refresh_preserves_destination_permissions(tmp_path):
    path = tmp_path/'driver.json'
    path.write_text('{}')
    path.chmod(0o640)
    s._atomic_json(path, {'fresh': True})
    assert stat.S_IMODE(path.stat().st_mode) == 0o640


@pytest.mark.parametrize('mask', [0o022, 0o027, 0o077])
def test_new_atomic_file_respects_process_umask(tmp_path, mask):
    path = tmp_path/'driver.json'
    old_mask = os.umask(mask)
    try:
        s._atomic_json(path, {'fresh': True})
    finally:
        os.umask(old_mask)
    assert stat.S_IMODE(path.stat().st_mode) == (0o666 & ~mask)


@pytest.mark.parametrize('refresh', [False, True])
def test_refresh_preserves_enriched_interface_and_parameters(monkeypatch, tmp_path, refresh):
    batch(monkeypatch)
    enriched = dict(driver_id='test-model', manufacturer='Test', driver_type='compression',
        interface_model='measured_adapter', usable_f_low_hz=500, usable_f_high_hz=6000,
        parameter_sources={'power_w': 'https://example.org/manufacturer/continuous-power-specification'},
        notes='Measured by the project owner',
        parameters={**PARAMS, 'fs_hz': 150, 'exit_area_m2': .0005, 'rear_load_mass_kg': .001,
                    'power_w': 80})
    s._save_driver(tmp_path, enriched)
    assert s.scrape_all(db_dir=tmp_path, manufacturer_filter=['Test'], refresh=refresh) == 1
    record = json.loads((tmp_path/'Test/test-model.json').read_text())
    assert record['parameters']['fs_hz'] == PARAMS['fs_hz']
    for field in ('exit_area_m2', 'rear_load_mass_kg', 'power_w'):
        assert record['parameters'][field] == enriched['parameters'][field]
    for field in ('interface_model', 'usable_f_low_hz', 'usable_f_high_hz', 'notes'):
        assert record[field] == enriched[field]


def test_origin_outage_aborts_batch_without_restarting_patience(monkeypatch, tmp_path, clock):
    batch(monkeypatch)
    entries = [ENTRY, {**ENTRY, 'name': 'Second', 'url': s.BASE_URL+'/Test/Second'},
               {**ENTRY, 'name': 'Third', 'url': s.BASE_URL+'/Test/Third'}]
    monkeypatch.setattr(s, 'discover_drivers', lambda *a, **k: entries)
    session = Session(*[response(522)]*10)
    visited = []
    def page(url, unused_session, **kwargs):
        visited.append(url)
        if url == ENTRY['url']:
            return dict(PARAMS)
        return s.request(session, url, delay=0, patience_s=kwargs['patience_s'])
    monkeypatch.setattr(s, 'scrape_driver_page', page)
    state = tmp_path/'state.json'
    with pytest.raises(s.OriginUnavailable):
        s.scrape_all(db_dir=tmp_path, manufacturer_filter=['Test'], patience_s=3, state_path=state)
    assert clock[0] == 3
    assert visited == [entry['url'] for entry in entries[:2]]
    assert (tmp_path/'Test/test-model.json').exists()
    assert not (tmp_path/'Test/test-third.json').exists()
    progress = json.loads(state.read_text())['Test']
    assert progress['complete'] is False and progress['scraped'] == 1


def test_schema_migration_replaces_malformed_parameter_mapping(monkeypatch, tmp_path):
    batch(monkeypatch)
    s._save_driver(tmp_path, dict(driver_id='test-model', manufacturer='Test', parameters=None))
    assert s.scrape_all(db_dir=tmp_path, manufacturer_filter=['Test']) == 1
    assert s._driver_is_current(tmp_path, 'Test', 'test-model')


@pytest.mark.parametrize('old_version', [None, 2])
def test_migration_removes_legacy_inferred_power(monkeypatch, tmp_path, old_version):
    batch(monkeypatch)
    s._save_driver(tmp_path, dict(driver_id='test-model', manufacturer='Test',
        scraper_schema_version=old_version,
        parameters={**PARAMS, 'power_w': 100, 'peak_power_w': 200}))
    assert s.scrape_all(db_dir=tmp_path, manufacturer_filter=['Test']) == 1
    record = json.loads((tmp_path/'Test/test-model.json').read_text())
    assert 'power_w' not in record['parameters']
    assert 'peak_power_w' not in record['parameters']  # Unproven stale rating must not survive
    assert s._driver_is_current(tmp_path, 'Test', 'test-model')


def test_failed_discovery_clears_old_complete_progress(monkeypatch, tmp_path):
    batch(monkeypatch)
    state = tmp_path/'state.json'
    state.write_text(json.dumps({'Test': {'complete': True}}))
    def fail(*args, **kwargs):
        assert json.loads(state.read_text())['Test']['complete'] is False
        raise s.ScrapeError('Incomplete pagination')
    monkeypatch.setattr(s, 'discover_drivers', fail)
    with pytest.raises(s.ScrapeError, match='pagination'):
        s.scrape_all(db_dir=tmp_path, manufacturer_filter=['Test'], refresh=True, state_path=state)
    progress = json.loads(state.read_text())['Test']
    assert progress['complete'] is False
    assert progress['failure_reason'] == 'Incomplete pagination'


def test_secondary_refresh_cannot_overwrite_verified_manufacturer(monkeypatch,tmp_path):
    batch(monkeypatch)
    original=dict(driver_id='test-model',manufacturer='Test',catalogue_status='manufacturer_verified',
                  parameter_source='https://manufacturer.example/model',parameters={**PARAMS,'fs_hz':123})
    s._save_driver(tmp_path,original)
    assert s.scrape_all(db_dir=tmp_path,manufacturer_filter=['Test'],refresh=True)==0
    assert json.loads((tmp_path/'Test/test-model.json').read_text())==original


def test_refresh_missing_inductance_cannot_complete(monkeypatch,tmp_path):
    batch(monkeypatch)
    incomplete={k:v for k,v in PARAMS.items() if k!='le_h'}
    monkeypatch.setattr(s,'scrape_driver_page',lambda *a,**kw:incomplete)
    with pytest.raises(s.ScrapeError,match='failed'):
        s.scrape_all(db_dir=tmp_path,manufacturer_filter=['Test'])
    assert not (tmp_path/'Test/test-model.json').exists()


def test_zero_inductance_can_complete_and_resume(monkeypatch,tmp_path):
    batch(monkeypatch)
    monkeypatch.setattr(s,'scrape_driver_page',lambda *a,**kw:{**PARAMS,'le_h':0})
    assert s.scrape_all(db_dir=tmp_path,manufacturer_filter=['Test'])==1
    assert s._driver_is_current(tmp_path,'Test','test-model')
