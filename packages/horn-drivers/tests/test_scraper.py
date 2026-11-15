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
PARAMS = {"fs_hz":100., "re_ohm":5., "bl_tm":8., "sd_m2":.01, "mms_kg":.012, "mmd_kg":.01}
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
def test_refresh_preserves_sourced_compression_category_without_cone_diameter(monkeypatch, tmp_path, sd):
    batch(monkeypatch)
    s._save_driver(tmp_path, dict(driver_id='test-model', manufacturer='Test',
        driver_type='compression', type_source='datasheet: https://example.org/ds.pdf',
        parameters={**PARAMS, 'sd_m2': sd}))
    monkeypatch.setattr(s, 'scrape_driver_page', lambda *a, **k: {**PARAMS, 'sd_m2': sd})
    assert s.scrape_all(db_dir=tmp_path, manufacturer_filter=['Test'], refresh=True) == 1
    record = json.loads((tmp_path/'Test/test-model.json').read_text())
    assert record['driver_type'] == 'compression'
    assert record['type_source'] == 'datasheet: https://example.org/ds.pdf'
    assert 'nominal_diameter' not in record


def test_refresh_replaces_unsourced_heuristic_category_with_site_category(monkeypatch, tmp_path):
    """Old records got 'compression' from Sd < 0.005; a small full-range cone
    must not keep that label through a refresh."""
    batch(monkeypatch)
    s._save_driver(tmp_path, dict(driver_id='test-model', manufacturer='Test',
        driver_type='compression', parameters={**PARAMS, 'sd_m2': .003}))
    monkeypatch.setattr(s, 'scrape_driver_page',
                        lambda *a, **k: {**PARAMS, 'sd_m2': .003, '_category': 'Full-range'})
    s.scrape_all(db_dir=tmp_path, manufacturer_filter=['Test'], refresh=True)
    record = json.loads((tmp_path/'Test/test-model.json').read_text())
    assert record['driver_type'] == 'cone'
    assert record['type_source'] == "loudspeakerdatabase.com category 'Full-range'"
    monkeypatch.setattr(s, 'scrape_driver_page', lambda *a, **k: {**PARAMS, 'sd_m2': .003})
    s._save_driver(tmp_path, {**record, 'driver_type': 'compression', 'type_source': None})
    s.scrape_all(db_dir=tmp_path, manufacturer_filter=['Test'], refresh=True)
    assert json.loads((tmp_path/'Test/test-model.json').read_text())['driver_type'] == 'unknown'


SITE_JSON = {"re": 3.3, "bl": 10.2, "mmd": 42.0, "rms": 2.7, "cms": 160.0, "sd": 549.0, "le": 0.44,
             "fs": 55.0, "qts": 0.52, "xmax": 5.42, "pmax": 400, "frmin": 55, "frmax": 4000,
             "spl1w": 94.4, "z": 4}
SIMILAR_JSON = {**SITE_JSON, "re": 5.5, "fs": 99.0, "frmin": 20, "frmax": 900}


def real_page(identity='/Test/Model_4\u03a9'):
    """Shape of a real driver page: own data first, 'similar drivers' after."""
    def hl(key, value):
        return (f'<li><span data-highlight="{key}">label</span>'
                f'<span class="symbol_value" data-highlight="{key}"><b>{value}</b></span></li>')
    return (f'<html><head><meta charset="utf-8"><meta property="og:url" content="{s.BASE_URL}{identity}"></head>'
            '<body><main><header><a class="type-ribbon" href="/Subwoofer"><h2>Subwoofer</h2></a></header>'
            '<table><tr><td class="value" data-highlight="size">12\u2033 / 300 mm</td></tr></table>'
            f"<div data-woofer='{json.dumps(SITE_JSON)}'></div><ul>"
            + hl('qms', '6.60') + hl('qes', '0.56') + hl('mms', '51') + hl('pw', '200') +
            '</ul></main><section class="similar">'
            f"<article data-woofer='{json.dumps(SIMILAR_JSON)}'></article>"
            + hl('pw', '999') + '</section></body></html>')


def test_driver_page_maps_band_power_damping_and_sensitivity(clock):
    """The site supplies Rms, Cms, band, power and sensitivity; the old parser
    threw them away, leaving every scraped driver with unknown band and power."""
    page = response(text=real_page())
    parsed = s.scrape_driver_page(s.BASE_URL + '/Test/Model_4\u03a9', Session(page))
    assert parsed['fs_hz'] == 55.0 and parsed['re_ohm'] == 3.3  # own block, not the similar card
    assert parsed['rms_kg_per_s'] == 2.7
    assert parsed['cms_m_per_n'] == pytest.approx(160e-6)
    assert parsed['sd_m2'] == 0.0549  # no binary noise from the unit conversion
    assert (parsed['_usable_f_low_hz'], parsed['_usable_f_high_hz']) == (55.0, 4000.0)
    assert parsed['sensitivity_db_1w1m'] == 94.4
    assert parsed['power_w'] == 200 and parsed['peak_power_w'] == 400
    assert parsed['qms'] == 6.6 and parsed['mms_kg'] == pytest.approx(0.051)
    assert parsed['_category'] == 'Subwoofer' and parsed['_nominal_diameter'] == '12in'


def test_page_without_charset_is_decoded_as_utf8():
    """requests falls back to Latin-1 without a charset, turning the ohm sign
    into mojibake and failing the identity check on every page."""
    class Resp:
        status_code = 200
        headers = {'Content-Type': 'text/html'}
        encoding = 'ISO-8859-1'
    resp = Resp()
    s._fix_encoding(resp)
    assert resp.encoding == 'utf-8'
    resp = Resp(); resp.headers = {'Content-Type': 'text/html; charset=windows-1252'}
    s._fix_encoding(resp)
    assert resp.encoding == 'ISO-8859-1'


def test_scraped_record_carries_band_provenance_and_schema(monkeypatch, tmp_path):
    batch(monkeypatch)
    monkeypatch.setattr(s, 'scrape_driver_page', lambda *a, **k: {
        **PARAMS, 'rms_kg_per_s': 2.7, '_usable_f_low_hz': 55.0, '_usable_f_high_hz': 4000.0,
        '_category': 'Woofer', '_nominal_diameter': '12in'})
    s.scrape_all(db_dir=tmp_path, manufacturer_filter=['Test'])
    record = json.loads((tmp_path/'Test/test-model.json').read_text())
    assert (record['usable_f_low_hz'], record['usable_f_high_hz']) == (55.0, 4000.0)
    assert record['driver_type'] == 'cone' and record['nominal_diameter'] == '12in'
    assert record['parameter_source'] == ENTRY['url']
    assert record['scraper_schema_version'] == s.SCRAPER_SCHEMA_VERSION
    assert record['record_schema_version'] == s.RECORD_SCHEMA_VERSION
    assert record['provenance']['retrieved_on']
    assert 'nominal' in record['provenance']['band_definition']
    assert not any(key.startswith('_') for key in record['parameters'])
    from horn_drivers.loader import _driver_from_dict
    driver = _driver_from_dict(record)
    assert driver.usable_f_low_hz == 55.0 and driver.rms_kg_per_s == 2.7


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
    assert record['parameters']['peak_power_w'] == 200
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


def test_refetch_drops_backfilled_legacy_provenance_notes(monkeypatch, tmp_path):
    batch(monkeypatch)
    s._save_driver(tmp_path, dict(
        driver_id='test-model', manufacturer='Test', model_name='Model', driver_type='unknown',
        parameters={**PARAMS, 'mms_kg': PARAMS['mmd_kg']}, parameter_source=ENTRY['url'],
        record_schema_version=1,
        provenance={'retrieved_on': '2026-02-26', 'retrieved_on_basis': 'commit date',
                    'known_limitations': ['Band, power handling, Rms and sensitivity were not captured.']}))
    assert s.scrape_all(db_dir=tmp_path, manufacturer_filter=['Test']) == 1
    provenance = json.loads((tmp_path/'Test/test-model.json').read_text())['provenance']
    assert 'retrieved_on_basis' not in provenance and 'known_limitations' not in provenance
    assert provenance['retrieved_on'] != '2026-02-26'
