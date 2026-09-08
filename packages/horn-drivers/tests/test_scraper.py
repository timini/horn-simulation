"""Offline regressions for failed sources, resume and non-destructive writes."""
import json
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
    assert s.request(session, "url", delay=0, patience_s=3) is None
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
