"""The polite fetch layer: robots.txt, pacing, cache, block handling, User-Agent."""

import datetime as dt
from types import SimpleNamespace

import pytest

from horn_drivers import scraper
from horn_drivers.fetch import url_key, DisallowedByRobots, PoliteFetcher, ResponseTooLarge, SourceBlocked

ROBOTS = b"User-agent: *\nDisallow: /private/\n\nUser-agent: horn-simulation-research\nDisallow: /nobots/\nCrawl-delay: 12\n"


class Session:
    def __init__(self, routes):
        self.routes = routes          # url -> list of (status, body, headers) served in order
        self.headers = {}
        self.calls = []

    def get(self, url, timeout=None, stream=False):
        self.calls.append(url)
        queue = self.routes.get(url)
        status, body, headers = (queue.pop(0) if len(queue) > 1 else queue[0]) if queue else (404, b"", {})
        resp = SimpleNamespace(status_code=status, content=body, headers=headers, encoding=None,
                               text=body.decode("utf-8", "replace"))
        if stream:
            resp.iter_content = lambda chunk_size: [body[i:i + chunk_size] for i in range(0, len(body), chunk_size)] or [b""]
            resp.close = lambda: None
        return resp


class Clock:
    def __init__(self):
        self.t = 1000.0
        self.sleeps = []
        self.day = dt.datetime(2026, 11, 17, 12, tzinfo=dt.timezone.utc)

    def monotonic(self):
        return self.t

    def sleep(self, seconds):
        self.sleeps.append(seconds)
        self.t += seconds

    def now(self):
        return self.day


@pytest.fixture(autouse=True)
def quiet(monkeypatch):
    monkeypatch.setattr(scraper, "_polite_sleep", lambda *a, **k: None)


def fetcher(tmp_path, routes, **kwargs):
    clock = Clock()
    session = Session(routes)
    f = PoliteFetcher("shop", min_delay_s=5.0, cache_root=tmp_path, session=session,
                      clock=clock.monotonic, sleep=clock.sleep, now=clock.now, **kwargs)
    return f, session, clock


def test_robots_disallow_is_refused_and_crawl_delay_applies(tmp_path):
    routes = {"https://shop.test/robots.txt": [(200, ROBOTS, {})],
              "https://shop.test/a": [(200, b"A", {})], "https://shop.test/b": [(200, b"B", {})]}
    f, session, clock = fetcher(tmp_path, routes)
    with pytest.raises(DisallowedByRobots):
        f.get("https://shop.test/nobots/x")
    assert f.get("https://shop.test/a").text == "A"
    assert f.get("https://shop.test/b").text == "B"
    assert session.calls == ["https://shop.test/robots.txt", "https://shop.test/a", "https://shop.test/b"]
    # Crawl-delay 12 s for our agent is longer than the 5 s floor.
    assert clock.sleeps and all(s >= 12.0 for s in clock.sleeps)
    assert session.headers["User-Agent"].startswith("horn-simulation-research/")
    assert "github.com/timini/horn-simulation" in session.headers["User-Agent"]
    assert (tmp_path / "shop").is_dir()


def test_missing_robots_allows_and_forbidden_robots_disallows(tmp_path):
    f, _, _ = fetcher(tmp_path, {"https://a.test/x": [(200, b"x", {})]})
    assert f.get("https://a.test/x").ok
    f, _, _ = fetcher(tmp_path / "2", {"https://b.test/robots.txt": [(403, b"", {})]})
    with pytest.raises(DisallowedByRobots):
        f.get("https://b.test/x")


def test_cache_serves_inside_window_and_refetches_after(tmp_path):
    routes = {"https://shop.test/a": [(200, b"one", {}), (200, b"two", {})]}
    f, session, clock = fetcher(tmp_path, routes, cache_days=7)
    first = f.get("https://shop.test/a")
    again = f.get("https://shop.test/a")
    assert (first.from_cache, again.from_cache, again.text) == (False, True, "one")
    assert session.calls.count("https://shop.test/a") == 1
    clock.day += dt.timedelta(days=8)
    fresh = PoliteFetcher("shop", min_delay_s=5.0, cache_root=tmp_path, session=session,
                          clock=clock.monotonic, sleep=clock.sleep, now=clock.now)
    assert fresh.get("https://shop.test/a").text == "two"


def test_404_is_cached_but_429_is_not(tmp_path):
    f, session, _ = fetcher(tmp_path, {"https://shop.test/gone": [(404, b"", {})]})
    assert f.get("https://shop.test/gone").status == 404
    assert f.get("https://shop.test/gone").from_cache


def test_two_consecutive_429_stop_the_source(tmp_path):
    routes = {"https://shop.test/a": [(429, b"", {"Retry-After": "20"})]}
    f, session, clock = fetcher(tmp_path, routes)
    with pytest.raises(SourceBlocked):
        f.get("https://shop.test/a")
    assert 20.0 in clock.sleeps            # Retry-After honoured before the one retry
    assert "429" in f.blocked_reason
    with pytest.raises(SourceBlocked):
        f.get("https://shop.test/other")
    assert session.calls.count("https://shop.test/other") == 0
    assert not (tmp_path / "shop" / f"{url_key('https://shop.test/a')}.body").exists()


def test_one_429_then_success_resets(tmp_path):
    routes = {"https://shop.test/a": [(429, b"", {}), (200, b"ok", {})],
              "https://shop.test/b": [(403, b"", {})], "https://shop.test/c": [(200, b"c", {})]}
    f, _, _ = fetcher(tmp_path, routes)
    assert f.get("https://shop.test/a").text == "ok"
    assert f.get("https://shop.test/b").status == 403
    assert f.get("https://shop.test/c").ok            # resets the run of blocks
    assert f.get("https://shop.test/b").status == 403
    assert f.blocked_reason is None


def test_two_403_in_a_row_stop_the_source(tmp_path):
    routes = {"https://shop.test/a": [(403, b"", {})], "https://shop.test/b": [(403, b"", {})]}
    f, _, _ = fetcher(tmp_path, routes)
    assert f.get("https://shop.test/a").status == 403
    with pytest.raises(SourceBlocked):
        f.get("https://shop.test/b")
    assert f.status()["blocked"]


def test_retry_after_beyond_limit_stops_without_waiting(tmp_path):
    f, _, clock = fetcher(tmp_path, {"https://shop.test/a": [(429, b"", {"Retry-After": "3600"})]})
    with pytest.raises(SourceBlocked, match="Retry-After"):
        f.get("https://shop.test/a")
    assert 3600 not in clock.sleeps


def test_oversized_body_is_refused(tmp_path):
    f, _, _ = fetcher(tmp_path, {"https://shop.test/big": [(200, b"x" * 5000, {})]}, max_bytes=1000)
    with pytest.raises(ResponseTooLarge):
        f.get("https://shop.test/big")


def test_request_returns_final_response_when_asked():
    session = Session({"u": [(404, b"", {})]})
    assert scraper.request(session, "u", delay=0) is None
    assert scraper.request(session, "u", delay=0, return_response=True).status_code == 404


def test_fetching_robots_txt_itself_is_one_request(tmp_path):
    f, session, _ = fetcher(tmp_path, {"https://shop.test/robots.txt": [(200, ROBOTS, {})]})
    assert f.get("https://shop.test/robots.txt").ok
    assert session.calls == ["https://shop.test/robots.txt"]


def test_robots_txt_is_fetchable_even_when_it_disallows_everything(tmp_path):
    f, session, _ = fetcher(tmp_path, {"https://shop.test/robots.txt": [(200, b"User-agent: *\nDisallow: /\n", {})]})
    assert f.get("https://shop.test/robots.txt").ok
    with pytest.raises(DisallowedByRobots):
        f.get("https://shop.test/product")
    assert session.calls == ["https://shop.test/robots.txt"]


def test_html_served_as_robots_means_no_rules(tmp_path):
    routes = {"https://shop.test/robots.txt": [(200, b"<!DOCTYPE html><html>Disallow: /</html>", {"Content-Type": "text/html"})],
              "https://shop.test/p": [(200, b"p", {})]}
    f, _, _ = fetcher(tmp_path, routes)
    assert f.get("https://shop.test/p").ok


def test_cli_robots_snapshot_skips_html(tmp_path):
    from horn_drivers.sources.cli import snapshot_robots
    from horn_drivers.sources.registry import SourceEntry
    entry = SourceEntry(id="shop", name="Shop", kind="manufacturer", base_url="https://shop.test",
                        approval_status="approved", min_delay_s=5, cache_days=90)
    routes = {"https://shop.test/robots.txt": [(200, b"<html>home</html>", {"Content-Type": "text/html"})]}
    f, _, _ = fetcher(tmp_path, routes)
    outcome = snapshot_robots(entry, f, tmp_path / "robots")
    assert outcome["note"].startswith("served an HTML page") and not (tmp_path / "robots" / "shop.txt").exists()
