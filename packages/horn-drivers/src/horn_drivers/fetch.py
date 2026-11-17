"""Shared polite HTTP layer for every driver and price source.

One ``PoliteFetcher`` per source per run. It:

- reads robots.txt once per host per run (stored in the cache) and refuses
  disallowed URLs with ``DisallowedByRobots``;
- sends one request at a time, at least ``min_delay_s`` (or the robots
  Crawl-delay for our agent, whichever is longer) plus jitter apart per host;
- caches every response it keeps (200, 404, 410) under
  ``<cache_root>/<source_id>/`` keyed by a hash of the URL, with the fetch
  time, and serves the cached copy while it is younger than ``cache_days``;
- stops the source for the run after two consecutive 429 or 403 responses
  (``SourceBlocked``), and records why;
- identifies itself with the project User-Agent (project name and public
  repository URL) through ``horn_drivers.scraper.request``.

It never retries around a block, rotates identity or imitates a browser.
"""

import datetime as _dt
import hashlib
import json
import random
import time
import urllib.robotparser
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Optional
from urllib.parse import urlsplit

from horn_drivers.scraper import USER_AGENT, ResponseTooLarge, _retry_after_seconds, build_session, request

ROBOTS_AGENT = "horn-simulation-research"
DEFAULT_CACHE_ROOT = Path("results/cache")
BLOCK_STATUSES = frozenset({403, 429})
BLOCK_AFTER = 2
CACHEABLE_STATUSES = frozenset({200, 404, 410})
MAX_RETRY_AFTER_S = 300.0
DEFAULT_MAX_BYTES = 20 * 1024 * 1024


class FetchError(RuntimeError):
    """Base class for fetch refusals."""


class DisallowedByRobots(FetchError):
    """robots.txt disallows this URL for our agent."""


class SourceBlocked(FetchError):
    """The source answered 429 or 403 twice in a row; it is stopped for this run."""


@dataclass
class FetchResult:
    url: str
    status: Optional[int]
    content: bytes
    content_type: str
    fetched_at: str  # ISO 8601 UTC
    from_cache: bool
    cache_path: Optional[Path] = None

    @property
    def ok(self) -> bool:
        return self.status == 200

    @property
    def text(self) -> str:
        return self.content.decode("utf-8", errors="replace")

    @property
    def retrieved_on(self) -> str:
        return self.fetched_at[:10]


def url_key(url: str) -> str:
    return hashlib.sha256(url.encode("utf-8")).hexdigest()[:40]


@dataclass
class PoliteFetcher:
    source_id: str
    min_delay_s: float = 5.0
    cache_days: float = 7.0
    cache_root: Path = DEFAULT_CACHE_ROOT
    jitter: float = 0.35
    timeout_s: float = 30.0
    max_bytes: int = DEFAULT_MAX_BYTES
    session: object = None
    clock: Callable[[], float] = time.monotonic
    sleep: Callable[[float], None] = time.sleep
    now: Callable[[], _dt.datetime] = lambda: _dt.datetime.now(_dt.timezone.utc)
    blocked_reason: Optional[str] = None
    requests_made: int = 0
    _consecutive_blocks: int = 0
    _last_request: Dict[str, float] = field(default_factory=dict)
    _robots: Dict[str, urllib.robotparser.RobotFileParser] = field(default_factory=dict)

    def __post_init__(self):
        if self.min_delay_s < 0 or self.cache_days < 0:
            raise ValueError("delay and cache window must be nonnegative")
        self.cache_root = Path(self.cache_root)
        if self.session is None:
            self.session = build_session()
        self.session.headers["User-Agent"] = USER_AGENT

    # -- cache -------------------------------------------------------------

    @property
    def cache_dir(self) -> Path:
        return self.cache_root / self.source_id

    def _paths(self, url: str):
        key = url_key(url)
        return self.cache_dir / f"{key}.json", self.cache_dir / f"{key}.body"

    def cached(self, url: str, max_age_days: Optional[float] = None) -> Optional[FetchResult]:
        """The cached response for *url* if it is inside the cache window."""
        meta_path, body_path = self._paths(url)
        if not meta_path.exists() or not body_path.exists():
            return None
        meta = json.loads(meta_path.read_text())
        window = self.cache_days if max_age_days is None else max_age_days
        age = self.now() - _dt.datetime.fromisoformat(meta["fetched_at"])
        if age.total_seconds() > window * 86400:
            return None
        return FetchResult(url=url, status=meta["status"], content=body_path.read_bytes(),
                           content_type=meta.get("content_type", ""), fetched_at=meta["fetched_at"],
                           from_cache=True, cache_path=body_path)

    def _store(self, result: FetchResult) -> None:
        meta_path, body_path = self._paths(result.url)
        meta_path.parent.mkdir(parents=True, exist_ok=True)
        body_path.write_bytes(result.content)
        meta_path.write_text(json.dumps({
            "url": result.url, "status": result.status, "content_type": result.content_type,
            "fetched_at": result.fetched_at, "bytes": len(result.content),
            "sha256": hashlib.sha256(result.content).hexdigest(), "source_id": self.source_id,
        }, indent=2))
        result.cache_path = body_path

    # -- pacing and robots ---------------------------------------------------

    def _wait_turn(self, host: str, delay: float) -> None:
        last = self._last_request.get(host)
        if last is not None:
            target = delay * (1.0 + random.uniform(0.0, self.jitter))
            remaining = target - (self.clock() - last)
            if remaining > 0:
                self.sleep(remaining)
        self._last_request[host] = self.clock()

    def _raw_get(self, url: str, delay: float):
        """One paced GET through scraper.request; returns the response or None."""
        self._wait_turn(urlsplit(url).netloc, delay)
        self.requests_made += 1
        return request(self.session, url, delay=0.0, max_attempts=1, timeout=self.timeout_s,
                       return_response=True, max_bytes=self.max_bytes)

    def robots(self, url: str) -> urllib.robotparser.RobotFileParser:
        """Parsed robots.txt for the host of *url*, fetched once per run and stored."""
        parts = urlsplit(url)
        host = parts.netloc
        if host in self._robots:
            return self._robots[host]
        robots_url = f"{parts.scheme}://{host}/robots.txt"
        parser = urllib.robotparser.RobotFileParser(robots_url)
        cached = self.cached(robots_url, max_age_days=1.0)
        if cached is None:
            resp = self._raw_get(robots_url, self.min_delay_s)
            status = getattr(resp, "status_code", None)
            content = resp.content if resp is not None else b""
            cached = FetchResult(robots_url, status, content, _content_type(resp), self._now_iso(), False)
            if status is not None:
                self._store(cached)
        if cached.status == 200 and _is_html(cached):
            parser.allow_all = True          # home page served for /robots.txt: no rules
        elif cached.status == 200:
            parser.parse(cached.text.splitlines())
        elif cached.status in (401, 403):
            parser.disallow_all = True       # access to robots.txt itself refused
        elif cached.status is not None and 400 <= cached.status < 500:
            parser.allow_all = True          # no robots.txt: everything allowed
        else:
            parser.disallow_all = True       # unreachable or 5xx: assume nothing allowed
        self._robots[host] = parser
        return parser

    def allowed(self, url: str) -> bool:
        return self.robots(url).can_fetch(ROBOTS_AGENT, url)

    def crawl_delay(self, url: str) -> float:
        robots = self.robots(url)
        try:
            value = robots.crawl_delay(ROBOTS_AGENT)
        except AttributeError:
            value = None
        return max(self.min_delay_s, float(value or 0.0))

    # -- the one public fetch -----------------------------------------------

    def get(self, url: str, *, max_age_days: Optional[float] = None) -> FetchResult:
        """Fetch *url* politely, or return the cached copy inside the window.

        Raises SourceBlocked once the source is stopped, DisallowedByRobots
        for a disallowed URL and ResponseTooLarge for an oversized body.
        A non-200 result is returned with its status (404s are cached too).
        """
        if self.blocked_reason:
            raise SourceBlocked(f"{self.source_id} stopped for this run: {self.blocked_reason}")
        if not url.startswith(("http://", "https://")):
            raise ValueError(f"not an http(s) URL: {url}")
        cached = self.cached(url, max_age_days)
        if cached is not None:
            return cached
        if urlsplit(url).path != "/robots.txt" and not self.allowed(url):
            raise DisallowedByRobots(f"robots.txt disallows {url} for {ROBOTS_AGENT}")
        self.robots(url)
        cached = self.cached(url, max_age_days)   # the robots check may have just fetched this URL
        if cached is not None:
            return cached
        delay = self.crawl_delay(url)

        resp = self._raw_get(url, delay)
        status = getattr(resp, "status_code", None)
        if status in BLOCK_STATUSES:
            self._count_block(url, status)
            if status == 429:
                # One polite retry after Retry-After; a second 429 stops the source.
                wait = _retry_after_seconds(resp, 2.0 * delay)
                if wait > MAX_RETRY_AFTER_S:
                    self._block(f"HTTP 429 with Retry-After {wait:.0f} s at {url}")
                self.sleep(max(wait, delay))
                resp = self._raw_get(url, delay)
                status = getattr(resp, "status_code", None)
                if status in BLOCK_STATUSES:
                    self._count_block(url, status)
        if status is not None and status not in BLOCK_STATUSES:
            self._consecutive_blocks = 0

        result = FetchResult(url=url, status=status, content=(resp.content if resp is not None else b""),
                             content_type=_content_type(resp), fetched_at=self._now_iso(), from_cache=False)
        if status in CACHEABLE_STATUSES:
            self._store(result)
        return result

    def _count_block(self, url: str, status: int) -> None:
        """Count a 429/403 and stop the source on the second in a row."""
        self._consecutive_blocks += 1
        if self._consecutive_blocks >= BLOCK_AFTER:
            self._block(f"{BLOCK_AFTER} consecutive HTTP 429/403 responses, last HTTP {status} at {url}")

    def _block(self, reason: str) -> None:
        self.blocked_reason = reason
        raise SourceBlocked(f"{self.source_id} stopped for this run: {reason}")

    def _now_iso(self) -> str:
        return self.now().astimezone(_dt.timezone.utc).isoformat(timespec="seconds")

    def status(self) -> dict:
        """Run summary for logs and the source report."""
        return {"source_id": self.source_id, "requests_made": self.requests_made,
                "blocked": self.blocked_reason is not None, "blocked_reason": self.blocked_reason}


def _is_html(result: FetchResult) -> bool:
    head = result.content[:512].lstrip().lower()
    return "html" in result.content_type.lower() or head.startswith((b"<!doctype", b"<html"))


def _content_type(resp) -> str:
    if resp is None:
        return ""
    return (getattr(resp, "headers", {}) or {}).get("Content-Type", "") or ""


__all__ = ["PoliteFetcher", "FetchResult", "FetchError", "DisallowedByRobots", "SourceBlocked",
           "ResponseTooLarge", "ROBOTS_AGENT", "url_key"]
