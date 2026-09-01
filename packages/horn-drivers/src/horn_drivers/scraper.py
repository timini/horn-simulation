"""Bulk scraper for loudspeaker T-S parameters from loudspeakerdatabase.com.

Enumerates manufacturers and drivers by browsing manufacturer pages,
then scrapes individual driver detail pages for Thiele-Small parameters.

The site embeds driver data in JSON ``data-woofer`` attributes, which we
extract directly — more reliable than HTML table parsing.

All values are converted to SI on extraction (mH -> H, cm^2 -> m^2,
mm -> m, g -> kg).  Drivers are classified as "compression" or "cone"
based on Sd, and nominal diameter is inferred geometrically.

Usage:
    horn-scrape-drivers --db data/drivers.json
    horn-scrape-drivers --db data/drivers.json --manufacturers Eminence,BC
    horn-scrape-drivers --db data/drivers.json --dry-run
"""

import argparse
import json
import math
import random
import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

BASE_URL = "https://loudspeakerdatabase.com"
MAX_RETRIES = 3

# Identify ourselves honestly, but send the rest of the headers a normal client
# would: a bare request with no Accept header is a common WAF rejection trigger.
USER_AGENT = (
    "horn-simulation-research/0.2 "
    "(+https://github.com/timini/horn-simulation; academic loudspeaker research)"
)
DEFAULT_HEADERS = {
    "User-Agent": USER_AGENT,
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    "Accept-Language": "en-GB,en;q=0.9",
    "Accept-Encoding": "gzip, deflate, br",
    "Connection": "keep-alive",
    "Upgrade-Insecure-Requests": "1",
}

# Cloudflare returns 52x when it reached the origin but the origin misbehaved.
# These are the site's problem, not ours, and rate limiting will not help: the
# only sane response is to wait a long time and try again.
ORIGIN_ERROR_CODES = frozenset({520, 521, 522, 523, 524, 525, 526, 527})
# Cloudflare/WAF codes that mean we are being asked to slow down or are blocked.
THROTTLE_CODES = frozenset({429, 503})

# Sd thresholds for nominal diameter inference (m^2).
DIAMETER_TABLE: List[Tuple[str, float]] = [
    # (label, nominal Sd in m^2 computed from pi*(d_eff/2)^2)
    ("2in", 0.0020),    # ~50 mm effective
    ("3in", 0.0046),    # ~76 mm effective
    ("4in", 0.0081),    # ~102 mm effective
    ("5in", 0.0110),    # ~118 mm effective
    ("6in", 0.0133),    # ~130 mm effective
    ("8in", 0.0214),    # ~165 mm
    ("10in", 0.0346),   # ~210 mm
    ("12in", 0.0531),   # ~260 mm
    ("15in", 0.0855),   # ~330 mm
    ("18in", 0.1195),   # ~390 mm
]

# Essential parameters — drivers missing any of these are skipped.
ESSENTIAL_PARAMS = ("fs_hz", "re_ohm", "bl_tm", "sd_m2", "mms_kg")


def slugify(manufacturer: str, model: str) -> str:
    """Generate a stable driver_id slug from manufacturer + model."""
    raw = f"{manufacturer}-{model}".lower()
    return re.sub(r"[^a-z0-9]+", "-", raw).strip("-")


def _parse_float(text: str) -> Optional[float]:
    """Extract a float from a string, returning None on failure."""
    if not text:
        return None
    cleaned = re.sub(r"[^\d.\-eE+]", "", text.strip())
    try:
        return float(cleaned)
    except (ValueError, TypeError):
        return None


def build_session():
    """A requests session with the headers a normal client would send."""
    import requests

    session = requests.Session()
    session.headers.update(DEFAULT_HEADERS)
    return session


def _polite_sleep(delay: float, jitter: float = 0.35) -> None:
    """Sleep for *delay* with jitter, so requests are not perfectly periodic."""
    if delay <= 0:
        return
    time.sleep(max(0.0, delay * (1.0 + random.uniform(-jitter, jitter))))


def _retry_after_seconds(resp, fallback: float) -> float:
    """Honour a Retry-After header when the server sends a usable one."""
    raw = resp.headers.get("Retry-After", "").strip()
    if raw.isdigit():
        return float(raw)
    return fallback


def request(
    session,
    url: str,
    *,
    delay: float = 2.0,
    max_attempts: int = 5,
    patience_s: float = 0.0,
    timeout: float = 30.0,
    max_backoff: float = 300.0,
):
    """GET *url* politely, returning a 200 response or None.

    Paces every request with a jittered delay, honours Retry-After on
    throttling responses, and backs off exponentially on transient failures.

    *patience_s* governs how long to keep waiting out a Cloudflare origin
    error (52x). Those mean the site's own backend is failing: retrying fast
    is pointless and rude, but an overnight run should survive an outage
    rather than abandoning the whole scrape. Set it to 0 to fail fast.
    """
    deadline = time.monotonic() + patience_s
    backoff = max(delay, 2.0)
    attempt = 0

    while True:
        attempt += 1
        _polite_sleep(delay)

        status = None
        try:
            resp = session.get(url, timeout=timeout)
            status = resp.status_code
            if status == 200:
                return resp
            reason = f"HTTP {status}"
        except Exception as exc:  # network-level failure
            reason = f"{type(exc).__name__}: {exc}"

        # Permanent: do not retry.
        if status in (400, 401, 403, 404, 410):
            print(f"    {reason} for {url} (not retrying)")
            return None

        waiting_out_origin = status in ORIGIN_ERROR_CODES and time.monotonic() < deadline

        if status in THROTTLE_CODES:
            wait = _retry_after_seconds(resp, backoff)
            print(f"    {reason}: throttled, waiting {wait:.0f}s")
        elif waiting_out_origin:
            wait = min(backoff, max_backoff)
            remaining = (deadline - time.monotonic()) / 60.0
            print(f"    {reason}: origin down, waiting {wait:.0f}s "
                  f"({remaining:.0f} min patience left)")
        elif attempt < max_attempts:
            wait = min(backoff, max_backoff)
            print(f"    {reason}: retry {attempt}/{max_attempts} in {wait:.0f}s")
        else:
            print(f"    {reason}: giving up on {url} after {attempt} attempts")
            return None

        time.sleep(wait)
        backoff = min(backoff * 2, max_backoff)


class OriginWedged(RuntimeError):
    """The site answers 200 but serves the wrong page for the URL requested.

    Seen for real: after an outage the origin came back returning one
    identical driver page for every path, including nonsense ones. A scrape
    that trusted those responses would have overwritten the whole database
    with copies of a single driver, so this is treated as fatal.
    """


def page_identity(soup) -> Optional[str]:
    """The path this page claims to be, from its og:url meta tag."""
    meta = soup.find("meta", attrs={"property": "og:url"})
    if not meta or not meta.get("content"):
        return None
    content = meta["content"].strip()
    if content.startswith(BASE_URL):
        content = content[len(BASE_URL):]
    return "/" + content.lstrip("/")


def _paths_match(requested_url: str, claimed_path: Optional[str]) -> bool:
    """True unless the page positively identifies itself as somewhere else."""
    if claimed_path is None:
        return True  # nothing to check against
    requested = requested_url[len(BASE_URL):] if requested_url.startswith(BASE_URL) else requested_url
    return requested.rstrip("/").lower() == claimed_path.rstrip("/").lower()


def origin_is_up(session, delay: float = 0.0) -> bool:
    """Cheap liveness probe that does not retry."""
    try:
        _polite_sleep(delay)
        return session.get(BASE_URL, timeout=20).status_code == 200
    except Exception:
        return False


def origin_is_healthy(session, delay: float = 0.0) -> bool:
    """Liveness plus routing sanity.

    A 200 is not enough: the origin has been observed serving one identical
    driver page for every path. Ask for a path that cannot exist. A healthy
    site rejects it; a wedged one answers 200 with a page whose og:url names
    somewhere else entirely.

    Comparing response bodies is not good enough, because the pages carry
    something dynamic and two fetches of the same URL differ by a byte.
    """
    from bs4 import BeautifulSoup

    try:
        _polite_sleep(delay)
        if session.get(BASE_URL, timeout=20).status_code != 200:
            return False
        nonsense = f"{BASE_URL}/__horn_sim_probe_{random.randint(10**6, 10**7)}"
        _polite_sleep(delay)
        probe = session.get(nonsense, timeout=20)
    except Exception:
        return False

    if probe.status_code in (404, 410):
        return True  # unknown paths are rejected, so routing works
    if probe.status_code != 200:
        # 429, 403 and friends tell us nothing about routing, and starting a
        # long scrape while being throttled is the wrong move either way.
        print(f"  probe returned HTTP {probe.status_code}; not starting yet")
        return False

    claimed = page_identity(BeautifulSoup(probe.text, "html.parser"))
    if not _paths_match(nonsense, claimed):
        print(f"  origin answers 200 for a nonexistent path and calls it "
              f"{claimed}; routing is broken, treating as down")
        return False
    return True


def wait_for_origin(session, patience_s: float, probe_interval: float = 120.0) -> bool:
    """Block until the site answers, or until *patience_s* runs out."""
    deadline = time.monotonic() + patience_s
    while True:
        if origin_is_healthy(session):
            return True
        if time.monotonic() >= deadline:
            return False
        remaining = (deadline - time.monotonic()) / 60.0
        print(f"  origin still down; re-probing in {probe_interval:.0f}s "
              f"({remaining:.0f} min patience left)", flush=True)
        time.sleep(min(probe_interval, max(1.0, deadline - time.monotonic())))


def infer_driver_type(sd_m2: float) -> str:
    """Classify driver as compression or cone based on Sd."""
    return "compression" if sd_m2 < 0.001 else "cone"


def infer_nominal_diameter(sd_m2: float) -> Optional[str]:
    """Infer nominal diameter from Sd using geometric lookup.

    For compression drivers (Sd < 0.001 m^2), returns None since
    compression driver naming doesn't follow cone diameter conventions.
    """
    if sd_m2 < 0.001:
        return None
    best_label = None
    best_dist = float("inf")
    for label, nominal_sd in DIAMETER_TABLE:
        dist = abs(sd_m2 - nominal_sd)
        if dist < best_dist:
            best_dist = dist
            best_label = label
    return best_label


def _parse_data_woofer(json_str: str) -> Optional[dict]:
    """Parse the data-woofer JSON attribute into SI parameters.

    The site embeds T-S data in data-woofer attributes as JSON like:
    {"re":5.39,"bl":19.39,"mmd":91.2,"rms":2.24,"cms":210.0,
     "sd":1159.0,"le":1.21,"fs":32.0,"qts":0.33,"xmax":8.0,
     "pmax":1600,"frmin":38,"frmax":700,"spl1w":96.8,"z":8}

    Units: sd=cm², le=mH, mmd=grams (dry mass), xmax=mm, z=ohms.
    Note: mmd is dry moving mass; the site separately shows Mms (with
    air load) in the detail page but the JSON only has mmd.
    """
    try:
        raw = json.loads(json_str)
    except (json.JSONDecodeError, TypeError):
        return None

    fs = raw.get("fs")
    re_ohm = raw.get("re")
    bl = raw.get("bl")
    sd_cm2 = raw.get("sd")
    mmd_g = raw.get("mmd")
    le_mh = raw.get("le")

    if not all(v is not None and v > 0 for v in [fs, re_ohm, bl, sd_cm2, mmd_g]):
        return None

    si: Dict[str, float] = {
        "fs_hz": float(fs),
        "re_ohm": float(re_ohm),
        "bl_tm": float(bl),
        "sd_m2": float(sd_cm2) * 1e-4,        # cm² -> m²
        # The JSON carries mmd, the dry moving mass. Keep it: a coupled
        # driver model needs mmd, because the air load is supplied by the
        # horn impedance and must not also be baked into the mass term.
        "mmd_kg": float(mmd_g) * 1e-3,        # g -> kg, dry moving mass
        "mms_kg": float(mmd_g) * 1e-3,        # provisional; HTML Mms overrides
    }

    if le_mh is not None and le_mh > 0:
        si["le_h"] = float(le_mh) * 1e-3      # mH -> H

    if raw.get("qts") is not None:
        si["qts"] = float(raw["qts"])
    if raw.get("z") is not None:
        si["nominal_impedance_ohm"] = float(raw["z"])
    if raw.get("xmax") is not None and raw["xmax"] > 0:
        si["xmax_m"] = float(raw["xmax"]) * 1e-3  # mm -> m
    if raw.get("pmax") is not None and raw["pmax"] > 0:
        pmax = float(raw["pmax"])
        si["peak_power_w"] = pmax             # pmax is program power
        si["power_w"] = pmax / 2.0            # RMS ≈ program / 2

    return si


def scrape_driver_page(
    url: str, session, delay: float = 0.0, patience_s: float = 0.0,
) -> Optional[dict]:
    """Scrape a single driver detail page for T-S parameters.

    Extracts data from the embedded data-woofer JSON attribute (primary),
    falling back to HTML parameter list parsing.  Also extracts Qms and
    Qes from the parameter list (not in the JSON).

    Returns a dict with SI-unit parameters, or None if essential params
    are missing.
    """
    from bs4 import BeautifulSoup

    resp = request(session, url, delay=delay, patience_s=patience_s)
    if resp is None:
        return None

    soup = BeautifulSoup(resp.text, "html.parser")

    # Refuse to trust a page that says it is a different driver: a wedged
    # origin serving one page for every URL would otherwise fill the whole
    # database with copies of it.
    claimed = page_identity(soup)
    if not _paths_match(url, claimed):
        raise OriginWedged(
            f"asked for {url} but the page identifies itself as {claimed}"
        )

    # --- Primary: extract from data-woofer JSON ---
    si = None
    for el in soup.find_all(attrs={"data-woofer": True}):
        candidate = _parse_data_woofer(el["data-woofer"])
        if candidate:
            si = candidate
            break

    if si is None:
        return None

    # --- Supplement: extract Qms, Qes, Mms from the HTML param list ---
    # These aren't in the data-woofer JSON but are in the <li> elements
    highlight_map = {
        "qms": "qms",
        "qes": "qes",
        "mms": "mms_g",  # Mms with air load, in grams
    }

    for key, field in highlight_map.items():
        el = soup.find(attrs={"data-highlight": key})
        if el:
            # Find the <b> tag inside the value span
            val_el = el.find_next(attrs={"data-highlight": key})
            if val_el:
                b_tag = val_el.find("b")
                if b_tag:
                    val = _parse_float(b_tag.get_text(strip=True))
                    if val is not None and val > 0:
                        if field == "mms_g":
                            # Mms includes the measurement air load; mmd_kg
                            # from the JSON is retained alongside it.
                            si["mms_kg"] = val * 1e-3
                        else:
                            si[field] = val

    return si


def discover_manufacturers(
    session, delay: float = 0.0, patience_s: float = 0.0,
) -> List[str]:
    """Discover manufacturer slugs from the homepage.

    The homepage shows driver cards with links like /Eminence/MODEL.
    We extract unique manufacturer slugs from these two-segment paths.
    """
    from bs4 import BeautifulSoup

    resp = request(session, BASE_URL, delay=delay, patience_s=patience_s)
    if resp is None:
        print("ERROR: could not fetch the homepage")
        return []

    soup = BeautifulSoup(resp.text, "html.parser")
    manufacturers = set()

    for link in soup.find_all("a", href=True):
        href = link["href"]
        # Driver pages are like /Manufacturer/Model
        if href.startswith("/") and href.count("/") == 2:
            parts = href.strip("/").split("/")
            if len(parts) == 2:
                mfr = parts[0]
                # Skip non-manufacturer paths
                if mfr.lower() not in (
                    "assets", "images", "simulators", "css", "js",
                ):
                    manufacturers.add(mfr)

    return sorted(manufacturers)


def _extract_driver_links(
    soup, manufacturer_slug: str, seen_urls: set,
) -> List[dict]:
    """Extract driver links from a parsed HTML page/fragment."""
    drivers = []
    for link in soup.find_all("a", href=True):
        href = link["href"]
        if href.startswith("/") and href.count("/") == 2:
            parts = href.strip("/").split("/")
            if len(parts) == 2 and parts[0] == manufacturer_slug:
                full_url = BASE_URL + href
                if full_url not in seen_urls:
                    seen_urls.add(full_url)
                    name = link.get_text(strip=True) or parts[1]
                    drivers.append({
                        "name": name,
                        "url": full_url,
                        "manufacturer": manufacturer_slug,
                    })
    return drivers


def discover_drivers(
    session, manufacturer_slug: str, delay: float = 1.0,
    patience_s: float = 0.0,
) -> List[dict]:
    """Fetch a manufacturer page and return all driver entries.

    Paginates through the site's infinite-scroll API to discover every
    driver, not just the first 40.

    Returns list of dicts with keys: name, url, manufacturer.
    """
    from bs4 import BeautifulSoup

    PAGE_SIZE = 40

    url = f"{BASE_URL}/{manufacturer_slug}"
    resp = request(session, url, delay=delay, patience_s=patience_s)
    if resp is None:
        print(f"  ERROR: could not fetch {url}")
        return []

    soup = BeautifulSoup(resp.text, "html.parser")

    # Extract total result count from the page
    total_count = PAGE_SIZE  # fallback: assume single page
    count_el = soup.find("script", class_="count")
    if count_el:
        try:
            total_count = int(count_el.get_text(strip=True))
        except (ValueError, TypeError):
            pass

    seen_urls: set = set()
    drivers = _extract_driver_links(soup, manufacturer_slug, seen_urls)

    # Paginate if there are more drivers than the initial page
    if total_count > PAGE_SIZE:
        for offset in range(PAGE_SIZE, total_count, PAGE_SIZE):
            page_url = (
                f"{BASE_URL}/next_page_api"
                f"/brand={manufacturer_slug}/offset={offset}"
            )
            resp = request(session, page_url, delay=delay, patience_s=patience_s)
            if resp is None:
                print(f"  WARN: pagination failed at offset={offset}")
                continue
            page_soup = BeautifulSoup(resp.text, "html.parser")
            drivers.extend(
                _extract_driver_links(page_soup, manufacturer_slug, seen_urls)
            )

    return drivers


def _driver_is_current(
    db_dir: Path,
    manufacturer: str,
    driver_id: str,
    required_fields: Tuple[str, ...] = (),
) -> bool:
    """True when this driver is on disk and already has *required_fields*.

    Skipping on existence alone makes a resume cheap, but it also means a
    database written by an older parser never picks up newly captured
    parameters. Naming the fields a current record must have turns the whole
    scrape into an idempotent migration: re-run it and only the stale records
    are re-fetched, and an interrupted run resumes exactly where it stopped.
    """
    path = db_dir / manufacturer / f"{driver_id}.json"
    if not path.exists():
        return False
    if not required_fields:
        return True
    try:
        params = json.loads(path.read_text()).get("parameters", {})
    except (json.JSONDecodeError, OSError):
        return False
    return all(field in params for field in required_fields)


def _load_state(state_path: Optional[Path]) -> dict:
    """Per-manufacturer progress from a previous run, for reporting and resume."""
    if state_path is None or not state_path.exists():
        return {}
    try:
        return json.loads(state_path.read_text())
    except (json.JSONDecodeError, OSError):
        return {}


def _save_state(state_path: Optional[Path], state: dict) -> None:
    if state_path is None:
        return
    state_path.parent.mkdir(parents=True, exist_ok=True)
    state_path.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n")


def _save_driver(db_dir: Path, driver: dict) -> None:
    """Write a single driver dict to its JSON file."""
    manufacturer = driver.get("manufacturer", "unknown")
    driver_id = driver["driver_id"]
    mfr_dir = db_dir / manufacturer
    mfr_dir.mkdir(parents=True, exist_ok=True)
    driver_file = mfr_dir / f"{driver_id}.json"
    driver_file.write_text(json.dumps(driver, indent=4) + "\n")


def scrape_all(
    db_dir: Optional[Path] = None,
    max_manufacturers: Optional[int] = None,
    manufacturer_filter: Optional[List[str]] = None,
    delay: float = 2.0,
    patience_s: float = 0.0,
    refresh: bool = False,
    state_path: Optional[Path] = None,
    required_fields: Tuple[str, ...] = (),
) -> int:
    """Scrape drivers from loudspeakerdatabase.com.

    When *db_dir* is provided, each driver is written to disk immediately
    after scraping (one JSON file per driver).  When *db_dir* is ``None``
    (dry-run mode), drivers are printed but not saved.

    Writing per driver makes the run resumable: by default a driver whose
    JSON already exists is skipped, so an interrupted scrape can simply be
    re-run. Pass *refresh* to re-fetch everything, which is what you want
    after changing how parameters are parsed or classified.

    *patience_s* is how long to keep waiting out a site outage before giving
    up, so an overnight run survives the origin going down.

    Returns the number of drivers successfully scraped.
    """
    session = build_session()
    state = _load_state(state_path)

    # Discover manufacturers or use provided list
    if manufacturer_filter:
        manufacturer_slugs = manufacturer_filter
        print(f"Using {len(manufacturer_slugs)} specified manufacturers: "
              f"{', '.join(manufacturer_slugs)}")
    else:
        print("Discovering manufacturers from homepage...")
        manufacturer_slugs = discover_manufacturers(
            session, delay=delay, patience_s=patience_s,
        )
        print(f"Found {len(manufacturer_slugs)} manufacturers: "
              f"{', '.join(manufacturer_slugs)}")

    if max_manufacturers:
        manufacturer_slugs = manufacturer_slugs[:max_manufacturers]

    total_scraped = 0

    for i, mfr_slug in enumerate(manufacturer_slugs):
        print(f"\n[{i + 1}/{len(manufacturer_slugs)}] {mfr_slug}")

        driver_entries = discover_drivers(
            session, mfr_slug, delay=delay, patience_s=patience_s,
        )
        print(f"  Found {len(driver_entries)} drivers")
        if not driver_entries:
            print(f"  {mfr_slug}: no drivers discovered, leaving for a later run")
            continue

        mfr_scraped = 0
        mfr_skipped = 0
        for entry in driver_entries:
            driver_id = slugify(entry["manufacturer"], entry["name"])
            if (not refresh) and db_dir is not None and _driver_is_current(
                db_dir, entry["manufacturer"], driver_id, required_fields
            ):
                mfr_skipped += 1
                continue

            print(f"  Scraping: {entry['name']}", end="", flush=True)
            try:
                params = scrape_driver_page(
                    entry["url"], session, delay=delay, patience_s=patience_s,
                )
            except OriginWedged:
                raise
            except Exception as e:
                print(f"\n    ERROR: {e}")
                params = None

            if params is None:
                print(" -> SKIP (no params)")
                continue

            # Check essential params
            missing = [p for p in ESSENTIAL_PARAMS if p not in params]
            if missing:
                print(f" -> SKIP (missing {', '.join(missing)})")
                continue

            sd = params["sd_m2"]
            driver_type = infer_driver_type(sd)
            nominal_diameter = infer_nominal_diameter(sd)

            model = entry["name"]
            manufacturer = entry["manufacturer"]

            driver = {
                "driver_id": driver_id,
                "manufacturer": manufacturer,
                "model_name": model,
                "driver_type": driver_type,
                "parameters": params,
            }
            if nominal_diameter:
                driver["nominal_diameter"] = nominal_diameter

            # Write immediately
            if db_dir is not None:
                _save_driver(db_dir, driver)

            mfr_scraped += 1
            total_scraped += 1
            tag = f"{driver_type}, {nominal_diameter}" if nominal_diameter else driver_type
            print(f" -> OK ({tag})")

        print(f"  {mfr_slug}: {mfr_scraped} new, {mfr_skipped} already present, "
              f"{len(driver_entries)} listed")
        state[mfr_slug] = {
            "listed": len(driver_entries),
            "scraped": mfr_scraped,
            "skipped": mfr_skipped,
            "completed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        }
        _save_state(state_path, state)

    return total_scraped


def main():
    parser = argparse.ArgumentParser(
        description="Scrape loudspeaker T-S parameters from loudspeakerdatabase.com"
    )
    parser.add_argument(
        "--db", type=str, default="data/drivers",
        help="Path to driver database directory (default: data/drivers)",
    )
    parser.add_argument(
        "--delay", type=float, default=2.0,
        help="Base delay between requests in seconds, jittered (default: 2.0)",
    )
    parser.add_argument(
        "--patience-hours", type=float, default=0.0,
        help="How long to keep waiting out a site outage (Cloudflare 52x) "
             "before giving up. 0 fails fast; use e.g. 10 for an overnight run.",
    )
    parser.add_argument(
        "--require-fields", type=str, default="mmd_kg",
        help="Comma-separated parameter names a stored driver must already "
             "have to be considered current. Records missing any of them are "
             "re-fetched, which makes a parser change an idempotent migration. "
             "Pass an empty string to skip purely on file existence.",
    )
    parser.add_argument(
        "--refresh", action="store_true",
        help="Re-fetch drivers that already exist in the database. Without "
             "this, existing files are skipped so a run resumes cheaply.",
    )
    parser.add_argument(
        "--state", type=str, default=None,
        help="Path to a JSON progress file recording per-manufacturer counts.",
    )
    parser.add_argument(
        "--wait-for-origin", action="store_true",
        help="Probe the site first and block until it responds, within the "
             "patience budget, instead of failing immediately.",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print scraped data without writing to DB",
    )
    parser.add_argument(
        "--max-manufacturers", type=int, default=None,
        help="Limit number of manufacturers to scrape",
    )
    parser.add_argument(
        "--manufacturers", type=str, default=None,
        help="Comma-separated manufacturer slugs as on the site (e.g. Eminence,BC,Oberton)",
    )
    args = parser.parse_args()

    mfr_filter = None
    if args.manufacturers:
        mfr_filter = [m.strip() for m in args.manufacturers.split(",")]

    db_dir = None if args.dry_run else Path(args.db)
    patience_s = max(0.0, args.patience_hours) * 3600.0

    if args.wait_for_origin:
        print(f"Probing {BASE_URL} (patience {args.patience_hours:.1f} h)...")
        if not wait_for_origin(build_session(), patience_s):
            print("ERROR: the site never came back within the patience budget.")
            return 1
        print("Site is responding; starting.")

    try:
        count = scrape_all(
            db_dir=db_dir,
            max_manufacturers=args.max_manufacturers,
            manufacturer_filter=mfr_filter,
            delay=args.delay,
            patience_s=patience_s,
            refresh=args.refresh,
            state_path=Path(args.state) if args.state else None,
            required_fields=tuple(
                f.strip() for f in args.require_fields.split(",") if f.strip()
            ),
        )
    except OriginWedged as exc:
        print(f"\nABORTED: {exc}")
        print("The site is answering but serving the wrong pages. Nothing was "
              "written from this point on; re-run when it recovers.")
        return 2
    print(f"\nScraped {count} drivers total.")
    if count == 0:
        print("Nothing was scraped. Treating this as a failed run so a retry "
              "does not mistake it for a completed database.")
    if db_dir:
        total = sum(1 for _ in db_dir.rglob("*.json"))
        print(f"Database: {db_dir} ({total} drivers)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
