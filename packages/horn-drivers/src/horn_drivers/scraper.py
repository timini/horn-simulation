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
import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

BASE_URL = "https://loudspeakerdatabase.com"
MAX_RETRIES = 3

# Sd thresholds for nominal diameter inference (m^2).
DIAMETER_TABLE: List[Tuple[str, float]] = [
    # (label, nominal Sd in m^2 computed from pi*(d_eff/2)^2)
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


def infer_driver_type(sd_m2: float) -> str:
    """Classify driver as compression or cone based on Sd."""
    return "compression" if sd_m2 < 0.005 else "cone"


def infer_nominal_diameter(sd_m2: float) -> Optional[str]:
    """Infer nominal diameter from Sd using geometric lookup.

    For compression drivers (Sd < 0.005 m^2), returns None since
    compression driver naming doesn't follow cone diameter conventions.
    """
    if sd_m2 < 0.005:
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
        "mms_kg": float(mmd_g) * 1e-3,        # g -> kg (mmd ≈ mms for scraping)
    }

    if le_mh is not None and le_mh > 0:
        si["le_h"] = float(le_mh) * 1e-3      # mH -> H

    if raw.get("qts") is not None:
        si["qts"] = float(raw["qts"])
    if raw.get("z") is not None:
        si["nominal_impedance_ohm"] = float(raw["z"])
    if raw.get("xmax") is not None and raw["xmax"] > 0:
        si["xmax_m"] = float(raw["xmax"]) * 1e-3  # mm -> m

    return si


def scrape_driver_page(url: str, session) -> Optional[dict]:
    """Scrape a single driver detail page for T-S parameters.

    Extracts data from the embedded data-woofer JSON attribute (primary),
    falling back to HTML parameter list parsing.  Also extracts Qms and
    Qes from the parameter list (not in the JSON).

    Returns a dict with SI-unit parameters, or None if essential params
    are missing.
    """
    from bs4 import BeautifulSoup

    resp = session.get(url, timeout=15)
    if resp.status_code != 200:
        print(f"  WARN: HTTP {resp.status_code} for {url}")
        return None

    soup = BeautifulSoup(resp.text, "html.parser")

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
                            si["mms_kg"] = val * 1e-3  # Override mmd with Mms
                        else:
                            si[field] = val

    return si


def discover_manufacturers(session) -> List[str]:
    """Discover manufacturer slugs from the homepage.

    The homepage shows driver cards with links like /Eminence/MODEL.
    We extract unique manufacturer slugs from these two-segment paths.
    """
    from bs4 import BeautifulSoup

    resp = session.get(BASE_URL, timeout=15)
    if resp.status_code != 200:
        print(f"ERROR: HTTP {resp.status_code} fetching homepage")
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
) -> List[dict]:
    """Fetch a manufacturer page and return all driver entries.

    Paginates through the site's infinite-scroll API to discover every
    driver, not just the first 40.

    Returns list of dicts with keys: name, url, manufacturer.
    """
    from bs4 import BeautifulSoup

    PAGE_SIZE = 40

    url = f"{BASE_URL}/{manufacturer_slug}"
    resp = None
    for attempt in range(MAX_RETRIES):
        try:
            resp = session.get(url, timeout=15)
            break
        except Exception as e:
            if attempt < MAX_RETRIES - 1:
                wait = delay * (attempt + 2)
                print(f"  RETRY manufacturer page ({attempt + 1}/{MAX_RETRIES}): {e}")
                time.sleep(wait)
            else:
                print(f"  ERROR: failed to fetch {url} after {MAX_RETRIES} attempts")
                return []

    if resp is None or resp.status_code != 200:
        print(f"  WARN: HTTP {resp.status_code if resp else '?'} for {url}")
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
            time.sleep(delay)
            page_url = (
                f"{BASE_URL}/next_page_api"
                f"/brand={manufacturer_slug}/offset={offset}"
            )
            for attempt in range(MAX_RETRIES):
                try:
                    resp = session.get(page_url, timeout=15)
                    if resp.status_code != 200:
                        print(f"  WARN: HTTP {resp.status_code} for page offset={offset}")
                        break
                    page_soup = BeautifulSoup(resp.text, "html.parser")
                    new = _extract_driver_links(
                        page_soup, manufacturer_slug, seen_urls,
                    )
                    drivers.extend(new)
                    break
                except Exception as e:
                    if attempt < MAX_RETRIES - 1:
                        wait = delay * (attempt + 2)
                        print(f"  RETRY page offset={offset} "
                              f"({attempt + 1}/{MAX_RETRIES}): {e}")
                        time.sleep(wait)
                    else:
                        print(f"  WARN: pagination failed at offset={offset} "
                              f"after {MAX_RETRIES} attempts: {e}")

    return drivers


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
    delay: float = 1.0,
) -> int:
    """Scrape drivers from loudspeakerdatabase.com.

    When *db_dir* is provided, each driver is written to disk immediately
    after scraping (one JSON file per driver).  When *db_dir* is ``None``
    (dry-run mode), drivers are printed but not saved.

    Returns the number of drivers successfully scraped.
    """
    import requests

    session = requests.Session()
    session.headers["User-Agent"] = (
        "horn-simulation-research/0.1 "
        "(+https://github.com/timini/horn; academic loudspeaker research)"
    )

    # Discover manufacturers or use provided list
    if manufacturer_filter:
        manufacturer_slugs = manufacturer_filter
        print(f"Using {len(manufacturer_slugs)} specified manufacturers: "
              f"{', '.join(manufacturer_slugs)}")
    else:
        print("Discovering manufacturers from homepage...")
        time.sleep(delay)
        manufacturer_slugs = discover_manufacturers(session)
        print(f"Found {len(manufacturer_slugs)} manufacturers: "
              f"{', '.join(manufacturer_slugs)}")

    if max_manufacturers:
        manufacturer_slugs = manufacturer_slugs[:max_manufacturers]

    total_scraped = 0

    for i, mfr_slug in enumerate(manufacturer_slugs):
        print(f"\n[{i + 1}/{len(manufacturer_slugs)}] {mfr_slug}")

        time.sleep(delay)
        driver_entries = discover_drivers(session, mfr_slug, delay=delay)
        print(f"  Found {len(driver_entries)} drivers")

        mfr_scraped = 0
        for entry in driver_entries:
            time.sleep(delay)
            print(f"  Scraping: {entry['name']}", end="", flush=True)

            params = None
            for attempt in range(MAX_RETRIES):
                try:
                    params = scrape_driver_page(entry["url"], session)
                    break
                except Exception as e:
                    if attempt < MAX_RETRIES - 1:
                        wait = delay * (attempt + 2)
                        print(f"\n    RETRY ({attempt + 1}/{MAX_RETRIES}): {e}, "
                              f"waiting {wait:.0f}s", end="", flush=True)
                        time.sleep(wait)
                    else:
                        print(f"\n    ERROR after {MAX_RETRIES} attempts: {e}")
                        break

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
            driver_id = slugify(manufacturer, model)

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

        print(f"  {mfr_slug}: {mfr_scraped}/{len(driver_entries)} scraped")

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
        "--delay", type=float, default=1.0,
        help="Delay between requests in seconds (default: 1.0)",
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
    count = scrape_all(
        db_dir=db_dir,
        max_manufacturers=args.max_manufacturers,
        manufacturer_filter=mfr_filter,
        delay=args.delay,
    )
    print(f"\nScraped {count} drivers total.")
    if db_dir:
        total = sum(1 for _ in db_dir.rglob("*.json"))
        print(f"Database: {db_dir} ({total} drivers)")


if __name__ == "__main__":
    main()
