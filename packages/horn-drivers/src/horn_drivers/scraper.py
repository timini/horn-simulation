"""Bulk scraper for loudspeaker T-S parameters from loudspeakerdatabase.com.

Enumerates manufacturers and drivers by browsing the site's index pages,
then scrapes individual driver detail pages for Thiele-Small parameters.

All values are converted to SI on extraction (mH -> H, cm^2 -> m^2,
mm -> m, g -> kg).  Drivers are classified as "compression" or "cone"
based on Sd, and nominal diameter is inferred geometrically.

Usage:
    horn-scrape-drivers --db data/drivers.json
    horn-scrape-drivers --db data/drivers.json --manufacturers Eminence,B&C
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

# Sd thresholds for nominal diameter inference (m^2).
# Each entry: (label, min_sd, max_sd) — ranges overlap slightly; we pick
# the closest geometric match.
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
    # Find closest match in diameter table
    best_label = None
    best_dist = float("inf")
    for label, nominal_sd in DIAMETER_TABLE:
        dist = abs(sd_m2 - nominal_sd)
        if dist < best_dist:
            best_dist = dist
            best_label = label
    return best_label


def scrape_driver_page(url: str, session) -> Optional[dict]:
    """Scrape a single driver detail page for T-S parameters.

    Returns a dict with SI-unit parameters, or None if essential params
    are missing.
    """
    from bs4 import BeautifulSoup

    resp = session.get(url, timeout=15)
    if resp.status_code != 200:
        print(f"  WARN: HTTP {resp.status_code} for {url}")
        return None

    soup = BeautifulSoup(resp.text, "html.parser")

    params: Dict[str, Optional[float]] = {}
    # Map page labels to intermediate field names.
    # Values will be converted to SI after collection.
    param_map = {
        "fs": "fs_hz",
        "re": "re_ohm",
        "bl": "bl_tm",
        "sd": "sd_cm2",
        "mms": "mms_kg",
        "mmd": "mmd_g",  # sometimes listed instead of Mms
        "le": "le_mh",
        "qms": "qms",
        "qes": "qes",
        "qts": "qts",
        "xmax": "xmax_mm",
        "vas": "vas_l",
        "nominal impedance": "z_nom",
        "impedance": "z_nom",
    }

    for row in soup.find_all("tr"):
        cells = row.find_all(["td", "th"])
        if len(cells) < 2:
            continue
        label_text = cells[0].get_text(strip=True).lower()
        value_text = cells[1].get_text(strip=True)

        for key, field_name in param_map.items():
            if key in label_text:
                val = _parse_float(value_text)
                if val is not None:
                    params[field_name] = val
                break

    if not params.get("fs_hz"):
        return None

    # Convert to SI
    si: Dict[str, Optional[float]] = {
        "fs_hz": params.get("fs_hz"),
        "re_ohm": params.get("re_ohm"),
        "bl_tm": params.get("bl_tm"),
        "mms_kg": params.get("mms_kg"),
        "qms": params.get("qms"),
        "qes": params.get("qes"),
        "qts": params.get("qts"),
    }

    # Sd: cm^2 -> m^2
    if params.get("sd_cm2"):
        si["sd_m2"] = params["sd_cm2"] * 1e-4

    # Le: mH -> H
    if params.get("le_mh"):
        si["le_h"] = params["le_mh"] * 1e-3

    # Xmax: mm -> m
    if params.get("xmax_mm"):
        si["xmax_m"] = params["xmax_mm"] * 1e-3

    # Mmd (dry mass in grams) -> Mms approximation if Mms not given
    if not si.get("mms_kg") and params.get("mmd_g"):
        # Mms ≈ Mmd + air load; for scraping purposes use Mmd directly
        # (air load is small and typically included in published Mms)
        si["mms_kg"] = params["mmd_g"] * 1e-3

    # Nominal impedance
    if params.get("z_nom"):
        si["nominal_impedance_ohm"] = params["z_nom"]

    return {k: v for k, v in si.items() if v is not None}


def discover_manufacturers(session) -> List[str]:
    """Fetch the site index and return a list of manufacturer page URLs.

    The main page at loudspeakerdatabase.com lists manufacturers as links.
    """
    from bs4 import BeautifulSoup

    resp = session.get(BASE_URL, timeout=15)
    if resp.status_code != 200:
        print(f"ERROR: HTTP {resp.status_code} fetching manufacturer index")
        return []

    soup = BeautifulSoup(resp.text, "html.parser")
    manufacturers = []

    for link in soup.find_all("a", href=True):
        href = link["href"]
        # Manufacturer pages are top-level paths like /Eminence, /B&C
        # Skip known non-manufacturer paths
        if href.startswith("/") and "/" not in href[1:]:
            name = href[1:]
            if name and name.lower() not in (
                "search", "about", "contact", "login", "register",
                "privacy", "terms", "faq", "api", "sitemap",
                "favicon.ico", "robots.txt",
            ):
                full_url = BASE_URL + href
                if full_url not in manufacturers:
                    manufacturers.append(full_url)

    return manufacturers


def discover_drivers(session, manufacturer_url: str) -> List[dict]:
    """Fetch a manufacturer page and return driver entries.

    Returns list of dicts with keys: name, url, manufacturer.
    """
    from bs4 import BeautifulSoup

    resp = session.get(manufacturer_url, timeout=15)
    if resp.status_code != 200:
        print(f"  WARN: HTTP {resp.status_code} for {manufacturer_url}")
        return []

    soup = BeautifulSoup(resp.text, "html.parser")
    manufacturer = manufacturer_url.rstrip("/").rsplit("/", 1)[-1]
    drivers = []

    for link in soup.find_all("a", href=True):
        href = link["href"]
        # Driver pages are like /Eminence/KAPPA_PRO-15A
        if href.startswith("/") and href.count("/") == 2:
            parts = href.strip("/").split("/")
            if len(parts) == 2 and parts[0].lower() == manufacturer.lower():
                name = link.get_text(strip=True) or parts[1]
                full_url = BASE_URL + href
                if not any(d["url"] == full_url for d in drivers):
                    drivers.append({
                        "name": name,
                        "url": full_url,
                        "manufacturer": manufacturer,
                    })

    return drivers


def scrape_all(
    max_manufacturers: Optional[int] = None,
    manufacturer_filter: Optional[List[str]] = None,
    delay: float = 1.0,
) -> List[dict]:
    """Scrape all drivers from loudspeakerdatabase.com.

    Args:
        max_manufacturers: Limit number of manufacturers to scrape.
        manufacturer_filter: Only scrape these manufacturers (case-insensitive).
        delay: Seconds between HTTP requests (rate limiting).

    Returns:
        List of driver dicts ready for database insertion.
    """
    import requests

    session = requests.Session()
    session.headers["User-Agent"] = (
        "horn-simulation-research/0.1 "
        "(+https://github.com/timini/horn; academic loudspeaker research)"
    )

    print("Discovering manufacturers...")
    time.sleep(delay)
    manufacturer_urls = discover_manufacturers(session)
    print(f"Found {len(manufacturer_urls)} manufacturers")

    # Apply filters
    if manufacturer_filter:
        filter_lower = [m.lower() for m in manufacturer_filter]
        manufacturer_urls = [
            url for url in manufacturer_urls
            if url.rstrip("/").rsplit("/", 1)[-1].lower() in filter_lower
        ]
        print(f"Filtered to {len(manufacturer_urls)} manufacturers: "
              f"{', '.join(url.rsplit('/', 1)[-1] for url in manufacturer_urls)}")

    if max_manufacturers:
        manufacturer_urls = manufacturer_urls[:max_manufacturers]

    drivers: List[dict] = []

    for i, mfr_url in enumerate(manufacturer_urls):
        mfr_name = mfr_url.rstrip("/").rsplit("/", 1)[-1]
        print(f"\n[{i + 1}/{len(manufacturer_urls)}] {mfr_name}")

        time.sleep(delay)
        driver_entries = discover_drivers(session, mfr_url)
        print(f"  Found {len(driver_entries)} drivers")

        for entry in driver_entries:
            time.sleep(delay)
            print(f"  Scraping: {entry['name']}")

            try:
                params = scrape_driver_page(entry["url"], session)
            except Exception as e:
                print(f"    ERROR: {e}")
                continue

            if params is None:
                print("    SKIP: no T-S parameters found")
                continue

            # Check essential params
            missing = [p for p in ESSENTIAL_PARAMS if p not in params]
            if missing:
                print(f"    SKIP: missing {', '.join(missing)}")
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

            drivers.append(driver)
            print(f"    OK: {driver_id} ({driver_type}"
                  f"{', ' + nominal_diameter if nominal_diameter else ''})")

    return drivers


def merge_into_db(db_path: Path, new_drivers: List[dict]) -> dict:
    """Merge scraped drivers into existing database (idempotent, keyed by driver_id)."""
    if db_path.exists():
        db = json.loads(db_path.read_text())
    else:
        db = {"schema_version": 2, "drivers": []}

    if "drivers" not in db:
        drivers_list = []
        for did, data in db.items():
            data.setdefault("driver_id", did)
            drivers_list.append(data)
        db = {"schema_version": 2, "drivers": drivers_list}

    existing_ids = {d["driver_id"] for d in db["drivers"]}
    added, updated = 0, 0

    for drv in new_drivers:
        if drv["driver_id"] in existing_ids:
            for i, existing in enumerate(db["drivers"]):
                if existing["driver_id"] == drv["driver_id"]:
                    db["drivers"][i] = drv
                    updated += 1
                    break
        else:
            db["drivers"].append(drv)
            added += 1

    print(f"\nMerge result: {added} added, {updated} updated, "
          f"{len(db['drivers'])} total drivers")
    return db


def main():
    parser = argparse.ArgumentParser(
        description="Scrape loudspeaker T-S parameters from loudspeakerdatabase.com"
    )
    parser.add_argument(
        "--db", type=str, default="data/drivers.json",
        help="Path to driver database JSON (default: data/drivers.json)",
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
        help="Comma-separated list of manufacturers to scrape (e.g. Eminence,B&C)",
    )
    args = parser.parse_args()

    mfr_filter = None
    if args.manufacturers:
        mfr_filter = [m.strip() for m in args.manufacturers.split(",")]

    drivers = scrape_all(
        max_manufacturers=args.max_manufacturers,
        manufacturer_filter=mfr_filter,
        delay=args.delay,
    )
    print(f"\nScraped {len(drivers)} drivers total.")

    if args.dry_run:
        print(json.dumps(drivers, indent=2))
        return

    db_path = Path(args.db)
    db = merge_into_db(db_path, drivers)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    db_path.write_text(json.dumps(db, indent=4) + "\n")
    print(f"Database written to {db_path}")


if __name__ == "__main__":
    main()
