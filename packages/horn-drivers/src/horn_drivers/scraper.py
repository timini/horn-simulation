"""Scrape compression driver T-S parameters from loudspeakerdatabase.com.

Idempotent: re-running merges new drivers into the existing database,
keyed by a ``manufacturer-model`` slug.  All values are converted to SI
on extraction (mH -> H, cm^2 -> m^2, mm -> m).

Usage:
    horn-scrape-drivers [--db data/drivers.json] [--max-pages 5]
"""

import argparse
import json
import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional


BASE_URL = "https://www.loudspeakerdatabase.com"
SEARCH_URL = f"{BASE_URL}/search"


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


def scrape_driver_page(url: str, session) -> Optional[dict]:
    """Scrape a single driver detail page for T-S parameters."""
    from bs4 import BeautifulSoup

    resp = session.get(url, timeout=15)
    if resp.status_code != 200:
        print(f"  WARN: HTTP {resp.status_code} for {url}")
        return None

    soup = BeautifulSoup(resp.text, "html.parser")

    params: Dict[str, Optional[float]] = {}
    param_map = {
        "fs": "fs_hz",
        "re": "re_ohm",
        "bl": "bl_tm",
        "sd": "sd_cm2",
        "mms": "mms_kg",
        "le": "le_mh",
        "qms": "qms",
        "qes": "qes",
        "qts": "qts",
        "xmax": "xmax_mm",
    }

    for row in soup.find_all("tr"):
        cells = row.find_all(["td", "th"])
        if len(cells) < 2:
            continue
        label_text = cells[0].get_text(strip=True).lower()
        value_text = cells[1].get_text(strip=True)

        for key, field in param_map.items():
            if key in label_text:
                val = _parse_float(value_text)
                if val is not None:
                    params[field] = val
                break

    if not params.get("fs_hz"):
        return None

    si = {
        "fs_hz": params.get("fs_hz", 0.0),
        "re_ohm": params.get("re_ohm", 0.0),
        "bl_tm": params.get("bl_tm", 0.0),
        "mms_kg": params.get("mms_kg", 0.0),
        "qms": params.get("qms"),
        "qes": params.get("qes"),
        "qts": params.get("qts"),
    }
    if params.get("sd_cm2"):
        si["sd_m2"] = params["sd_cm2"] * 1e-4
    if params.get("le_mh"):
        si["le_h"] = params["le_mh"] * 1e-3
    if params.get("xmax_mm"):
        si["xmax_m"] = params["xmax_mm"] * 1e-3

    return {k: v for k, v in si.items() if v is not None}


def scrape_compression_drivers(
    max_pages: int = 5,
    delay: float = 1.0,
) -> List[dict]:
    """Search for compression drivers and scrape their parameters."""
    import requests

    session = requests.Session()
    session.headers["User-Agent"] = (
        "Mozilla/5.0 (horn-simulation-research; +https://github.com/timini/horn)"
    )

    drivers = []

    for page in range(1, max_pages + 1):
        print(f"Fetching search page {page}...")
        try:
            resp = session.get(
                SEARCH_URL,
                params={"q": "compression driver", "page": str(page)},
                timeout=15,
            )
        except requests.RequestException as e:
            print(f"  ERROR: {e}")
            break

        if resp.status_code != 200:
            print(f"  WARN: HTTP {resp.status_code}, stopping.")
            break

        from bs4 import BeautifulSoup

        soup = BeautifulSoup(resp.text, "html.parser")
        links = soup.find_all("a", href=re.compile(r"/speaker/"))
        if not links:
            print("  No more results.")
            break

        for link in links:
            href = link.get("href", "")
            name = link.get_text(strip=True)
            if not href or not name:
                continue

            url = href if href.startswith("http") else BASE_URL + href
            print(f"  Scraping: {name} ({url})")
            time.sleep(delay)

            try:
                params = scrape_driver_page(url, session)
            except Exception as e:
                print(f"    ERROR: {e}")
                continue

            if params is None:
                print("    SKIP: missing essential params")
                continue

            parts = name.split(maxsplit=1)
            manufacturer = parts[0] if parts else "Unknown"
            model = parts[1] if len(parts) > 1 else name

            driver_id = slugify(manufacturer, model)
            drivers.append({
                "driver_id": driver_id,
                "manufacturer": manufacturer,
                "model_name": model,
                "driver_type": "compression",
                "parameters": params,
            })
            print(f"    OK: {driver_id}")

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
    parser = argparse.ArgumentParser(description="Scrape compression driver T-S params.")
    parser.add_argument(
        "--db", type=str, default="data/drivers.json",
        help="Path to driver database JSON (default: data/drivers.json)",
    )
    parser.add_argument(
        "--max-pages", type=int, default=5,
        help="Maximum search result pages to scrape (default: 5)",
    )
    parser.add_argument(
        "--delay", type=float, default=1.0,
        help="Delay between requests in seconds (default: 1.0)",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print scraped data without writing to DB",
    )
    args = parser.parse_args()

    drivers = scrape_compression_drivers(max_pages=args.max_pages, delay=args.delay)
    print(f"\nScraped {len(drivers)} drivers.")

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
