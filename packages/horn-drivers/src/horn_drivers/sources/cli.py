"""horn-sources: list sources, snapshot robots.txt, fetch ECB rates, run an adapter.

    horn-sources list
    horn-sources robots [--only bc,thomann-uk]     # fetch and store robots.txt, update snapshot dates
    horn-sources fx                                 # ECB daily rates into data/prices/fx/
    horn-sources run <source_id> [--limit N]        # one approved adapter
    horn-sources coverage [--manufacturers 18Sound,BC] [--driver-type cone]
                                                    # drivers with at least one exact price

All network access goes through PoliteFetcher with the registry rate limit.
"""

import argparse
import datetime as _dt
import json
import sys
from pathlib import Path
from urllib.parse import urlsplit

from horn_drivers.fetch import DEFAULT_CACHE_ROOT, SourceBlocked, _is_html
from horn_drivers.sources.registry import (DEFAULT_REGISTRY, SourceEntry, adapter_ids, create_adapter,
                                           load_source_registry, make_fetcher, save_source_registry)

ROBOTS_DIR = Path("data/sources/robots")


def snapshot_robots(entry: SourceEntry, fetcher, robots_dir: Path = ROBOTS_DIR) -> dict:
    """Fetch robots.txt for one source, store a copy, return what happened."""
    parts = urlsplit(entry.base_url)
    url = entry.robots_url or f"{parts.scheme}://{parts.netloc}/robots.txt"
    result = fetcher.get(url, max_age_days=1.0)
    robots_dir.mkdir(parents=True, exist_ok=True)
    outcome = {"id": entry.id, "url": url, "status": result.status, "date": result.retrieved_on}
    if result.ok and _is_html(result):
        # Some sites answer /robots.txt with their home page. That is not a
        # robots file (no rules apply) and its content is not stored.
        outcome["note"] = "served an HTML page, not a robots file"
    elif result.ok:
        (robots_dir / f"{entry.id}.txt").write_bytes(result.content)
    return outcome



def price_coverage(records, offers, manufacturers=None, driver_type=None) -> dict:
    """Drivers in *records* with at least one exact or manual offer: overall, per region, per source."""
    from horn_drivers.identity import canonical_manufacturer
    from horn_drivers.prices import RANKABLE_MATCHES
    wanted = {canonical_manufacturer(m) or m for m in manufacturers} if manufacturers else None
    pilot = {}
    for record in records:
        maker = canonical_manufacturer(record.get("manufacturer")) or record.get("manufacturer")
        if wanted and maker not in wanted:
            continue
        if driver_type and record.get("driver_type") != driver_type:
            continue
        pilot[record["driver_id"]] = maker
    priced, by_region, by_source, by_maker = set(), {}, {}, {}
    for offer in offers:
        if offer.driver_id not in pilot or offer.match_confidence not in RANKABLE_MATCHES:
            continue
        priced.add(offer.driver_id)
        by_region.setdefault(offer.region, set()).add(offer.driver_id)
        by_source.setdefault(offer.source_id, set()).add(offer.driver_id)
        by_maker.setdefault(pilot[offer.driver_id], set()).add(offer.driver_id)
    makers = {}
    for maker in sorted(set(pilot.values())):
        total = sum(1 for m in pilot.values() if m == maker)
        makers[maker] = {"drivers": total, "priced": len(by_maker.get(maker, ()))}
    regions_per_driver = {}
    for region, ids in by_region.items():
        for driver_id in ids:
            regions_per_driver[driver_id] = regions_per_driver.get(driver_id, 0) + 1
    return {"drivers": len(pilot), "priced": len(priced),
            "per_region": {k: len(v) for k, v in sorted(by_region.items())},
            "per_source": {k: len(v) for k, v in sorted(by_source.items())},
            "per_manufacturer": makers,
            "in_two_or_more_regions": sum(1 for n in regions_per_driver.values() if n >= 2)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--cache-root", type=Path, default=DEFAULT_CACHE_ROOT)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list")
    robots = sub.add_parser("robots")
    robots.add_argument("--only", default=None, help="Comma list of source ids.")
    robots.add_argument("--robots-dir", type=Path, default=ROBOTS_DIR)
    fx = sub.add_parser("fx")
    fx.add_argument("--prices-dir", type=Path, default=Path("data/prices"))
    run = sub.add_parser("run")
    run.add_argument("source_id")
    run.add_argument("--limit", type=int, default=None)
    run.add_argument("--db", type=Path, default=Path("data/drivers"))
    run.add_argument("--incomplete", type=Path, default=Path("data/drivers-incomplete"))
    run.add_argument("--index-db", type=Path, action="append", default=None,
                     help="Driver folders to match retailer listings against (repeatable).")
    run.add_argument("--prices-dir", type=Path, default=Path("data/prices"))
    run.add_argument("--review-csv", type=Path, default=None)
    cover = sub.add_parser("coverage")
    cover.add_argument("--prices-dir", type=Path, default=Path("data/prices"))
    cover.add_argument("--db", type=Path, action="append", default=None)
    cover.add_argument("--manufacturers", default=None, help="Comma list, for example 18Sound,BC,Celestion.")
    cover.add_argument("--driver-type", default=None)
    args = parser.parse_args(argv)

    entries = load_source_registry(args.registry)
    if args.command == "list":
        known = set(adapter_ids())
        for entry in entries.values():
            print(f"{entry.id:24} {entry.kind:12} {entry.approval_status:9} delay {entry.min_delay_s:>4g} s  "
                  f"adapter {'yes' if entry.id in known else 'no '}  {entry.base_url}")
        return 0

    if args.command == "robots":
        wanted = set(args.only.split(",")) if args.only else None
        report = []
        for entry in entries.values():
            if (wanted and entry.id not in wanted) or (not wanted and entry.approval_status == "blocked"):
                continue
            fetcher = make_fetcher(entry, args.cache_root)
            try:
                outcome = snapshot_robots(entry, fetcher, args.robots_dir)
            except SourceBlocked as error:
                outcome = {"id": entry.id, "status": None, "error": str(error)}
            except Exception as error:  # one unreachable host must not stop the others
                outcome = {"id": entry.id, "status": None, "error": f"{type(error).__name__}: {error}"}
            if outcome.get("status") is not None:
                entry.robots_snapshot_date = outcome["date"]
            report.append(outcome)
            print(json.dumps(outcome))
        save_source_registry(entries, args.registry)
        return 0

    if args.command == "coverage":
        from horn_drivers.loader import load_drivers_raw
        from horn_drivers.prices import load_offers
        records = [r for folder in (args.db or [Path("data/drivers"), Path("data/drivers-curated")])
                   for r in load_drivers_raw(str(folder))]
        makers = args.manufacturers.split(",") if args.manufacturers else None
        print(json.dumps(price_coverage(records, load_offers(args.prices_dir), makers, args.driver_type), indent=2))
        return 0

    if args.command == "fx":
        from horn_drivers.prices import fetch_ecb_rates
        entry = entries["ecb"]
        path = fetch_ecb_rates(make_fetcher(entry, args.cache_root), args.prices_dir)
        print(f"Stored {path}")
        return 0

    adapter = create_adapter(args.source_id, args.registry, args.cache_root)
    if adapter.kind == "retailer":
        from horn_drivers.loader import load_drivers_raw
        from horn_drivers.sources.base import DriverIndex
        records = [r for folder in (args.index_db or [Path("data/drivers"), Path("data/drivers-curated")])
                   for r in load_drivers_raw(str(folder))]
        review = args.review_csv or Path("results/price-review") / f"{args.source_id}-{_dt.date.today()}.csv"
        summary = adapter.run(DriverIndex(records), args.prices_dir, review, limit=args.limit)
    else:
        summary = adapter.run(args.db, args.incomplete, limit=args.limit)
    print(json.dumps({**summary.to_dict(), "fetch": adapter.fetcher.status()}, indent=2))
    return 0 if summary.blocked_reason is None else 2


if __name__ == "__main__":
    sys.exit(main())
