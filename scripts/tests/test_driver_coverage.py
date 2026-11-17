"""Tests for scripts/driver_coverage.py on a small made-up repository tree."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import driver_coverage  # noqa: E402

PARAMS = {"fs_hz": 50.0, "re_ohm": 5.5, "qes": 0.3, "qms": 5.0, "qts": 0.28, "sd_m2": 0.05,
          "mms_kg": 0.06, "bl_tm": 18.0, "le_h": 0.001, "nominal_impedance_ohm": 8.0}


def write(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(data if isinstance(data, str) else json.dumps(data))


def make_tree(root: Path) -> None:
    write(root / "data/sources.json", {"sources": [
        {"id": "acme", "kind": "manufacturer", "base_url": "https://www.acme.example",
         "approval_status": "approved", "status_note": ""},
        {"id": "legacy", "kind": "aggregator", "base_url": "https://legacy.example",
         "approval_status": "approved", "status_note": ""},
        {"id": "shop", "kind": "retailer", "base_url": "https://shop.example",
         "approval_status": "blocked", "status_note": "403 on the first two requests"},
    ]})
    prov = {"source_id": "acme", "url": "https://www.acme.example/a10", "retrieved_on": "2026-11-17",
            "method": "page"}
    v2 = {"driver_id": "acme-a10", "record_schema_version": 2, "manufacturer": "Acme", "model_name": "A10",
          "driver_type": "cone", "parameters": dict(PARAMS, power_w=300.0, sensitivity_db_1w1m=95.0),
          "field_provenance": {k: dict(prov) for k in list(PARAMS) + ["power_w", "sensitivity_db_1w1m"]}}
    v2["field_provenance"]["fs_hz"]["conflicts"] = [{"source_id": "legacy", "value": 48.0}]
    v1 = {"driver_id": "acme-a12", "record_schema_version": 1, "manufacturer": "Acme", "model_name": "A12",
          "driver_type": "cone", "parameters": dict(PARAMS), "nominal_diameter": "12in",
          "parameter_source": "https://legacy.example/Acme/A12", "provenance": {"retrieved_on": "2026-02-26"}}
    write(root / "data/drivers/Acme/acme-a10.json", v2)
    write(root / "data/drivers/Acme/acme-a12.json", v1)
    write(root / "data/prices/shop/2026-11-17.jsonl",
          json.dumps({"driver_id": "acme-a10", "source_id": "shop", "price": 100.0, "currency": "GBP",
                      "region": "UK", "vat_included": True, "match_confidence": "exact",
                      "retrieved_on": "2026-11-17"}) + "\n")


def test_report_counts_records_fields_prices_and_blocked_sources(tmp_path):
    make_tree(tmp_path)
    text = driver_coverage.build_report(tmp_path)
    assert "Newest data date: 2026-11-17" in text
    assert "| Acme | 1 | 1 | 2 |" in text  # one record from acme, one from legacy
    assert "| Power (rated) | 1 | 50 % | 1 | acme 1 |" in text
    assert "| UK | 1 | 1 |" in text
    assert "Stored conflicts (old values kept under `field_provenance.<field>.conflicts`): 1 on 1 records." in text
    assert "| shop | retailer | blocked | 403 on the first two requests |" in text
    assert chr(0x2014) not in text and chr(0x2013) not in text


def test_main_writes_the_same_file_twice(tmp_path):
    make_tree(tmp_path)
    driver_coverage.main(["--root", str(tmp_path), "--out", "docs/DRIVER_COVERAGE.md"])
    first = (tmp_path / "docs/DRIVER_COVERAGE.md").read_text()
    driver_coverage.main(["--root", str(tmp_path), "--out", "docs/DRIVER_COVERAGE.md"])
    assert (tmp_path / "docs/DRIVER_COVERAGE.md").read_text() == first
