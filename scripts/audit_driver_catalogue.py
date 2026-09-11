"""Audit source integrity; retain suspect records with an explicit quarantine.

Run with --apply to annotate shared catalogue and its legacy subset. Reports
are separate from numerical validation, and do not claim unreviewed rows good.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
from horn_drivers.catalogue import audit_records
from horn_drivers.scraper import _atomic_json


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--db',type=Path,default=Path('data/drivers'))
    ap.add_argument('--legacy',type=Path,default=Path('data/drivers.json'))
    ap.add_argument('--output',type=Path,default=Path('data/catalogue-audit.json'))
    ap.add_argument('--apply',action='store_true');a=ap.parse_args()
    paths=sorted(a.db.glob('*/*.json'));records=[json.loads(p.read_text()) for p in paths]
    groups=audit_records(records)
    if a.apply:
        for p,r in zip(paths,records):_atomic_json(p,r)
        if a.legacy.exists():
            old=json.loads(a.legacy.read_text());byid={r['driver_id']:r for r in records}
            old['drivers']=[byid.get(r['driver_id'],r) for r in old['drivers']]
            audit_records(old['drivers']);_atomic_json(a.legacy,old)
    report=dict(record_count=len(records),statuses=dict(Counter(r.get('catalogue_status','unknown') for r in records)),
                repeated_parameter_groups=groups,verified_nominal_size_count=sum(r.get('nominal_diameter_verified') is True for r in records),
                note='Automated identity anomaly audit. Legacy-unverified does not mean validated. Quarantined rows cannot load for simulation.')
    _atomic_json(a.output,report);print(json.dumps({k:v for k,v in report.items() if k!='repeated_parameter_groups'},indent=2))

if __name__=='__main__':main()
