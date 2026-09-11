"""Refresh a selected nominal size from REDCATT's public engineering API.

Defaults to the 6.5-inch research subset. API license permits attributed
research/comparison, not bulk republication as a competing catalogue.
"""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import requests
from horn_drivers.catalogue import redcatt_record
from horn_drivers.loader import _driver_from_dict
from horn_drivers.scraper import _atomic_json


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--db',type=Path,default=Path('data/drivers'))
    ap.add_argument('--size',type=float,default=6.5)
    a=ap.parse_args()
    with requests.Session() as session:
        schema=session.get('https://www.redcatt.net/api/v1/products/schema',timeout=30);schema.raise_for_status()
        if not schema.json().get('api_version','').startswith('1.'):
            raise ValueError('Review changed API unit contract before importing')
        result=session.get('https://www.redcatt.net/api/v1/products',params={'per_page':500,'size':a.size,'category':'LF'},timeout=30)
        result.raise_for_status();payload=result.json()
    if payload.get('links',{}).get('next'):
        raise ValueError('Incomplete result; pagination must be handled before writing')
    products=payload['data']
    if not products or any(p.get('size_inches')!=a.size for p in products):
        raise ValueError('Source did not return the requested size')
    rows=[redcatt_record(p,datetime.now(timezone.utc).date().isoformat()) for p in products]
    if len({r['driver_id'] for r in rows})!=len(rows):raise ValueError('Duplicate source ordering codes')
    for r in rows:_driver_from_dict(r)
    # No files change until every returned record has passed conversion.
    for r in rows:_atomic_json(a.db/'REDCATT'/(r['driver_id']+'.json'),r)
    print(f'Updated {len(rows)} attributed REDCATT research records. Run catalogue audit afterward.')
    print('No deletions inferred from this filtered feed; periodically review the complete lifecycle feed.')

if __name__=='__main__':main()
