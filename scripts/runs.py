#!/usr/bin/env python3
"""List preserved runs or print the latest completed run's absolute path."""
import argparse
from datetime import datetime
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
STATUSES = ('completed', 'running', 'failed')


def parse_timestamp(value):
    if not isinstance(value, str):
        raise ValueError('Run timestamp must be a string')
    # Python 3.10 needs an explicit UTC offset; 3.11 also accepts Z.
    stamp = datetime.fromisoformat(value.replace('Z', '+00:00'))
    if stamp.tzinfo is None:
        raise ValueError('Run timestamp has no timezone')
    return stamp


def inventory(roots):
    """Read one level under each run root; never infer success from file names."""
    records = {}
    for root in roots:
        root = Path(root).resolve()
        if not root.is_dir():
            raise ValueError(f'Run root is not a directory: {root}')
        for path in sorted(root.iterdir()):
            if not path.is_dir() or path.is_symlink():
                continue
            path = path.resolve()
            row = {'path': str(path), 'status': 'unmanaged', 'started_at': None, 'finished_at': None}
            manifest = path/'manifest.json'
            if manifest.exists():
                try:
                    data = json.loads(manifest.read_text())
                    if not isinstance(data, dict) or data.get('status') not in STATUSES:
                        raise ValueError('Unrecognized run manifest')
                    for field in ('started_at', 'finished_at'):
                        value = data.get(field)
                        if value is not None:
                            parse_timestamp(value)
                        row[field] = value
                    if row['started_at'] is None or (data['status'] != 'running' and row['finished_at'] is None):
                        raise ValueError('Missing run timestamp')
                    row['status'] = data['status']
                except (OSError, ValueError, TypeError) as error:
                    row['status'] = 'invalid_manifest'
                    row['error'] = str(error)
            records[str(path)] = row
    return list(records.values())


def latest(records, status='completed'):
    candidates = [r for r in records if r['status'] in STATUSES and (status == 'any' or r['status'] == status)]
    if not candidates:
        raise ValueError(f'No {status} runs found')
    def key(row):
        # Completed/failed runs are ordered by finish time; running by start.
        stamp = row['finished_at'] if row['status'] != 'running' else row['started_at']
        return parse_timestamp(stamp).timestamp(), row['path']
    return max(candidates, key=key)['path']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('list', 'latest'))
    parser.add_argument('--root', type=Path, action='append', help='Parent of run folders; repeat for custom locations')
    parser.add_argument('--status', choices=(*STATUSES, 'any'), default='completed', help='Latest selector only; default completed')
    args = parser.parse_args()
    try:
        default_root = ROOT/'results'
        rows = [] if args.root is None and not default_root.exists() else inventory(args.root or [default_root])
        print(json.dumps(rows, indent=2) if args.action == 'list' else latest(rows, args.status))
    except ValueError as error:
        parser.error(str(error))
    return 0


if __name__ == '__main__':
    sys.exit(main())
