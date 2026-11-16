"""Merge validate_horn_reflectance.py runs into the committed summary JSON.

The base run supplies the protocol, gates, comparisons and diagnostics. Each
FEM run supplies only its "fem" block, keyed by a label such as mesh_4mm. A FEM
run must share the base run's nominal predictions and imported manifest, so the
merged file never mixes inputs. Values are copied unrounded.
"""
import argparse
import json
from pathlib import Path

CATALOG_ID = 'hajicek-neely-horn-reflectance'


def load(path):
    return json.loads((Path(path)/'validation.json').read_text())


def summarise(base_dir, fem_dirs, catalog):
    base = load(base_dir)
    record = next(r for r in catalog['references'] if r['id'] == CATALOG_ID)
    protocol = base['protocol']
    fem = {}
    for label, path in fem_dirs:
        run = load(path)
        if run.get('fem') is None:
            raise ValueError(f'{path} has no FEM block')
        for key in ('imported_manifest_sha256', 'prediction_sha256'):
            if run['protocol'][key] != protocol[key]:
                raise ValueError(f'{path} {key} differs from the base run')
        fem[label] = {**run['fem'], 'run_code_sha256': run['protocol']['code_sha256']}
    failed = any(len(set(count.split('/'))) > 1 for count in base['gate_summary'].values())
    failed |= any(not f['vs_same_model_tmm']['passed'] for f in fem.values())
    return {
        'reference_id': CATALOG_ID, 'source_url': record['source_url'], 'doi': record['doi'],
        'license': record['reuse'].split(',')[0],
        'status': 'frozen_comparison_failed_gates' if failed else 'frozen_comparison_passed_gates',
        'physical_horn_workflow_validated': base['physical_horn_workflow_validated'],
        'script': 'scripts/validate_horn_reflectance.py',
        'protocol': protocol, 'prediction_sha256': protocol['prediction_sha256'],
        'gate_summary': base['gate_summary'], 'probe_plane_estimates': base['probe_plane_estimates'],
        'numerical_checks': base['numerical_checks'], 'comparisons': base['comparisons'],
        'fem': fem or None, 'post_hoc_notes': base.get('post_hoc_notes', []),
        'local_results': {'base': str(base_dir), **{label: str(path) for label, path in fem_dirs}},
        'scope': base['scope'],
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base-run', required=True, help='validate_horn_reflectance.py output directory')
    p.add_argument('--fem-run', action='append', default=[], metavar='LABEL=DIR')
    p.add_argument('--catalog', default='data/validation/references.json')
    p.add_argument('--output', default='data/validation/hajicek_horn_reflectance.json')
    a = p.parse_args()
    fem_dirs = [tuple(item.split('=', 1)) for item in a.fem_run]
    result = summarise(a.base_run, fem_dirs, json.loads(Path(a.catalog).read_text()))
    Path(a.output).write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
