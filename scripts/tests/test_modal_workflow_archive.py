"""Check that archived automatic results retain their completion-time seals."""
import hashlib
import json
from pathlib import Path
import tarfile

DATA = Path(__file__).resolve().parents[2]/'data/validation'


def test_modal_workflow_and_empty_outcome_preserve_complete_evidence():
    identity = json.loads((DATA/'worked_example_modal_800_1600_manifest.json').read_text())
    for filename, digest in identity['files'].items():
        assert hashlib.sha256((DATA/filename).read_bytes()).hexdigest() == digest
    for filename, prefix, empty in (
        ('worked_example_modal_800_1600.tar.gz', 'worked-example-modal/', False),
        ('modal_empty_workflow.tar.gz', 'modal-empty/', True),
    ):
        with tarfile.open(DATA/filename) as archive:
            def read(name): return archive.extractfile(prefix+name).read()
            manifest = json.loads(read('manifest.json'))
            assert manifest['source_revision'] == identity['reproduction_source_commit']
            assert manifest['status'] == 'completed' and manifest['exit_code'] == 0
            assert manifest['output_sha256']
            for name, digest in manifest['output_sha256'].items():
                assert hashlib.sha256(read(name)).hexdigest() == digest
            ranking = json.loads(read('outputs/auto/report/auto_ranking.json'))
            if empty:
                assert ranking == []
                assert not any(name.endswith('_results.csv') for name in manifest['output_sha256'])
            else:
                assert ranking[0]['radiation_model'] == 'modal_baffled'
                assert ranking[0]['validation_status'] == 'experimental_prediction'
                assert abs(ranking[0]['passband_ripple_db']-1.4847110595181192) < 1e-9
                assert b'nonuniform axisymmetric aperture velocity' in read('outputs/auto/report/auto_report.html')
                assert b'cached' in read('resume.log')
