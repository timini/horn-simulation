"""Keep the executed independent aperture evidence auditable."""
import importlib.util
import json
from pathlib import Path
import tarfile
import numpy as np
import pytest

path=Path(__file__).resolve().parents[1]/'validate_modal_mouth.py'
spec=importlib.util.spec_from_file_location('modal_mouth_validation',path)
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)


def test_archive_has_complete_raw_sweeps_and_matching_identity():
    directory=v.ROOT/'data/validation'
    manifest=json.loads((directory/'modal_mouth_manifest.json').read_text())
    for name,digest in manifest['files'].items():assert v.sha(directory/name)==digest
    result=json.loads((directory/'modal_mouth_reference.json').read_text())
    assert result['passed'] and len(result['comparisons'])==3
    assert result['physical_validation_status']=='experimental_prediction'
    with tarfile.open(directory/'modal_mouth_artifacts.tar.gz') as archive:
        def read(name):return archive.extractfile('modal-mouth/'+name).read()
        protocol=json.loads(read('protocol.json'));evidence=json.loads(read('solve-evidence.json'))
        assert protocol['source_revision']==manifest['reproduction_source_commit']
        assert protocol['limits']==v.LIMITS
        assert evidence['protocol_sha256']==v.hashlib.sha256(read('protocol.json')).hexdigest()
        for name,digest in {**protocol['inputs'],**evidence['files']}.items():
            assert v.hashlib.sha256(read(name)).hexdigest()==digest
        for name in v.CASES:
            rows=json.loads(read(name+'.json'))
            assert len(rows)==33
            np.testing.assert_allclose([r['frequency'] for r in rows],protocol['frequencies_hz'])


@pytest.mark.parametrize('bad',[0.,np.nan,np.inf])
def test_nonfinite_or_zero_complex_reference_is_rejected(bad):
    with pytest.raises(ValueError,match='Invalid comparison'):v.change(np.array([1.+0j]),np.array([bad+0j]))


def test_magnitude_and_phase_disagreements_are_not_fitted_away():
    delta=v.change(np.array([1.+0j]),np.array([2.j]))
    assert delta['level_db']==pytest.approx(20*np.log10(2))
    assert delta['phase_deg']==pytest.approx(90.)
