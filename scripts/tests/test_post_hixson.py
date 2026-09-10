"""Measured component evidence must stay distinct from assembly qualification."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

path = Path(__file__).resolve().parents[1]/'validate_post_hixson.py'
spec = importlib.util.spec_from_file_location('post_hixson', path)
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)


def test_no_fitted_amplitude_or_phase_shift():
    result = v.metrics(np.array([2j, 2j]), np.array([1+0j, 1+0j]))
    assert result['p95_magnitude_db'] == pytest.approx(20*np.log10(2))
    assert result['p95_phase_deg'] == pytest.approx(90)
    assert result['max_complex_normalized'] == pytest.approx(np.sqrt(5))


@pytest.mark.parametrize('bad', [0., np.nan, np.inf])
def test_invalid_reference_cannot_pass(bad):
    with pytest.raises(ValueError, match='Invalid impedance'):
        v.metrics(np.array([1+0j]), np.array([complex(bad)]))


def test_archived_measurement_comparison_replays_from_raw_curves():
    import hashlib
    import io
    import json
    import tarfile
    import pandas as pd

    directory = v.ROOT/'data/validation'
    manifest = json.loads((directory/'post_hixson_manifest.json').read_text())
    for name, digest in manifest['files'].items():
        assert v.sha(directory/name) == digest
    result = json.loads((directory/'post_hixson_reference.json').read_text())
    assert result['physical_assembly_qualified'] is False
    with tarfile.open(directory/'post_hixson_artifacts.tar.gz') as archive:
        def read(name):
            return archive.extractfile('post-hixson/'+name).read()
        protocol = json.loads(read('protocol.json'))
        assert protocol['source_revision'] == manifest['reproduction_source_commit']
        assert protocol['limits'] == v.LIMITS
        for name, digest in json.loads(read('outputs.json')).items():
            assert hashlib.sha256(read(name)).hexdigest() == digest
        # Independently rebuild complex measured impedance and error percentiles.
        components = []
        for component in ('resistance', 'reactance'):
            data = np.loadtxt(io.BytesIO(read(f'measured_{component}.csv')), delimiter=',', skiprows=1)
            assert len(data) > 400
            components.append(np.interp(v.KA, data[:, 0], data[:, 1]))
        reference = components[0]+1j*components[1]
        for case in v.CASES:
            frame = pd.read_csv(io.BytesIO(read(case+'.csv')))
            assert len(frame) == 81
            assert frame.radiation_model.eq('modal_baffled').all()
            assert frame.bc_mode.eq('velocity').all()
            np.testing.assert_allclose(frame.frequency*2*np.pi*.271/343, v.KA)
            assert np.isfinite(frame[['z_real','z_imag','relative_residual']]).all().all()
            assert frame.z_real.gt(0).all()
            assert frame.relative_residual.max() < 1e-8
            predicted = (frame.z_real.to_numpy()+1j*frame.z_imag.to_numpy())/(1.225*343)
            level = np.abs(20*np.log10(abs(predicted)/abs(reference)))
            phase = np.abs(np.angle(predicted/reference, deg=True))
            assert result['comparisons'][case]['p95_magnitude_db'] == pytest.approx(np.percentile(level, 95))
            assert result['comparisons'][case]['p95_phase_deg'] == pytest.approx(np.percentile(phase, 95))
    assert result['measurement_comparison_pass']
    assert result['mesh_refinement_pass']
