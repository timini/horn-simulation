import hashlib
from io import BytesIO
from zipfile import ZipFile
import numpy as np
import pandas as pd
import pytest
from horn_analysis.reference_data import read_response, import_archive, import_pipe_archive
from horn_analysis.validation import compare_curves


def test_impedance_domain_and_phase_are_preserved():
    frame, phase = read_response(b'100 8 90\n200 4 -90\n', 'zma')
    assert phase
    np.testing.assert_allclose(frame.z_real_ohm, 0, atol=1e-12)
    np.testing.assert_allclose(frame.z_imag_ohm, [8,-4])
    _, phase = read_response(b'100 80 0\n200 82 0\n', 'frd')
    assert not phase


def test_import_verifies_archive_and_ignores_unrelated_files(tmp_path):
    content = BytesIO()
    with ZipFile(content, 'w') as z:
        z.writestr('horn/FreqResp 7.5.frd', '100 80 0\n200 82 0\n')
        z.writestr('../../escape.exe', 'unrelated file')
    raw = content.getvalue()
    record = {'id':'fixture', 'sha256':hashlib.sha256(raw).hexdigest()}
    report = import_archive(raw, record, tmp_path)
    assert len(report['curves']) == 1
    assert report['curves'][0]['angle_deg'] == 7.5
    assert not report['physical_validation_passed']
    assert len(list(tmp_path.iterdir())) == 1
    with pytest.raises(ValueError, match='checksum'):
        import_archive(raw + b'changed', record, tmp_path)


@pytest.mark.parametrize('data', [b'100 1 0\n100 2 0\n', b'100 1 0\n200 nan 0\n', b'100 1 0\n200 garbage 0\n'])
def test_bad_measurements_are_not_silently_cleaned(data):
    with pytest.raises(ValueError):
        read_response(data, 'frd')


def test_validation_cannot_hide_notch_between_reference_samples():
    ref = pd.DataFrame({'frequency':[100,1000], 'spl':[80,80]})
    prediction = pd.DataFrame({'frequency':[100,300,1000], 'spl':[80,60,80]})
    comparison = compare_curves(ref, prediction, 100, 1000)
    assert comparison['max_absolute_error_db'] == 20
    assert not comparison['level_gate_passed']
    shifted = ref.copy(); shifted.spl += 10
    assert not compare_curves(ref, shifted, 100, 1000)['level_gate_passed']
    relative = compare_curves(ref, shifted, 100, 1000, relative=True)
    assert not relative['physical_validation_passed']
    assert not relative['level_gate_passed']
    assert relative['shape_gate_passed']
    assert relative['mode'] == 'relative_shape'


def test_pipe_import_preserves_measurement_units_and_extensionless_exports(tmp_path):
    content = BytesIO()
    with ZipFile(content, 'w') as z:
        z.writestr('Raw_data/Measured_Impedance/Experimenter_O2/Brass_O/repeat1', '100 -0.01 0.4\n200 0.03 0.8\n')
        z.writestr('Raw_data/Simulated_Impedance/Cylinder_unflanged/reference.txt', '100 1000 4000\n200 3000 8000\n')
        z.writestr('Raw_data/Measured_Impedance/Experimenter_O2/Brass_O/run.py', 'do not execute')
    raw = content.getvalue()
    record = {'id': 'fixture', 'sha256': hashlib.sha256(raw).hexdigest()}
    result = import_pipe_archive(raw, record, tmp_path)
    measured, simulated = result['curves']
    assert len(result['curves']) == 2
    assert measured['reference_kind'] == 'measurement'
    assert simulated['reference_kind'] == 'independent_simulation'
    # Small negative measured resistance is retained as evidence, not clipped.
    assert pd.read_csv(tmp_path/measured['csv']).z_real_normalized.iloc[0] == -0.01
    assert pd.read_csv(tmp_path/simulated['csv']).z_real_pa_s_per_m3.iloc[0] == 1000
    assert not result['physical_validation_passed']
    with pytest.raises(ValueError, match='checksum'):
        import_pipe_archive(raw + b'changed', record, tmp_path)
