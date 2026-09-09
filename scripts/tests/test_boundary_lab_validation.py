"""Fail-closed checks for the independent numerical comparison protocol."""
import json
from pathlib import Path
import sys
import numpy as np
import pytest

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
import validate_boundary_lab as validation


def test_reference_grid_requires_complete_exact_frequencies():
    expected = np.geomspace(250,2000,13)
    validation.require_grid(expected, expected)
    for actual in (expected[:-1], expected[::-1], expected*1.001, np.full(13,np.nan)):
        with pytest.raises(ValueError, match='frequency grid'):
            validation.require_grid(actual, expected)


def test_complex_comparison_detects_phase_and_rejects_nonfinite():
    assert validation.relative_error([1j], [1]) == pytest.approx(np.sqrt(2))
    for a,b in (([],[]),([1],[1,2]),([np.nan],[1]),([1],[np.inf])):
        with pytest.raises(ValueError, match='finite nonempty'):
            validation.relative_error(a,b)


def test_changed_case_inputs_and_protocol_are_rejected(tmp_path):
    protocol = json.loads(validation.PROTOCOL.read_text())
    for name,_,_ in validation.case_definitions(protocol):
        (tmp_path/name).mkdir()
        for filename in validation.INPUT_NAMES:
            (tmp_path/name/filename).write_text('frozen input')
    manifest = dict(definition=protocol, inputs=validation.protocol_inputs(tmp_path,protocol))
    validation.write_json(tmp_path/'protocol.json',manifest)
    validation.verify_inputs(tmp_path)
    name = next(validation.case_definitions(protocol))[0]
    (tmp_path/name/'horn.msh').write_text('changed input')
    with pytest.raises(ValueError,match='Inputs or protocol'):
        validation.verify_inputs(tmp_path)
    (tmp_path/name/'horn.msh').write_text('frozen input')
    manifest['definition']['voltage_rms'] = 1.
    validation.write_json(tmp_path/'protocol.json',manifest)
    with pytest.raises(ValueError,match='Inputs or protocol'):
        validation.verify_inputs(tmp_path)


def test_changed_or_empty_stage_evidence_is_rejected(tmp_path):
    output = tmp_path/'response.csv'
    output.write_text('frozen response')
    validation.seal_stage(tmp_path,'horn',[output])
    validation.verify_stage(tmp_path,'horn')
    output.write_text('changed response')
    with pytest.raises(ValueError,match='evidence changed'):
        validation.verify_stage(tmp_path,'horn')
    validation.seal_stage(tmp_path,'horn',[])
    with pytest.raises(ValueError,match='evidence changed'):
        validation.verify_stage(tmp_path,'horn')


@pytest.mark.skipif(sys.platform == 'win32', reason='External runner requires POSIX')
@pytest.mark.parametrize('parent_waits', [False, True])
def test_external_runner_reaps_worker_after_parent_exit_or_timeout(tmp_path, parent_waits):
    import os
    import subprocess
    import time
    pid_file = tmp_path/'child.pid'
    child = ("import os,signal,time; from pathlib import Path; "
             "signal.signal(signal.SIGTERM,signal.SIG_IGN); "
             f"Path({str(pid_file)!r}).write_text(str(os.getpid())); time.sleep(60)")
    parent = ("import subprocess,sys,time; from pathlib import Path; "
              f"subprocess.Popen([sys.executable,'-c',{child!r}]); "
              f"p=Path({str(pid_file)!r}); "
              "\nwhile not p.exists(): time.sleep(.01)\n" +
              ("time.sleep(60)" if parent_waits else ""))
    try:
        if parent_waits:
            with pytest.raises(subprocess.TimeoutExpired):
                validation.run_logged([sys.executable,'-c',parent],tmp_path/'log',timeout=2)
        else:
            validation.run_logged([sys.executable,'-c',parent],tmp_path/'log',timeout=10)
        pid = int(pid_file.read_text())
        for _ in range(50):
            state = subprocess.run(['ps','-o','stat=','-p',str(pid)],capture_output=True,text=True).stdout.strip()
            if not state or state.startswith('Z'):
                break
            time.sleep(.02)
        else:
            pytest.fail('External solver worker survived its parent')
    finally:
        if pid_file.exists():
            try:
                os.kill(int(pid_file.read_text()),9)
            except ProcessLookupError:
                pass


def test_committed_reference_archive_and_numeric_fixture_have_recorded_identity():
    import hashlib
    directory=validation.ROOT/'data/validation'
    manifest=json.loads((directory/'boundary_lab_reference_manifest.json').read_text())
    for name,digest in manifest['files'].items():
        assert hashlib.sha256((directory/name).read_bytes()).hexdigest()==digest
    result=json.loads((directory/'boundary_lab_reference.json').read_text())
    assert result['passed'] and len(result['cases'])==6
    assert all(case['passed'] and case['frequencies']==13 for case in result['cases'])
