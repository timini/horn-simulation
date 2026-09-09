"""Run provenance must not consume the run's own generated files."""
import importlib.util
from pathlib import Path
import subprocess
import json
import pytest


def test_custom_run_folder_does_not_invalidate_source_snapshot(tmp_path, monkeypatch):
    script = Path(__file__).resolve().parents[3]/'scripts/run_pipeline.py'
    spec = importlib.util.spec_from_file_location('run_pipeline', script)
    launcher = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(launcher)
    subprocess.run(['git','init',str(tmp_path)],check=True,capture_output=True)
    monkeypatch.setattr(launcher,'ROOT',tmp_path)
    (tmp_path/'source.py').write_text('original source\n')
    run = tmp_path/'custom runs'/'one'
    before = launcher.source_hashes(run)
    run.mkdir(parents=True)
    (run/'manifest.json').write_text('{}')
    (run/'output.csv').write_text('generated result\n')
    assert launcher.source_hashes(run) == before
    # Similar names outside this run are still source inputs.
    neighbor = tmp_path/'custom runs'/'one-more'
    neighbor.mkdir()
    (neighbor/'settings.json').write_text('{}')
    assert launcher.source_hashes(run) != before
    assert 'custom runs/one-more/settings.json' in launcher.source_hashes(run)
    (tmp_path/'source.py').write_text('changed source\n')
    assert launcher.source_hashes(run)['source.py'] != before['source.py']


def _launcher_module():
    script = Path(__file__).resolve().parents[3]/'scripts/run_pipeline.py'
    spec = importlib.util.spec_from_file_location('engine_launcher', script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_engine_identity_detects_version_and_launcher_changes(tmp_path):
    import os
    launcher = _launcher_module()
    executable = tmp_path/'nextflow'
    executable.write_text('#!/bin/sh\necho "version 24.10.5 build 5935"\n')
    executable.chmod(0o755)
    environment = {**os.environ, 'PATH':str(tmp_path)}
    first = launcher.nextflow_identity(environment)
    assert first['version'] == '24.10.5' and first['build'] == '5935'
    executable.write_text('#!/bin/sh\necho "version 24.10.5 build 5935"\n# replacement launcher\n')
    second = launcher.nextflow_identity(environment)
    assert second['launcher_sha256'] != first['launcher_sha256']
    executable.write_text('#!/bin/sh\necho "version 25.10.0 build 6000"\n')
    assert launcher.nextflow_identity(environment)['version'] != first['version']


def test_resume_rejects_changed_engine_before_launch(tmp_path, monkeypatch, capsys):
    import json
    import pytest
    launcher = _launcher_module()
    old_engine = dict(executable='/old/nextflow', launcher_sha256='old', version='24.10.5',build='5935')
    (tmp_path/'manifest.json').write_text(json.dumps({'nextflow_engine':old_engine}))
    monkeypatch.setattr(launcher, 'inspect_images', lambda: {})
    monkeypatch.setattr(launcher, 'source_hashes', lambda *_: {})
    monkeypatch.setattr(launcher, 'java_environment', lambda: {})
    monkeypatch.setattr(launcher, 'nextflow_identity', lambda _: {**old_engine, 'version':'25.10.0'})
    monkeypatch.setattr(launcher.sys, 'argv', ['run_pipeline.py','--run-dir',str(tmp_path),'-resume'])
    monkeypatch.setattr(launcher.subprocess,'call',lambda *a,**kw: pytest.fail('Changed engine must not launch'))
    with pytest.raises(SystemExit) as error:
        launcher.main()
    assert error.value.code == 2
    assert 'Nextflow engine or launcher changed' in capsys.readouterr().err


@pytest.mark.parametrize('change', ['edit', 'remove', 'add', 'missing_seal'])
def test_resume_preserves_original_seal_and_rejects_changed_outputs(tmp_path, monkeypatch, capsys, change):
    launcher = _launcher_module()
    output = tmp_path/'outputs'/'ranking.json'
    output.parent.mkdir()
    output.write_text('[{"score": 1}]')
    engine = dict(version='24.10.5')
    manifest = dict(nextflow_engine=engine, source_sha256={}, containers={}, status='completed',
                    output_sha256=launcher.output_hashes(tmp_path))
    if change == 'edit': output.write_text('[{"score": 2}]')
    elif change == 'remove': output.unlink()
    elif change == 'add': (output.parent/'extra.csv').write_text('replacement')
    else: del manifest['output_sha256']
    path = tmp_path/'manifest.json'
    original = json.dumps(manifest)
    path.write_text(original)
    monkeypatch.setattr(launcher, 'inspect_images', lambda: {})
    monkeypatch.setattr(launcher, 'source_hashes', lambda *_: {})
    monkeypatch.setattr(launcher, 'java_environment', lambda: {})
    monkeypatch.setattr(launcher, 'nextflow_identity', lambda _: engine)
    monkeypatch.setattr(launcher.sys, 'argv', ['run_pipeline.py', '--run-dir', str(tmp_path), '-resume'])
    monkeypatch.setattr(launcher.subprocess, 'call', lambda *a, **kw: pytest.fail('Altered evidence must not launch'))
    with pytest.raises(SystemExit) as error: launcher.main()
    assert error.value.code == 2
    assert 'output' in capsys.readouterr().err.lower()
    assert path.read_text() == original


def test_java_selection_overrides_incompatible_inherited_command(tmp_path, monkeypatch):
    launcher = _launcher_module()
    executable = tmp_path/'bin'/'java'
    executable.parent.mkdir()
    executable.write_text('#!/bin/sh\necho \'openjdk version "21.0.11"\' >&2\n')
    executable.chmod(0o755)
    monkeypatch.setenv('NXF_JAVA_HOME', str(tmp_path))
    monkeypatch.setenv('JAVA_CMD', '/incompatible/java26')
    environment = launcher.java_environment()
    assert environment['JAVA_CMD'] == str(executable.resolve())
    assert environment['JAVA_HOME'] == environment['NXF_JAVA_HOME'] == str(tmp_path.resolve())
