"""Run provenance must not consume the run's own generated files."""
import importlib.util
from pathlib import Path
import subprocess


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
