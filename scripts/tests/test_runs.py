import importlib.util
import json
from pathlib import Path
import pytest

spec = importlib.util.spec_from_file_location('runs',Path(__file__).resolve().parents[1]/'runs.py')
runs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runs)


def record(root, name, status, started, finished=None):
    path = root/name
    path.mkdir()
    (path/'manifest.json').write_text(json.dumps(dict(status=status,started_at=started,finished_at=finished)))
    return str(path.resolve())


def test_latest_completed_excludes_newer_failed_running_and_unmanaged(tmp_path):
    old = record(tmp_path,'old','completed','2026-01-01T10:00:00Z','2026-01-01T11:00:00Z')
    selected = record(tmp_path,'selected','completed','2026-01-01T10:00:00Z','2026-01-01T13:00:00+01:00')
    failed = record(tmp_path,'failed','failed','2026-01-01T13:00:00Z','2026-01-01T14:00:00Z')
    running = record(tmp_path,'running','running','2026-01-01T15:00:00Z')
    (tmp_path/'historical').mkdir()
    rows = runs.inventory([tmp_path])
    assert runs.latest(rows) == selected
    assert runs.latest(rows,'failed') == failed
    assert runs.latest(rows,'running') == running
    assert runs.latest(rows,'any') == running
    assert next(r for r in rows if r['path'].endswith('/historical'))['status'] == 'unmanaged'


def test_malformed_manifest_never_becomes_latest(tmp_path):
    for name, data in [('badjson','{'),('array','[]'),('missingtime','{"status":"completed"}'),
                       ('naive','{"status":"completed","started_at":"2026-01-01","finished_at":"2026-01-02"}')]:
        path=tmp_path/name
        path.mkdir()
        (path/'manifest.json').write_text(data)
    rows=runs.inventory([tmp_path])
    assert all(r['status']=='invalid_manifest' for r in rows)
    with pytest.raises(ValueError,match='No completed runs'):
        runs.latest(rows)


def test_multiple_roots_deduplicate_and_do_not_follow_symlinks(tmp_path):
    root=tmp_path/'a'; root.mkdir()
    other=tmp_path/'b'; other.mkdir()
    selected=record(root,'run','completed','2026-01-01T10:00:00Z','2026-01-01T11:00:00Z')
    (other/'link').symlink_to(root/'run',target_is_directory=True)
    rows=runs.inventory([root,root,other])
    assert len(rows)==1 and runs.latest(rows)==selected
    with pytest.raises(ValueError,match='not a directory'):
        runs.inventory([tmp_path/'missing'])


def test_fresh_checkout_has_empty_inventory_and_no_latest(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(runs,'ROOT',tmp_path)
    monkeypatch.setattr(runs.sys,'argv',['runs.py','list'])
    assert runs.main()==0
    assert json.loads(capsys.readouterr().out)==[]
    monkeypatch.setattr(runs.sys,'argv',['runs.py','latest'])
    with pytest.raises(SystemExit) as error:
        runs.main()
    assert error.value.code==2
    assert 'No completed runs found' in capsys.readouterr().err
    monkeypatch.setattr(runs.sys,'argv',['runs.py','list','--root',str(tmp_path/'missing')])
    with pytest.raises(SystemExit) as error:
        runs.main()
    assert error.value.code==2
    assert 'Run root is not a directory' in capsys.readouterr().err
