"""Fetch pinned FRD/ZMA and acoustic-impedance archives without altering values.

Downloads stay outside the source tree. Import is evidence preparation, not a
physical-validation pass. Archive executables, workbooks, macros and MATLAB
scripts are ignored; they are never executed.
"""
import argparse
import hashlib
from io import BytesIO
import json
from pathlib import Path
import re
from urllib.request import urlopen
from zipfile import ZipFile

import numpy as np
import pandas as pd
from scipy.io import loadmat
from horn_core.acoustics import validate_response


def read_response(data, kind):
    """Read three-column FRD (dB/degrees) or ZMA (ohms/degrees) text.

    Comment headers are allowed; malformed numeric rows fail. A zero-filled
    phase column is recorded as unavailable rather than measured zero phase.
    """
    if kind not in ('frd', 'zma'):
        raise ValueError('Expected frd or zma')
    rows = []
    for number, line in enumerate(data.decode('utf-8-sig', errors='strict').splitlines(), 1):
        line = line.strip()
        if not line or line.startswith(('*', '#', ';')):
            continue
        try:
            row = [float(v) for v in line.split()]
        except ValueError as exc:
            raise ValueError(f'Malformed response at line {number}') from exc
        if len(row) != 3:
            raise ValueError(f'Expected three columns at line {number}')
        rows.append(row)
    if len(rows) < 2:
        raise ValueError('Response needs at least two samples')
    values = np.asarray(rows)
    validate_response(values[:, 0], values[:, 1], values[:, 2])
    phase_available = bool(np.any(values[:, 2] != 0))
    if kind == 'frd':
        frame = pd.DataFrame(values, columns=['frequency', 'spl', 'phase_deg'])
    else:
        if np.any(values[:, 1] <= 0):
            raise ValueError('Impedance magnitude must be positive')
        frame = pd.DataFrame(values, columns=['frequency', 'magnitude_ohm', 'phase_deg'])
        z = values[:, 1]*np.exp(1j*np.deg2rad(values[:, 2]))
        frame['z_real_ohm'], frame['z_imag_ohm'] = z.real, z.imag
    return frame, phase_available


def import_archive(archive, record, output):
    """Verify the source checksum and import only numeric response members."""
    digest = hashlib.sha256(archive).hexdigest()
    if digest != record['sha256']:
        raise ValueError(f"Source checksum changed: {record['id']}")
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    imported = []
    with ZipFile(BytesIO(archive)) as z:
        for name in sorted(z.namelist()):
            suffix = Path(name).suffix.lower()
            if suffix not in ('.frd', '.zma'):
                continue
            raw = z.read(name)
            frame, phase_available = read_response(raw, suffix[1:])
            member_hash = hashlib.sha256(raw).hexdigest()
            filename = re.sub(r'[^A-Za-z0-9_.-]+', '_', name) + '.csv'
            frame.to_csv(output/filename, index=False)
            match = re.search(r'FreqResp ([\d.]+)\.frd$', name)
            imported.append({
                'source_member': name, 'source_sha256': member_hash,
                'csv': filename, 'csv_sha256': hashlib.sha256((output/filename).read_bytes()).hexdigest(),
                'quantity': 'sound_pressure_level' if suffix == '.frd' else 'electrical_impedance',
                'phase_available': phase_available,
                'angle_deg': float(match[1]) if match else None,
                'plane': ('vertical' if '/SE V/' in name else 'horizontal') if match else None,
                'samples': len(frame), 'frequency_range_hz': [float(frame.frequency.iloc[0]), float(frame.frequency.iloc[-1])],
            })
    return {'reference': record, 'curves': imported, 'physical_validation_passed': False}


def import_pipe_archive(archive, record, output):
    """Import the Ernoult benchmark without confusing its two impedance units.

    Measured Z is dimensionless Z/(rho*c/S); simulation Z is Pa*s/m^3.
    Extensionless measurement files are intentional. No clipping, smoothing,
    temperature correction or phase-convention conversion is applied here.
    """
    if hashlib.sha256(archive).hexdigest() != record['sha256']:
        raise ValueError(f"Source checksum changed: {record['id']}")
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    curves = []
    with ZipFile(BytesIO(archive)) as z:
        for name in sorted(z.namelist()):
            parts = Path(name).parts
            if name.endswith('/') or len(parts) < 4 or parts[0] != 'Raw_data':
                continue
            measured = parts[1] == 'Measured_Impedance'
            if not measured and parts[1] != 'Simulated_Impedance':
                continue
            if Path(name).suffix not in ('', '.txt'):
                continue
            raw = z.read(name)
            values = np.loadtxt(BytesIO(raw), ndmin=2)
            if values.shape[1] != 3:
                raise ValueError(f'Expected three impedance columns: {name}')
            validate_response(values[:, 0], values[:, 1], values[:, 2])
            unit = 'normalized' if measured else 'pa_s_per_m3'
            frame = pd.DataFrame(values, columns=['frequency', f'z_real_{unit}', f'z_imag_{unit}'])
            # A path hash prevents distinct archive paths collapsing to one name.
            filename = re.sub(r'[^A-Za-z0-9_.-]+', '_', name) + '.' + hashlib.sha256(name.encode()).hexdigest()[:12] + '.csv'
            frame.to_csv(output/filename, index=False)
            curves.append({
                'source_member': name, 'source_sha256': hashlib.sha256(raw).hexdigest(),
                'csv': filename, 'csv_sha256': hashlib.sha256((output/filename).read_bytes()).hexdigest(),
                'quantity': 'input_acoustic_impedance', 'impedance_units': unit,
                'reference_kind': 'measurement' if measured else 'independent_simulation',
                'configuration': parts[3] if measured else parts[2],
                'operator': parts[2] if measured else None,
                'phase_available': True,
                'samples': len(frame),
                'frequency_range_hz': [float(frame.frequency.iloc[0]), float(frame.frequency.iloc[-1])],
            })
    if not curves:
        raise ValueError('No recognized pipe impedance files')
    return {'reference': record, 'curves': curves, 'physical_validation_passed': False}

THL_CHANNELS = ('zl1', 'zl2')
CGS_ACOUSTIC_OHM = 'cgs_acoustic_ohm'  # dyn*s/cm^5; 1 cgs acoustic ohm = 1e5 Pa*s/m^3


def read_thl(raw):
    """Read a Neely/Hajicek THL measurement (MATLAB level-4 MAT data).

    Returns the uniform frequency grid without its DC bin and the two probe
    channels of measured load impedance in cgs acoustic ohms. The DC bin is not
    a measurement (the authors replace it by a causal estimate), so it is
    omitted rather than repaired. No other value is changed.
    """
    data = loadmat(BytesIO(raw))
    missing = {'df', 'rate', *THL_CHANNELS} - set(data)
    if missing:
        raise ValueError(f'THL file lacks {sorted(missing)}')
    df, rate = float(np.squeeze(data['df'])), float(np.squeeze(data['rate']))
    channels = {name: np.asarray(data[name], dtype=complex).ravel() for name in THL_CHANNELS}
    lengths = {len(v) for v in channels.values()}
    if len(lengths) != 1 or not np.isclose(rate, 2*lengths.pop()*df, rtol=1e-9):
        raise ValueError('THL channels must share a grid ending one bin below Nyquist')
    n = len(channels['zl1'])
    frequency = np.arange(1, n)*df
    channels = {k: v[1:] for k, v in channels.items()}
    validate_response(frequency, *[v.real for v in channels.values()], *[v.imag for v in channels.values()])
    z0 = float(np.ravel(data['z0'])[0]) if 'z0' in data else None
    return frequency, channels, {'df_hz': df, 'sample_rate_hz': rate, 'author_z0_cgs': z0}


def import_thl_record(files, record, output):
    """Import a multi-file Zenodo record whose measurements are THL files.

    ``files`` maps catalog file names to downloaded bytes. Every listed file is
    checked against its pinned SHA-256 before anything is written.
    """
    expected = {f['name']: f for f in record['files']}
    if set(files) != set(expected):
        raise ValueError(f"File set differs from catalog: {record['id']}")
    for name, raw in files.items():
        if hashlib.sha256(raw).hexdigest() != expected[name]['sha256']:
            raise ValueError(f"Source checksum changed: {record['id']}/{name}")
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    curves = []
    for name in sorted(files):
        entry = expected[name]
        if entry.get('role') != 'measurement':
            continue
        frequency, channels, info = read_thl(files[name])
        for channel, z in channels.items():
            filename = re.sub(r'[^A-Za-z0-9_.-]+', '_', f'{name}.{channel}') + '.csv'
            pd.DataFrame({'frequency': frequency, f'z_real_{CGS_ACOUSTIC_OHM}': z.real,
                          f'z_imag_{CGS_ACOUSTIC_OHM}': z.imag}).to_csv(output/filename, index=False)
            curves.append({
                'source_member': name, 'source_sha256': entry['sha256'], 'channel': channel,
                'csv': filename, 'csv_sha256': hashlib.sha256((output/filename).read_bytes()).hexdigest(),
                'quantity': 'input_acoustic_impedance', 'impedance_units': CGS_ACOUSTIC_OHM,
                'reference_kind': 'measurement', 'configuration': entry.get('configuration'),
                'phase_available': True, 'samples': len(frequency),
                'frequency_range_hz': [float(frequency[0]), float(frequency[-1])], **info,
            })
    if not curves:
        raise ValueError('No THL measurement files listed')
    return {'reference': record, 'curves': curves, 'physical_validation_passed': False}


def _fetch(url):
    return urlopen(url, timeout=120).read()


def _load_record_files(record, cache_dir, out):
    """Fetch each catalogued file of a multi-file record, preferring a cache."""
    files = {}
    for entry in record['files']:
        name = entry['name']
        candidates = [Path(cache_dir)/record['id']/name, Path(cache_dir)/name] if cache_dir else []
        cached = next((c for c in candidates if c.exists()), None)
        files[name] = cached.read_bytes() if cached else _fetch(f"{record['download_base_url']}/{name}/content")
    source = out/record['id']/'source'
    source.mkdir(parents=True, exist_ok=True)
    return files, source


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog', required=True)
    p.add_argument('--output-dir', required=True)
    p.add_argument('--cache-dir', help='Reuse previously downloaded verified archives')
    a = p.parse_args()
    catalog = json.loads(Path(a.catalog).read_text())
    out = Path(a.output_dir); out.mkdir(parents=True, exist_ok=True)
    summary = []
    seen = {}
    for record in catalog['references']:
        importers = {'zip_frd_zma': import_archive, 'zip_pipe_impedance': import_pipe_archive}
        if record.get('format') == 'zenodo_files_thl':
            files, source = _load_record_files(record, a.cache_dir, out)
            result = import_thl_record(files, record, out/record['id'])
            for name, raw in files.items():
                (source/name).write_bytes(raw)
        elif not record.get('download_url') or record.get('format') not in importers:
            continue
        else:
            name = record['archive_name']
            cached = Path(a.cache_dir)/name if a.cache_dir else out/name
            archive = cached.read_bytes() if cached.exists() else urlopen(record['download_url'], timeout=60).read()
            result = importers[record['format']](archive, record, out/record['id'])
            (out/name).write_bytes(archive)
        for curve in result['curves']:
            # Two probe channels of one THL file are distinct curves.
            key = (curve['source_sha256'], curve.get('channel'))
            if key in seen:
                curve['duplicate_of'] = seen[key]
            else:
                seen[key] = f"{record['id']}/{curve['source_member']}"
        (out/record['id']/'manifest.json').write_text(json.dumps(result, indent=2))
        summary.append(result)
    (out/'inventory.json').write_text(json.dumps({'datasets': summary, 'unique_curve_files': len(seen), 'physical_validation_passed': False}, indent=2))
    print(f'Imported {sum(len(r["curves"]) for r in summary)} curves ({len(seen)} unique files) from {len(summary)} archives')


if __name__ == '__main__':
    main()
