"""Fetch pinned FRD/ZMA and acoustic-impedance archives without altering values.

Downloads stay outside the source tree. Import is evidence preparation, not a
physical-validation pass. Archive executables, workbooks and macros are ignored.
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
        if not record.get('download_url') or record.get('format') not in importers:
            continue
        name = record['archive_name']
        cached = Path(a.cache_dir)/name if a.cache_dir else out/name
        archive = cached.read_bytes() if cached.exists() else urlopen(record['download_url'], timeout=60).read()
        result = importers[record['format']](archive, record, out/record['id'])
        (out/name).write_bytes(archive)
        for curve in result['curves']:
            key = curve['source_sha256']
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
