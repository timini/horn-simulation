"""Catalogue source audit and manufacturer API conversion.

Numeric consistency cannot prove product identity. A repeated parameter set
across many products is quarantined pending source verification, not repaired
by guessing a model from Sd. Small families may legitimately share parameters.
"""
from collections import defaultdict
import math
import re

FINGERPRINT_FIELDS = ('fs_hz', 're_ohm', 'bl_tm', 'sd_m2', 'mms_kg')


def repeated_parameter_groups(records, minimum=10):
    groups = defaultdict(list)
    for r in records:
        p = r.get('parameters', {})
        values = [p.get(k) for k in FINGERPRINT_FIELDS]
        if all(isinstance(v, (float, int)) and not isinstance(v, bool) and math.isfinite(v) and v > 0 for v in values):
            groups[(r.get('manufacturer'), *[round(v, 10) for v in values])].append(r['driver_id'])
    return [dict(manufacturer=k[0], parameters=dict(zip(FINGERPRINT_FIELDS, k[1:])),
                 driver_ids=sorted(ids), count=len(ids)) for k, ids in groups.items() if len(ids) >= minimum]


def audit_records(records):
    groups = repeated_parameter_groups(records)
    suspect = {i for g in groups for i in g['driver_ids']}
    for r in records:
        # Repair the known UTF-8-as-Latin-1 display error without altering IDs.
        if isinstance(r.get('model_name'), str):
            r['model_name'] = r['model_name'].replace('Î©', 'Ω')
        if r.get('catalogue_status') == 'manufacturer_verified':
            continue
        if r['driver_id'] in suspect:
            r['catalogue_status'] = 'quarantined'
            r['quarantine_reason'] = 'Mass repeated T/S fingerprint across unrelated products; source identity requires manufacturer verification.'
        else:
            r.setdefault('catalogue_status', 'legacy_unverified')
        r.setdefault('nominal_diameter_verified', False)
    return groups


def _number(value, *, zero=False):
    if isinstance(value, bool):
        raise ValueError('Boolean is not an engineering value')
    if isinstance(value, str) and not re.fullmatch(r'\d+(?:\.\d+)?', value.strip()):
        raise ValueError('Compound or nonnumeric engineering value')
    value = float(value)
    if not math.isfinite(value) or value < 0 or (value == 0 and not zero):
        raise ValueError('Engineering values must be positive and finite')
    return value


def redcatt_record(product, retrieved_on):
    """Convert one active LF product using API v1 units; no AI prose ingestion.

    Missing values stay missing; coax dual ratings cannot become scalar inputs.
    Nominal band is deliberately unknown until a formal response source exists.
    """
    p = product
    if p['category'] != 'LF' or p['status'] != 'active':
        raise ValueError('Only active single-section LF drivers supported')
    ts = p['thiele_small_parameters']
    params = {key: _number(ts[key]) for key in ('fs_hz', 're_ohm', 'bl_tm', 'qms', 'qes', 'qts')}
    params.update(sd_m2=_number(ts['sd_cm2'])*1e-4, mms_kg=_number(ts['mms_grams'])*1e-3,
                  le_h=_number(ts['le_mh'], zero=True)*1e-3,
                  nominal_impedance_ohm=_number(p['general']['nominal_impedance']))
    for key, value, scale in [('power_w',p['power_handling'].get('aes_watts'),1),
                              ('peak_power_w',p['power_handling'].get('peak_watts'),1),
                              ('xmax_m',p['physical'].get('xmax_mm'),.001)]:
        if value is not None:
            params[key] = _number(value)*scale
    code = p['code']
    if not re.fullmatch(r'[A-Za-z0-9.\-]+', code):
        raise ValueError('Invalid ordering code')
    size = _number(p['size_inches'])
    r = dict(driver_id='redcatt-'+code.lower().replace('.','p'), manufacturer='REDCATT',
             model_name=p['name'], ordering_code=code, driver_type='cone',
             nominal_diameter=f'{size:g}in', nominal_diameter_inches=size, nominal_diameter_verified=True,
             catalogue_status='manufacturer_verified', parameters=params, parameter_source=p['url'],
             usable_f_low_hz=None, usable_f_high_hz=None,
             physical_envelope_mm={'maximum_width':None,'depth':p['physical'].get('depth_mm')},
             provenance={'retrieved_on':retrieved_on,'source_updated_at':p['updated_at'],
                         'capture_method':'REDCATT API v1 engineering fields; units per products/schema',
                         'band_definition':'Unknown; free-text AI summaries are not a response specification',
                         'limitations':['Frame maximum including ears must be checked against drawing','No cone transfer function or phase-plug interface measurement']})
    # An overall diameter below the bolt circle cannot describe the full envelope.
    outer = p['physical'].get('overall_diameter_mm'); bolt = p['physical'].get('bolt_circle_diameter_mm')
    if outer is not None and bolt is not None and outer > bolt:
        r['physical_envelope_mm']['maximum_width'] = outer
        r['overall_diameter_m'] = outer*.001
    return r
