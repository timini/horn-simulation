"""Backward-compatibility shim — driver loading now lives in horn_drivers.loader.

Import from horn_drivers.loader directly for new code.
"""


def load_drivers(db_path):
    from horn_drivers.loader import load_drivers as _load
    return _load(db_path)


def load_driver(db_path, driver_id):
    from horn_drivers.loader import load_driver as _load
    return _load(db_path, driver_id)
