from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

import numpy as np


class FlareProfile(str, Enum):
    """Enumeration for the horn flare profile types."""
    EXPONENTIAL = "exponential"
    HYPERBOLIC = "hyperbolic"
    CONICAL = "conical"


@dataclass
class HornParameters:
    """
    A structured representation of a horn's geometric parameters.
    All dimensions are in meters.
    """
    flare_profile: FlareProfile
    throat_radius: float
    mouth_radius: float
    length: float

    # Flare constant for exponential or hyperbolic horns
    m: float = 1.0


@dataclass
class DriverParameters:
    """Thiele-Small parameters for a compression driver. All values in SI units."""

    driver_id: str
    manufacturer: str
    model_name: str

    # Essential T-S parameters (SI)
    fs_hz: float        # Resonance frequency (Hz)
    re_ohm: float       # DC voice coil resistance (Ohm)
    bl_tm: float        # Force factor (T·m)
    sd_m2: float        # Effective diaphragm area (m²)
    mms_kg: float       # Moving mass (kg)
    le_h: float         # Voice coil inductance (H)

    # Q factors — provide at least two; third is derived in __post_init__
    qms: Optional[float] = None
    qes: Optional[float] = None
    qts: Optional[float] = None

    # Derived mechanical params (computed in __post_init__)
    cms_m_per_n: Optional[float] = field(default=None, repr=False)
    rms_kg_per_s: Optional[float] = field(default=None, repr=False)

    # Optional metadata
    driver_type: Optional[str] = None
    xmax_m: Optional[float] = None
    nominal_impedance_ohm: Optional[float] = None

    def __post_init__(self):
        omega_s = 2.0 * np.pi * self.fs_hz

        # Derive Cms from Mms and fs: Cms = 1 / (Mms * omega_s^2)
        if self.cms_m_per_n is None:
            self.cms_m_per_n = 1.0 / (self.mms_kg * omega_s ** 2)

        # Derive Q factors — need at least two of three
        q_count = sum(x is not None for x in [self.qms, self.qes, self.qts])
        if q_count >= 2:
            if self.qts is None:
                self.qts = (self.qms * self.qes) / (self.qms + self.qes)
            elif self.qes is None:
                self.qes = (self.qms * self.qts) / (self.qms - self.qts)
            elif self.qms is None:
                self.qms = (self.qes * self.qts) / (self.qes - self.qts)
        elif q_count == 0:
            # Derive Qes from electrical/mechanical parameters
            if self.bl_tm > 0:
                self.qes = (omega_s * self.mms_kg * self.re_ohm) / (self.bl_tm ** 2)
        elif q_count == 1:
            # With one Q and electrical params, derive others
            if self.qes is None and self.bl_tm > 0:
                self.qes = (omega_s * self.mms_kg * self.re_ohm) / (self.bl_tm ** 2)
            if self.qms is not None and self.qes is not None and self.qts is None:
                self.qts = (self.qms * self.qes) / (self.qms + self.qes)
            elif self.qts is not None and self.qes is not None and self.qms is None:
                self.qms = (self.qes * self.qts) / (self.qes - self.qts)
            elif self.qts is not None and self.qms is not None and self.qes is None:
                self.qes = (self.qms * self.qts) / (self.qms - self.qts)

        # Derive Rms from Qms: Rms = omega_s * Mms / Qms
        if self.rms_kg_per_s is None and self.qms is not None and self.qms > 0:
            self.rms_kg_per_s = omega_s * self.mms_kg / self.qms
