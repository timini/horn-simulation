from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

import numpy as np


class FlareProfile(str, Enum):
    """Enumeration for the horn flare profile types."""
    EXPONENTIAL = "exponential"
    HYPERBOLIC = "hyperbolic"
    CONICAL = "conical"
    TRACTRIX = "tractrix"
    OS = "os"
    LECLEACH = "lecleach"
    CD = "cd"


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
    """Thiele-Small parameters for a loudspeaker driver. All values in SI units."""

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

    # Phase plug exit area (m²) — for drivers like 8" cones with phase plugs
    # where the effective throat area is much smaller than Sd
    exit_area_m2: Optional[float] = None

    # Optional metadata
    driver_type: Optional[str] = None
    nominal_diameter: Optional[str] = None  # e.g. "18in", "15in", "1in"
    xmax_m: Optional[float] = None
    nominal_impedance_ohm: Optional[float] = None
    power_w: Optional[float] = None           # RMS / AES / continuous power (W)
    peak_power_w: Optional[float] = None      # Peak / program power (W)
    usable_f_low_hz: Optional[float] = None
    usable_f_high_hz: Optional[float] = None
    overall_diameter_m: Optional[float] = None  # Maximum frame span, including mounting ears
    parameter_source: Optional[str] = None
    interface_model: Optional[str] = None
    mmd_kg: Optional[float] = None  # diaphragm mass without the measured free-air load
    rear_load_mass_kg: Optional[float] = None

    @property
    def coupled_moving_mass_kg(self) -> float:
        """Characterized diaphragm plus rear air load, or disclosed Mms fallback."""
        if self.mmd_kg is not None and self.rear_load_mass_kg is not None:
            return self.mmd_kg + self.rear_load_mass_kg
        return self.mms_kg

    @property
    def effective_throat_area(self) -> float:
        """Return the effective throat area: exit_area_m2 if set, else sd_m2."""
        return self.exit_area_m2 if self.exit_area_m2 is not None else self.sd_m2

    def _validate_optional_numbers(self):
        positive = ("overall_diameter_m", "qms", "qes", "qts", "cms_m_per_n", "exit_area_m2", "xmax_m",
                    "nominal_impedance_ohm", "power_w", "peak_power_w", "usable_f_low_hz",
                    "usable_f_high_hz", "mmd_kg")
        # Zero resistance or rear air mass is a meaningful limiting case.
        nonnegative = ("rms_kg_per_s", "rear_load_mass_kg")
        for name in positive + nonnegative:
            value = getattr(self, name)
            if value is not None and (not np.isfinite(value) or value < 0 or (name in positive and value == 0)):
                requirement = "positive" if name in positive else "nonnegative"
                raise ValueError(f"{self.driver_id}: {name} must be finite and {requirement}")
        if self.usable_f_low_hz is not None and self.usable_f_high_hz is not None and self.usable_f_low_hz >= self.usable_f_high_hz:
            raise ValueError(f"{self.driver_id}: usable frequency band must be increasing")

    def __post_init__(self):
        for name in ("fs_hz", "re_ohm", "bl_tm", "sd_m2", "mms_kg"):
            value = getattr(self, name)
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{self.driver_id}: {name} must be finite and positive")
        if not np.isfinite(self.le_h) or self.le_h < 0:
            raise ValueError(f"{self.driver_id}: le_h must be finite and nonnegative")
        self._validate_optional_numbers()
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
                if self.qms <= self.qts:
                    raise ValueError(f"{self.driver_id}: deriving qes requires qms > qts")
                self.qes = (self.qms * self.qts) / (self.qms - self.qts)
            elif self.qms is None:
                if self.qes <= self.qts:
                    raise ValueError(f"{self.driver_id}: deriving qms requires qes > qts")
                self.qms = (self.qes * self.qts) / (self.qes - self.qts)
        elif q_count == 0:
            # Derive Qes and Qms from electrical/mechanical parameters
            # Qes = (2*pi*fs*Mms*Re) / BL^2
            if self.bl_tm > 0:
                self.qes = (omega_s * self.mms_kg * self.re_ohm) / (self.bl_tm ** 2)
            # Qms requires Rms — if not available, leave None
        elif q_count == 1:
            # With one Q and electrical params, derive others
            if self.qes is None and self.bl_tm > 0:
                self.qes = (omega_s * self.mms_kg * self.re_ohm) / (self.bl_tm ** 2)
            if self.qms is not None and self.qes is not None and self.qts is None:
                self.qts = (self.qms * self.qes) / (self.qms + self.qes)
            elif self.qts is not None and self.qes is not None and self.qms is None:
                self._validate_optional_numbers()
                if self.qes <= self.qts:
                    raise ValueError(f"{self.driver_id}: deriving qms requires qes > qts")
                self.qms = (self.qes * self.qts) / (self.qes - self.qts)
            elif self.qts is not None and self.qms is not None and self.qes is None:
                if self.qms <= self.qts:
                    raise ValueError(f"{self.driver_id}: deriving qes requires qms > qts")
                self.qes = (self.qms * self.qts) / (self.qms - self.qts)

        # Derive Rms from Qms: Rms = omega_s * Mms / Qms
        if self.rms_kg_per_s is None and self.qms is not None and self.qms > 0:
            self.rms_kg_per_s = omega_s * self.mms_kg / self.qms

        self._validate_optional_numbers()
