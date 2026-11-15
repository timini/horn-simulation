"""Driver motor and driver-to-throat interface model (RMS phasors, exp(+i w t)).

Network, from the amplifier to the horn throat::

    V_g -- Ze = Re + j w Le -- gyrator Bl -- Zm = Rms + j w M + 1/(j w Cms)
                                           -- rear chamber: Sd^2 rho c^2 / (j w Vrc)   (series, optional)
                                           -- area transformer Sd (diaphragm volume velocity U_d = Sd v)
                                           -- front chamber two-port (shunt compliance or short duct)
                                           -- horn throat: Z_A = z_horn / S_t     (acoustic, Pa s/m^3)

Front (throat) chamber, volume Vtc between the diaphragm and the horn throat:

* Without a chamber area it is a lumped shunt compliance C = Vtc / (rho c^2)
  (the acoustic compliance of a small enclosed volume; Beranek and Mellow,
  *Acoustics: Sound Fields and Transducers*, 2012, chapter 4). The throat
  volume velocity is U_t = U_d / (1 + j w C Z_A). On a resistive throat
  load Z_A = rho c / S_t the magnitude falls by 3 dB where w C Z_A = 1, that
  is f = c S_t / (2 pi Vtc).
* With a chamber area Atc it is a uniform lossless duct of area Atc and
  length Vtc / Atc (plane-wave transfer matrix). For k Vtc / Atc << 1 this
  reduces to the same shunt compliance plus a small series air mass.

Rear chamber, a sealed volume Vrc behind the diaphragm, adds the air-spring
stiffness Sd^2 rho c^2 / Vrc in series with the suspension. Leaving it unset
keeps the free-air rear load already contained in the published Mms (open
back), which is the previous behaviour.

All interface elements are lossless, so electrical input power equals copper
loss plus suspension loss plus the acoustic power delivered to the throat.

A chamber volume of zero, with no area, uses exactly the previous ideal area
transformer arithmetic, so results are bit-for-bit unchanged.

Cone drivers with no chamber data get an estimated front chamber equal to the
air volume of the cone itself: a right circular cone of base Sd and depth
0.25 times the effective diameter, V = Sd * 0.25 * D_eff / 3. Any adapter
between the cone rim and the throat adds to this, so the estimate is a lower
bound and is reported as "estimated". Compression drivers and drivers of
unknown type are not estimated: their front chamber is treated as zero and
flagged as unknown.
"""
import math

import numpy as np

from horn_core.acoustics import C0, RHO0

CONE_DEPTH_TO_DIAMETER = 0.25
# First radial (0,1) mode of a thin cylindrical cavity, j'_1 zero: k a = 3.8317.
FIRST_RADIAL_CAVITY_ROOT = 3.8317


def estimate_cone_front_chamber_m3(sd_m2: float) -> float:
    """Air volume of a cone of base area Sd and depth 0.25 x effective diameter."""
    diameter = 2.0 * math.sqrt(sd_m2 / math.pi)
    return sd_m2 * CONE_DEPTH_TO_DIAMETER * diameter / 3.0


def resolve_front_chamber(driver):
    """Return (volume_m3, area_m2 or None, source) for the driver's front chamber.

    source is "specified", "estimated" (cone geometric estimate) or "unknown"
    (zero volume used, the ideal area transformer).
    """
    if driver.front_chamber_volume_m3 is not None:
        return float(driver.front_chamber_volume_m3), driver.front_chamber_area_m2, "specified"
    if (driver.driver_type or "").lower() == "cone":
        return estimate_cone_front_chamber_m3(driver.sd_m2), None, "estimated"
    return 0.0, None, "unknown"


def front_chamber_matrix(omega, volume_m3, area_m2=None, *, rho=RHO0, c=C0):
    """Transfer matrix (A, B, C, D) with [p_d, U_d] = [[A, B], [C, D]] [p_t, U_t]."""
    omega = np.asarray(omega, dtype=float)
    one = np.ones_like(omega, dtype=complex)
    if area_m2 is None:
        return one, 0.0 * one, 1j * omega * volume_m3 / (rho * c**2), one
    length = volume_m3 / area_m2
    kl = omega * length / c
    impedance = rho * c / area_m2
    return np.cos(kl) * one, 1j * impedance * np.sin(kl), 1j * np.sin(kl) / impedance, np.cos(kl) * one


def driver_operating_point(driver, frequencies, z_horn, throat_area, v_g, *, rho=RHO0, c=C0):
    """Solve the motor, chamber and horn network. Inputs are assumed validated.

    z_horn is the complex specific acoustic throat impedance (Pa s/m).
    """
    frequencies = np.asarray(frequencies, dtype=float)
    omega = 2.0 * np.pi * frequencies
    sd = driver.sd_m2
    moving_mass = driver.coupled_moving_mass_kg
    rms = driver.rms_kg_per_s if driver.rms_kg_per_s is not None else 0.0

    z_e = driver.re_ohm + 1j * omega * driver.le_h
    z_mech_driver = rms + 1j * omega * moving_mass + 1.0 / (1j * omega * driver.cms_m_per_n)
    if driver.rear_chamber_volume_m3 is not None:
        z_mech_driver = z_mech_driver + sd**2 * rho * c**2 / (1j * omega * driver.rear_chamber_volume_m3)

    volume, area, source = resolve_front_chamber(driver)
    ideal = volume == 0.0
    if ideal:
        # Previous ideal area transformer, kept verbatim for exact reproduction.
        z_mech_load = z_horn * (sd ** 2) / throat_area
    else:
        z_acoustic = z_horn / throat_area
        a, b, cc, d = front_chamber_matrix(omega, volume, area, rho=rho, c=c)
        throat_ratio = cc * z_acoustic + d  # U_d / U_t
        z_mech_load = sd**2 * (a * z_acoustic + b) / throat_ratio

    z_mech_total = z_mech_driver + z_mech_load
    # Solve the two motor equations without dividing by mechanical impedance.
    # The undamped, unloaded resonance has zero current and finite velocity.
    denominator = z_e * z_mech_total + driver.bl_tm**2
    current = v_g * z_mech_total / denominator
    velocity = driver.bl_tm * v_g / denominator
    with np.errstate(divide="ignore", invalid="ignore"):
        electrical_impedance = np.divide(v_g, current)
    electrical_impedance = np.where(current == 0, complex(np.inf, 0), electrical_impedance)

    if ideal:
        p_throat = z_horn * velocity * (sd / throat_area)
        u_throat = velocity * sd
        horn_power = np.abs(velocity)**2 * np.real(z_mech_load)
    else:
        u_throat = velocity * sd / throat_ratio
        p_throat = z_acoustic * u_throat
        horn_power = np.real(z_acoustic) * np.abs(u_throat)**2

    return {"drive_voltage_rms": v_g, "throat_pressure": p_throat, "velocity_rms": velocity,
            "throat_volume_velocity": u_throat,
            "current_rms": current, "electrical_impedance": electrical_impedance,
            "displacement_peak_m": np.sqrt(2) * np.abs(velocity) / omega,
            "input_power_w": np.real(v_g * np.conjugate(current)),
            "copper_power_w": np.abs(current)**2 * driver.re_ohm,
            "mechanical_loss_w": np.abs(velocity)**2 * rms,
            "horn_power_w": horn_power,
            "front_chamber_volume_m3": volume,
            "front_chamber_area_m2": area,
            "front_chamber_source": source,
            "rear_chamber_volume_m3": driver.rear_chamber_volume_m3}


def interface_summary(driver, throat_area):
    """Chamber assumption actually used, for rankings and reports."""
    volume, area, source = resolve_front_chamber(driver)
    corner = C0 * throat_area / (2 * math.pi * volume) if volume > 0 else None
    return {"front_chamber_volume_cm3": volume * 1e6, "front_chamber_source": source,
            "front_chamber_area_cm2": None if area is None else area * 1e4,
            "front_chamber_corner_hz": corner,
            "rear_chamber_volume_cm3": None if driver.rear_chamber_volume_m3 is None else driver.rear_chamber_volume_m3 * 1e6,
            "compression_ratio": driver.sd_m2 / throat_area}


def interface_evidence_gaps(driver, throat_area, f_high_hz):
    """Evidence gaps for the chamber assumptions. These flag, they do not reject.

    * front_chamber_volume_estimated: cone geometric estimate in use.
    * front_chamber_volume_unknown: no chamber data and no estimate (compression
      or unknown driver type); the ideal transformer is optimistic when Sd > St.
    * diaphragm_smaller_than_throat: Sd < St, the interface expands instead of
      compressing, which the plane-wave transformer does not describe.
    * front_chamber_not_lumped_in_band: the target band reaches the first radial
      mode of a cavity of the diaphragm's radius (k a = 3.83), above which a
      single uniform chamber pressure is not a valid assumption.
    """
    volume, _, source = resolve_front_chamber(driver)
    gaps = []
    if source == "estimated":
        gaps.append("front_chamber_volume_estimated")
    elif source == "unknown" and driver.sd_m2 > throat_area * (1 + 1e-12):
        gaps.append("front_chamber_volume_unknown")
    if driver.sd_m2 < throat_area * (1 - 1e-12):
        gaps.append("diaphragm_smaller_than_throat")
    radius = math.sqrt(driver.sd_m2 / math.pi)
    if volume > 0 and f_high_hz > FIRST_RADIAL_CAVITY_ROOT * C0 / (2 * math.pi * radius):
        gaps.append("front_chamber_not_lumped_in_band")
    return gaps
