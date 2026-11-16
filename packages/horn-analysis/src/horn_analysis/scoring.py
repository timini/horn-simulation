"""Target specification and the preference score used to rank feasible designs.

Hard constraints (ripple limit, minimum output, size caps, driver limits)
decide feasibility in ``evaluation.evaluate_response``. Infeasible designs
score 0. Feasible designs are ordered by a weighted sum of four terms, each
mapped linearly onto [0, 1] with fixed reference scales:

* sensitivity: ``clip((L_mean - 60 dB) / 60 dB, 0, 1)``, where ``L_mean`` is
  the log-frequency mean observer level over the target band;
* ripple: ``clip(1 - ripple / 12 dB, 0, 1)``, peak-to-peak level over the
  target band;
* extension: ``clip(0.5 + E / 2 octaves, 0, 1)``, where
  ``E = log2(f_low / f3_low) + log2(f3_high / f_high)`` and the f3 points are
  3 dB below ``L_mean``, searched outward from the band edges over the full
  simulated sweep (negative E means the response is already 3 dB down
  inside the band);
* size: ``clip(0.5 - log2(V / V_ref) / 8, 0, 1)`` with ``V = pi r_mouth^2 L``
  and ``V_ref = pi (c / (2 pi f_low))^2 * c / (4 f_low)``, the cylinder of a
  ka = 1 mouth and a quarter-wave length at f_low.

The weights are normalised to sum to 1. Because every term is linear in its
quantity, the weights set fixed exchange rates. With the defaults
(sensitivity 0.625, ripple 0.25, extension 0.125, size 0):

* 1 dB less ripple is worth 2 dB of sensitivity;
* 1 octave more extension is worth 6 dB of sensitivity;
* size does not affect the order.

So a design 15 dB more sensitive is not beaten by one that is 1 dB (or even
3 dB) flatter. ``exchange_rates`` reports the rates for any weights.
"""
from dataclasses import dataclass, field
import math
from typing import Optional

import numpy as np

SENSITIVITY_FLOOR_DB = 60.0
SENSITIVITY_SPAN_DB = 60.0
RIPPLE_SPAN_DB = 12.0
EXTENSION_SPAN_OCTAVES = 2.0
SIZE_SPAN_DOUBLINGS = 8.0
SPEED_OF_SOUND_M_S = 343.0

DEFAULT_SCORE_WEIGHTS = {"sensitivity": 0.625, "ripple": 0.25, "extension": 0.125, "size": 0.0}

SCORE_FORMULA = (
    "score = w_sensitivity*clip((L_mean-60)/60,0,1) + w_ripple*clip(1-ripple/12,0,1)"
    " + w_extension*clip(0.5+E/2,0,1) + w_size*clip(0.5-log2(V/V_ref)/8,0,1); weights sum to 1;"
    " infeasible designs score 0. L_mean: log-frequency mean level over the target band (dB);"
    " ripple: its peak-to-peak range (dB); E = log2(f_low/f3_low)+log2(f3_high/f_high) (octaves),"
    " f3 3 dB below L_mean searched outward from the band edges over the full sweep;"
    " V = pi*r_mouth^2*L; V_ref = pi*(c/(2*pi*f_low))^2*c/(4*f_low)."
)


def parse_score_weights(weights):
    """Return complete, normalised weights from None, a dict or 'key=value,...'."""
    if weights is None:
        return dict(DEFAULT_SCORE_WEIGHTS)
    if isinstance(weights, str):
        parsed = {}
        for item in filter(None, (part.strip() for part in weights.split(","))):
            key, sep, value = item.partition("=")
            if not sep:
                raise ValueError(f"Score weight '{item}' must be written key=value")
            parsed[key.strip()] = float(value)
        weights = parsed
    unknown = set(weights) - set(DEFAULT_SCORE_WEIGHTS)
    if unknown:
        raise ValueError(f"Unknown score weight(s): {sorted(unknown)}; use {sorted(DEFAULT_SCORE_WEIGHTS)}")
    merged = {**DEFAULT_SCORE_WEIGHTS, **{k: float(v) for k, v in weights.items()}}
    if any(not np.isfinite(v) or v < 0 for v in merged.values()):
        raise ValueError("Score weights must be finite and nonnegative")
    total = sum(merged.values())
    if total <= 0:
        raise ValueError("At least one score weight must be positive")
    return {k: v / total for k, v in merged.items()}


@dataclass
class TargetSpec:
    """User-specified target band, constraints and score weights."""
    f_low_hz: float
    f_high_hz: float
    max_length_m: Optional[float] = None
    max_mouth_radius_m: Optional[float] = None
    voltage_rms: float = 2.83
    observation_distance_m: float = 1.0
    max_ripple_db: float = 6.0
    max_compression_ratio: float = 10.0
    min_output_db: Optional[float] = None
    score_weights: Optional[dict] = field(default=None)
    # Separate stress level (W into nominal impedance) for excursion checks over
    # the full sweep; the SPL used for scoring stays at voltage_rms.
    test_power_w: Optional[float] = None

    def __post_init__(self):
        from horn_core.acoustics import validate_band
        validate_band(self.f_low_hz, self.f_high_hz)
        for name in ("voltage_rms", "observation_distance_m", "max_ripple_db", "max_compression_ratio", "max_length_m", "max_mouth_radius_m", "test_power_w"):
            value = getattr(self, name)
            if value is not None and (not np.isfinite(value) or value <= 0):
                raise ValueError(f"{name} must be finite and positive")
        if self.min_output_db is not None and not np.isfinite(self.min_output_db):
            raise ValueError("Minimum output must be finite")
        self.score_weights = parse_score_weights(self.score_weights)


def reference_volume_m3(f_low_hz):
    """Cylinder of a ka = 1 mouth and quarter-wave length at f_low."""
    radius = SPEED_OF_SOUND_M_S / (2 * math.pi * f_low_hz)
    return math.pi * radius**2 * SPEED_OF_SOUND_M_S / (4 * f_low_hz)


def score_terms(target, average_db, ripple_db, extension_octaves, mouth_radius=None, length=None):
    """Each preference term on [0, 1]; size is neutral (0.5) when unknown."""
    terms = {
        "sensitivity": float(np.clip((average_db - SENSITIVITY_FLOOR_DB) / SENSITIVITY_SPAN_DB, 0, 1)),
        "ripple": float(np.clip(1 - ripple_db / RIPPLE_SPAN_DB, 0, 1)),
        "extension": float(np.clip(0.5 + extension_octaves / EXTENSION_SPAN_OCTAVES, 0, 1)),
        "size": 0.5,
    }
    if mouth_radius is not None and length is not None and mouth_radius > 0 and length > 0:
        ratio = math.pi * mouth_radius**2 * length / reference_volume_m3(target.f_low_hz)
        terms["size"] = float(np.clip(0.5 - math.log2(ratio) / SIZE_SPAN_DOUBLINGS, 0, 1))
    return terms


def weighted_score(target, terms):
    return float(sum(target.score_weights[k] * terms[k] for k in DEFAULT_SCORE_WEIGHTS))


def exchange_rates(weights):
    """Score-equivalent dB of sensitivity per unit of each other quantity."""
    per_db = weights["sensitivity"] / SENSITIVITY_SPAN_DB
    if per_db == 0:
        return {"note": "sensitivity weight is zero; no dB exchange rate is defined"}
    return {
        "sensitivity_db_per_db_of_ripple": weights["ripple"] / RIPPLE_SPAN_DB / per_db,
        "sensitivity_db_per_octave_of_extension": weights["extension"] / EXTENSION_SPAN_OCTAVES / per_db,
        "sensitivity_db_per_doubling_of_horn_volume": weights["size"] / SIZE_SPAN_DOUBLINGS / per_db,
    }


FRONT_CHAMBER_POLICY = (
    "Levels and scores include each driver's front chamber. A cone driver without chamber data gets a"
    " geometric estimate (front_chamber_estimated); a constraint that fails only with that estimate is"
    " listed in estimated_chamber_only_reasons and flagged, not rejected. A supplied volume, including 0,"
    " is used as given.")

_PIPELINE_DEFAULTS = {"max_ripple_db": 6.0, "max_compression_ratio": 10.0}
_BOUND_KEYS = {
    "throat_radius_m": ("throat_radius", None, None, "throat_radius_range"),
    "mouth_radius_m": ("mouth_radius", "min_mouth_radius", "max_mouth_radius", "mouth_radius_range"),
    "length_m": ("length", "min_length", "max_length", "length_range"),
}


def describe_scoring(target, design=None, user_bounds=None):
    """Record the formula, weights and the origin of every bound for outputs."""
    user_bounds = {k: v for k, v in (user_bounds or {}).items() if v is not None}
    constraints = {}
    for name in ("max_ripple_db", "max_compression_ratio"):
        value = getattr(target, name)
        constraints[name] = {"value": value, "source": "pipeline_default" if value == _PIPELINE_DEFAULTS[name] else "user"}
    constraints["min_output_db"] = {"value": target.min_output_db,
                                    "source": "user" if target.min_output_db is not None else "not_set"}
    bounds = {}
    for label, (fixed, low, high, design_key) in _BOUND_KEYS.items():
        searched = (design or {}).get(design_key) or None
        if fixed in user_bounds:
            source = "user_fixed"
        elif (low and low in user_bounds) or (high and high in user_bounds):
            source = "heuristic_with_user_limit"
        elif label == "throat_radius_m":
            source = "driver_prescreen"
        else:
            source = "heuristic"
        bounds[label] = {"searched_range": searched, "source": source,
                         "user_values": {k: user_bounds[k] for k in (fixed, low, high) if k and k in user_bounds}}
    return {
        "formula": SCORE_FORMULA,
        "weights": dict(target.score_weights),
        "weights_source": "pipeline_default" if target.score_weights == parse_score_weights(None) else "user",
        "exchange_rates": exchange_rates(target.score_weights),
        "target_band_hz": [target.f_low_hz, target.f_high_hz],
        "constraints": constraints,
        "geometry_bounds": bounds,
        "size_compromise": (design or {}).get("size_compromise"),
        "front_chamber_policy": FRONT_CHAMBER_POLICY,
    }
