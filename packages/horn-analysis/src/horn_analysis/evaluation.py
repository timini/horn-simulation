"""Common target-band acceptance for analytical and FEM candidates.

Preference scores never override hard constraints. Feasibility in an ideal
model does not certify a physical design; evidence gaps remain explicit.
"""
import numpy as np
from horn_core.acoustics import validate_response, inlet_area_from_frame, baffled_piston_on_axis, pressure_level
from horn_analysis.transfer_function import compute_driver_response, scale_solver_spl


def coupled_output(frame, driver, target, legacy_radius=None):
    from horn_core.duct import DEFAULT_AIR
    for name,value in [("air_c_m_s",DEFAULT_AIR.c),("air_rho_kg_m3",DEFAULT_AIR.rho)]:
        if name in frame and not np.allclose(frame[name],value,rtol=1e-12):
            raise ValueError("Driver/observer coupling currently requires the default air properties")
    f = validate_response(frame.frequency, frame.spl, frame.z_real, frame.z_imag)
    area = inlet_area_from_frame(frame, legacy_radius)
    p = compute_driver_response(driver, f, frame.z_real.to_numpy(), frame.z_imag.to_numpy(), area, target.voltage_rms)
    if "schema_version" not in frame:
        return scale_solver_spl(frame.spl.to_numpy(), p), area, "legacy_mouth_pressure"
    if not (frame.schema_version == 2).all() or not (frame.phasor_convention == "exp(+iwt)_rms").all():
        raise ValueError("Unsupported solver phase convention/schema")
    if not (frame.bc_mode == "dirichlet").all():
        raise ValueError("Driver coupling requires a unit-pressure Dirichlet transfer run")
    mouth_area = frame.mouth_area_m2.to_numpy()
    if np.any(mouth_area <= 0) or not np.isfinite(mouth_area).all() or not np.allclose(mouth_area, mouth_area[0], rtol=1e-6):
        raise ValueError("Inconsistent mouth area")
    u = frame.mouth_u_real.to_numpy() + 1j*frame.mouth_u_imag.to_numpy()
    validate_response(f, u)
    if 'radiation_model' in frame and (frame.radiation_model=='modal_baffled').any():
        if not (frame.radiation_model=='modal_baffled').all() or 'modal_mode_count' not in frame or not (frame.modal_mode_count==16).all():
            raise ValueError('Inconsistent modal aperture contract')
        from horn_core.modal_radiation import modal_baffled_on_axis
        columns=[f'modal_v_{n}_{part}' for n in range(16) for part in ('real','imag')]
        if any(name not in frame for name in columns):
            raise ValueError('Missing modal aperture velocity coefficients')
        velocity=np.column_stack([frame[f'modal_v_{n}_real'].to_numpy()+1j*frame[f'modal_v_{n}_imag'].to_numpy() for n in range(16)])
        if not np.isfinite(velocity).all() or not np.allclose(u,mouth_area*velocity[:,0],rtol=1e-8,atol=1e-15):
            raise ValueError('Inconsistent modal aperture volume velocity')
        observer=modal_baffled_on_axis(f,velocity*p[:,None],np.sqrt(mouth_area[0]/np.pi),target.observation_distance_m)
        return pressure_level(observer),area,'modal_baffled_on_axis'
    observer = baffled_piston_on_axis(f, u*p, float(mouth_area[0]), target.observation_distance_m)
    return pressure_level(observer), area, "uniform_baffled_piston_on_axis"


def evaluate_response(frequencies, levels, target, driver=None, throat_area=None, length=None, mouth_radius=None, operating_point=None):
    f = validate_response(frequencies, levels)
    levels = np.asarray(levels)
    if f[0] > target.f_low_hz or f[-1] < target.f_high_hz:
        raise ValueError("Incomplete response: target band is not bracketed")
    # Include every original sample and both band edges; no smoothing of notches.
    selected = f[(f > target.f_low_hz) & (f < target.f_high_hz)]
    grid = np.r_[target.f_low_hz, selected, target.f_high_hz]
    x = np.log(grid)
    y = np.interp(x, np.log(f), levels)
    ripple = float(np.ptp(y))
    integrate = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    average = float(integrate(y, x)/(x[-1]-x[0]))
    # f3 is referenced to the same target-band mean as the ripple and mean
    # level columns, and searched over the full sweep so extension counts.
    log_f = np.log(f)
    f3_threshold = average - 3.0
    f3_low, low_bound = _f3_edge(log_f, levels, x[0], x[-1], f3_threshold, -1)
    f3_high, high_bound = _f3_edge(log_f, levels, x[-1], x[0], f3_threshold, 1)
    extension = float(np.log2(target.f_low_hz/f3_low) + np.log2(f3_high/target.f_high_hz))
    reasons, missing = [], []
    if ripple > target.max_ripple_db + 1e-9:
        reasons.append("target_band_ripple_exceeded")
    if target.min_output_db is not None and y.min() < target.min_output_db:
        reasons.append("minimum_output_not_met")
    if length is not None and target.max_length_m is not None and length > target.max_length_m:
        reasons.append("maximum_length_exceeded")
    if mouth_radius is not None and target.max_mouth_radius_m is not None and mouth_radius > target.max_mouth_radius_m:
        reasons.append("maximum_mouth_radius_exceeded")
    if driver is not None:
        if throat_area is None or throat_area <= 0:
            raise ValueError("Driver assessment requires positive throat area")
        if mouth_radius is not None and driver.sd_m2 > np.pi * mouth_radius**2 * (1 + 1e-12):
            reasons.append("driver_larger_than_mouth")
        if driver.sd_m2 / throat_area > target.max_compression_ratio:
            reasons.append("compression_ratio_exceeded")
        if driver.usable_f_low_hz is None or driver.usable_f_high_hz is None:
            missing.append("driver_usable_band_unknown")
        if ((driver.usable_f_low_hz is not None and driver.usable_f_low_hz > target.f_low_hz)
                or (driver.usable_f_high_hz is not None and driver.usable_f_high_hz < target.f_high_hz)):
            reasons.append("outside_driver_usable_band")
        if driver.rms_kg_per_s is None:
            missing.append("mechanical_damping_unknown")
        if not driver.parameter_source:
            missing.append("driver_parameter_provenance_missing")
        if not driver.interface_model:
            missing.append("driver_interface_unverified")
        if driver.mmd_kg is None or driver.rear_load_mass_kg is None:
            missing.append("moving_mass_and_rear_load_not_separated")
        if driver.xmax_m is None or driver.power_w is None:
            missing.append("maximum_output_limits_unknown")
    operating_limits = {}
    if driver is not None and operating_point is not None:
        for quantity, limit, reason in [("displacement_peak_m",driver.xmax_m,"driver_excursion_exceeded"),
                                         ("input_power_w",driver.power_w,"driver_nominal_power_exceeded")]:
            values = np.asarray(operating_point[quantity])
            validate_response(f, values)
            peak = float(np.max(np.interp(x,np.log(f),values)))
            operating_limits["max_"+quantity] = peak
            if limit is not None and peak > limit:
                reasons.append(reason)
        operating_limits["operating_limit_scope"] = "Linear model at requested voltage; not a validated maximum-output or thermal prediction"
    from horn_analysis.scoring import score_terms, weighted_score
    terms = score_terms(target, average, ripple, extension, mouth_radius, length)
    score = weighted_score(target, terms)
    return {
        **operating_limits,
        "passband_ripple_db": ripple,
        "f3_low_hz": f3_low, "f3_high_hz": f3_high,
        "f3_low_is_bound": low_bound, "f3_high_is_bound": high_bound,
        "f3_reference": "target_band_mean_minus_3db",
        "extension_octaves": extension,
        "score_terms": terms,
        "avg_sensitivity_db": average,  # retained JSON key; reports label actual quantity
        "composite_score": 0.0 if reasons else score,
        "model_feasible": not reasons,
        "eligibility_status": "infeasible" if reasons else ("insufficient_evidence" if missing else "model_feasible"),
        "rejection_reasons": reasons, "evidence_gaps": missing,
        "validation_status": "experimental_prediction",
        "target_min_level_db": float(y.min()),
        "target_max_level_db": float(y.max()),
    }


def _f3_edge(log_f, levels, log_edge, log_other_edge, threshold, direction):
    """First crossing of ``threshold`` measured from one target-band edge.

    If the level at the edge is at or above the threshold, walk outward
    (``direction`` -1 below the band, +1 above) over the full sweep to the
    first sample below it. Otherwise walk inward towards the other band edge
    to the first sample at or above it. Crossings are interpolated linearly
    against log frequency. Returns (f3_hz, is_sweep_bound).
    """
    level_at_edge = float(np.interp(log_edge, log_f, levels))
    outward = level_at_edge >= threshold
    step = direction if outward else -direction
    if step > 0:
        mask = (log_f > log_edge) if outward else (log_f > log_edge) & (log_f < log_other_edge)
        path_x, path_y = log_f[mask], levels[mask]
    else:
        mask = (log_f < log_edge) if outward else (log_f < log_edge) & (log_f > log_other_edge)
        path_x, path_y = log_f[mask][::-1], levels[mask][::-1]
    xs = np.r_[log_edge, path_x]
    ys = np.r_[level_at_edge, path_y]
    for i in range(1, len(xs)):
        crossed = ys[i] < threshold if outward else ys[i] >= threshold
        if crossed:
            fraction = (ys[i-1]-threshold)/(ys[i-1]-ys[i])
            return float(np.exp(xs[i-1] + fraction*(xs[i]-xs[i-1]))), False
    if outward:
        return float(np.exp(xs[-1])), True
    return float(np.exp(log_other_edge)), False


def radiation_domain_rejection(frequencies, mouth_radius, radiation_model, flange_width=0.):
    """Reject unsupported geometries without masking invalid inputs/solver errors."""
    if radiation_model=='modal_baffled' and np.max(2*np.pi*np.asarray(frequencies)*mouth_radius/343.)>30:
        return {"model_feasible":False,"simulation_eligible":False,
                "eligibility_status":"infeasible","composite_score":0.,
                "rejection_reasons":["radiation_model_out_of_domain"],
                "rejection_detail":"Modal aperture quadrature requires ka <= 30",
                "evidence_gaps":[],"validation_status":"outside_supported_model_domain",
                "extension_octaves":None,"passband_ripple_db":None,"avg_sensitivity_db":None}
    if radiation_model not in {"finite_flange", "unflanged", "unflanged_piston"}:
        return None
    from horn_core.duct import DEFAULT_AIR, RadiationDomainError, circular_pipe_radiation
    try:
        circular_pipe_radiation(2*np.pi*np.asarray(frequencies)/DEFAULT_AIR.c, mouth_radius, flange_width if radiation_model == "finite_flange" else 0.)
    except RadiationDomainError as error:
        return {"model_feasible": False, "simulation_eligible": False,
                "eligibility_status": "infeasible", "composite_score": 0.,
                "rejection_reasons": ["radiation_model_out_of_domain"],
                "rejection_detail": str(error), "evidence_gaps": [],
                "validation_status": "outside_supported_model_domain",
                "extension_octaves": None, "passband_ripple_db": None,
                "avg_sensitivity_db": None}
    return None
