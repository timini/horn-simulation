"""Infinite-baffle off-axis pressure from modal mouth velocities."""
import numpy as np
import pytest
from scipy.special import j1

from horn_core.modal_radiation import (modal_baffled_directivity_db, modal_baffled_far_field,
                                       modal_baffled_on_axis, radial_roots)


def test_uniform_mode_alone_gives_the_piston_pattern():
    radius = 0.1
    frequencies = np.array([200., 1000., 3000., 8000.])
    angles = np.linspace(0., 90., 91)
    velocities = np.zeros((len(frequencies), 16), dtype=complex)
    velocities[:, 0] = 0.3 - 0.2j
    pressure = modal_baffled_far_field(frequencies, velocities, radius, angles)
    k = 2 * np.pi * frequencies / 343.
    x = k[:, None] * radius * np.sin(np.radians(angles))[None, :]
    with np.errstate(invalid="ignore", divide="ignore"):
        piston = np.where(x > 0, 2 * j1(x) / x, 1.)
    ratio = pressure / pressure[:, :1]
    assert np.max(np.abs(ratio - piston)) < 1e-6
    levels = modal_baffled_directivity_db(frequencies, velocities, radius, angles[1:])
    expected = 20 * np.log10(np.maximum(np.abs(piston[:, 1:]), 1e-15))
    assert np.max(np.abs(levels - expected)) < 1e-6


def test_far_field_matches_the_on_axis_rayleigh_integral_at_long_range():
    rng = np.random.default_rng(3)
    radius, distance = 0.08, 400.
    frequencies = np.array([500., 2000.])
    velocities = rng.normal(size=(2, 16)) + 1j * rng.normal(size=(2, 16))
    k = 2 * np.pi * frequencies / 343.
    near = modal_baffled_on_axis(frequencies, velocities, radius, distance, order=512)
    far = modal_baffled_far_field(frequencies, velocities, radius, [0.])[:, 0]
    assert np.allclose(near * distance * np.exp(1j * k * distance), far, rtol=1e-3)


def test_higher_modes_radiate_nothing_on_axis_but_do_off_axis():
    velocities = np.zeros((1, 4), dtype=complex)
    velocities[0, 2] = 1.
    pressure = modal_baffled_far_field([4000.], velocities, 0.15, [0., 45., 90.])
    assert abs(pressure[0, 0]) < 1e-12
    assert abs(pressure[0, 2]) > 1e-3
    assert len(radial_roots(4)) == 4


def test_angles_outside_the_front_half_space_are_rejected():
    with pytest.raises(ValueError, match="0 and 90"):
        modal_baffled_far_field([1000.], np.ones((1, 2)), 0.1, [120.])
