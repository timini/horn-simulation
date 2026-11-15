"""Neumann inlet velocity must use the same driver/chamber network as the analysis coupling."""
import numpy as np
import pytest

pytest.importorskip("dolfinx")

from horn_core.parameters import DriverParameters
from horn_solver.solver import _compute_neumann_velocity


@pytest.mark.parametrize("volume", [0., 2e-4])
def test_inlet_velocity_is_throat_volume_velocity_over_area(volume):
    driver = DriverParameters("d", "T", "M", 110., 5.3, 9., .013, .0085, .0001, qms=2.7, qes=.38,
                              driver_type="cone", front_chamber_volume_m3=volume)
    area, rho_c = .0033, 343. * 1.225
    w = 2 * np.pi * 1500.
    zm = driver.rms_kg_per_s + 1j * w * driver.mms_kg + 1 / (1j * w * driver.cms_m_per_n)
    z_acoustic = rho_c / area
    shunt = 1 + 1j * w * volume / (1.225 * 343.**2) * z_acoustic  # U_d / U_t
    load = driver.sd_m2**2 * z_acoustic / shunt
    velocity = driver.bl_tm * 2.83 / ((driver.re_ohm + 1j * w * driver.le_h) * (zm + load) + driver.bl_tm**2)
    expected = velocity * driver.sd_m2 / shunt / area
    assert _compute_neumann_velocity(1500., driver, area, rho_c, 0., 2.83) == pytest.approx(expected, rel=1e-10)
