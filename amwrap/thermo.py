# -*- coding: utf-8 -*-
"""
Thermodynamic and vertical-profile helper functions.
"""

import numpy as np

from astropy import units as u
from astropy import constants as c

from .constants import (
    MASS_DRY_AIR,
    MASS_WATER,
    STD_TEMPERATURE,
    STD_PRESSURE,
    STD_LAPSE_RATE,
    R_DRY_AIR,
)


@u.quantity_input
def mixing_ratio_from_relative_humidity(
        pressure: u.Quantity["pressure"],  # noqa: F821
        temperature: u.Quantity[u.deg_C] | u.Quantity["temperature"],  # noqa: F821
        relative_humidity: u.Quantity["dimensionless"],  # noqa: F821
    ):
    from pint import Quantity
    from metpy import calc
    p  = pressure.to("hPa").value * Quantity("hPa")
    t  = temperature.to("deg_C", equivalencies=u.temperature()).value * Quantity("degC")
    rh = relative_humidity.to("").value * Quantity("dimensionless")
    mr = calc.mixing_ratio_from_relative_humidity(p, t, rh)
    # Returned value is mass mixing ratio, convert to volumetric mixing ratio using masses
    return mr.m * MASS_DRY_AIR / MASS_WATER * u.dimensionless_unscaled


@u.quantity_input
def precipitable_water(
        pressure: u.Quantity["pressure"],  # noqa: F821
        temperature: u.Quantity[u.deg_C] | u.Quantity["temperature"],  # noqa: F821
        relative_humidity: u.Quantity["dimensionless"],  # noqa: F821
    ):
    from pint import Quantity
    from metpy import calc
    dewpoint = calc.dewpoint_from_relative_humidity(
            temperature.to("deg_C", equivalencies=u.temperature()).value * Quantity("degC"),
            relative_humidity.to("").value * Quantity("dimensionless"),
    )
    pwv = calc.precipitable_water(
            pressure.to("hPa").value * Quantity("hPa"),
            dewpoint,
    )
    return pwv.m_as("mm") * u.mm


@u.quantity_input
def altitude_from_pressure(pressure: u.Quantity["pressure"]):  # noqa: F821
    """
    Convert pressure to height using the U.S. Standard Atmosphere (NOAA 1976).
    Implementation taken from `metpy.calc.pressure_height_std` itself from
    the formula outlined in Hobbs & Wallace (1977) pg. 60-61.
    """
    prefactor = STD_TEMPERATURE / STD_LAPSE_RATE
    exponent = R_DRY_AIR * STD_LAPSE_RATE / c.g0
    altitude = prefactor * (1 - (pressure / STD_PRESSURE).to("")**exponent)
    return altitude.to("km")


@u.quantity_input
def interp_by_pressure(
            values,
            pressure: u.Quantity["pressure"]|None=None,  # noqa: F821
            pressure_base: u.Quantity["pressure"]|None=None,  # noqa: F821
    ):
    if pressure is not None:
        if values.shape != pressure.shape:
            raise ValueError(f"Invalid shapes: {values.shape=} != {pressure.shape=}")
    if pressure_base is not None:
        if pressure is None:
            raise ValueError("Pressure must be assigned when clipping to base.")
        if pressure_base > pressure.max() or pressure_base < pressure.min():
            raise ValueError(f"Pressure base outside bounds: {pressure_base=}")
    if pressure_base is None:
        return values
    else:
        ix_base = (pressure[pressure > pressure_base]).argmin()
        v_base = np.interp(pressure_base, pressure[::-1], values[::-1])
        _values = values.copy()
        _values[ix_base] = v_base
        return _values[ix_base:]
