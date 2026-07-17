# -*- coding: utf-8 -*-
"""
Wrapper for Scott Paine's AM atmospheric radiative transfer code.

The ``Model`` class provides an interface for configuring and running AM for
specified model atmospheric conditions. A model can be initialized from arrays
assigned units from `astropy.units` or constructed from a standard climatology.

Typical usage example:
    .. code-block:: python

        m = Model.from_climatology("midlatitude_winter")
        m.troposphere_h2o_scaling = 0.8
        df = m.run()
        print(df.attrs["warnings?"])
"""

# TODO
# - Perform interpolations/extrapolation on vertical profiles
# - Configure the water vapor mixing ratio given a total PWV and a scale
#   height using an exponential function.
# - Convert all units from Astropy to Metpy.

__all__ = [
    "AmExecutable",
    "AM_PARALLEL",
    "AM_SERIAL",
    "BOTH_AM_CALLABLE",
    "NO_AM_CALLABLE",
    "ENV",
    "CACHE_DIR",
    "MOD_DIR",
    "MASS_DRY_AIR",
    "MASS_DRY_AIR_MOL",
    "MASS_WATER",
    "RHO_WATER",
    "STD_TEMPERATURE",
    "STD_PRESSURE",
    "STD_LAPSE_RATE",
    "R_DRY_AIR",
    "mixing_ratio_from_relative_humidity",
    "precipitable_water",
    "altitude_from_pressure",
    "interp_by_pressure",
    "Climatology",
    "CLIMATOLOGIES",
    "Model",
]

from .driver import (
    AmExecutable,
    AM_PARALLEL,
    AM_SERIAL,
    NO_AM_CALLABLE,
    BOTH_AM_CALLABLE,
    ENV,
    CACHE_DIR,
)
from .constants import (
    MOD_DIR,
    MASS_DRY_AIR,
    MASS_DRY_AIR_MOL,
    MASS_WATER,
    RHO_WATER,
    STD_TEMPERATURE,
    STD_PRESSURE,
    STD_LAPSE_RATE,
    R_DRY_AIR,
)
from .thermo import (
    mixing_ratio_from_relative_humidity,
    precipitable_water,
    altitude_from_pressure,
    interp_by_pressure,
)
from .climatology import Climatology, CLIMATOLOGIES
from .model import Model
