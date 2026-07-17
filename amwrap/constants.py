# -*- coding: utf-8 -*-
"""
Physical constants and package paths shared across amwrap modules.
"""

import os
from pathlib import Path

from astropy import units as u
from astropy import constants as c


MOD_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
MASS_DRY_AIR =     28.96546 * u.u  # amu
MASS_DRY_AIR_MOL = MASS_DRY_AIR.value * 1e-3 * u.kg / u.mol
MASS_WATER =       18.01528 * u.u
RHO_WATER =         0.9998395 * u.g / u.cm**3
STD_TEMPERATURE = 288.0 * u.K
STD_PRESSURE =   1013.25 * u.hPa
STD_LAPSE_RATE =    6.5 * u.K / u.km
R_DRY_AIR = (c.R / MASS_DRY_AIR_MOL).to("J/(K kg)")
