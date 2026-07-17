# -*- coding: utf-8 -*-
"""
Standard atmospheric climatology profiles and derived column quantities.
"""

import numpy as np

from astropy import units as u
from astropy import constants as c

from .constants import MOD_DIR, MASS_WATER, RHO_WATER
from .thermo import interp_by_pressure


class Climatology:
    names = (
            "midlatitude_summer",
            "midlatitude_winter",
            "subarctic_summer",
            "subarctic_winter",
            "tropical",
            "us_standard",
    )

    @u.quantity_input
    def __init__(
                self,
                name: str="midlatitude_winter",
                pressure_base: u.Quantity["pressure"]|None=None,  # noqa: F821
        ):
        """
        Standard atmospheric vertical profiles for the US and North America as
        reported in Anderson et al. (1986) "AFGL Atmospheric Constituent
        Profiles (0-120km)", AFGL-TR-86-0110.

        Args:
            name: str
                Name of the climatology, e.g., "midlatitude_winter".
            pressure_base:
                Mask the profiles for pressure values less than the given
                value. Profiles range from 1018 to 3.6e-5 hPa.
        """
        if name not in self.names:
            raise ValueError(f"Invalid name: {name}")
        self.name = name
        climatology_path = MOD_DIR / "climatology"
        data = np.loadtxt(climatology_path / f"{name}.dat")
        gas_minor = np.loadtxt(climatology_path / "gas_minor.dat")
        gas_trace = np.loadtxt(climatology_path / "gas_trace.dat")
        # Clip the pressure levels to a value if given.
        pressure = data[:, 1] * u.mbar

        def clip(arr):
            return interp_by_pressure(arr, pressure, pressure_base)
        # The slices will create views but multiplying by units will create copies.
        # Using the "<<" syntax to add units "in place" also creates copies in
        # this case, so just use the "*" for better readibility.
        self.altitude    = clip(data[:, 0] * u.km)
        self.pressure    = clip(data[:, 1] * u.mbar)
        self.density     = clip(data[:, 2] * u.cm**-3)
        self.temperature = clip(data[:, 3] * u.K)
        # Technically volumetric mixing ratio has units of "mol/mol"
        # The units are in parts-per-million, normalize to unity as the "vmr"
        # parameter expects in AM.
        from_ppm = 1e-6 * u.dimensionless_unscaled
        self.mixing_ratio = {
                # Major species
                "h2o":   clip(data[:, 4] * from_ppm),
                "co2":   clip(data[:, 5] * from_ppm),
                "o3":    clip(data[:, 6] * from_ppm),
                "n2o":   clip(data[:, 7] * from_ppm),
                "co":    clip(data[:, 8] * from_ppm),
                "ch4":   clip(data[:, 9] * from_ppm),
                "o2":    clip(data[:,10] * from_ppm),
                # Minor species
                "no":    clip(gas_minor[:, 0] * from_ppm),
                "so2":   clip(gas_minor[:, 1] * from_ppm),
                "no2":   clip(gas_minor[:, 2] * from_ppm),
                "nh3":   clip(gas_minor[:, 3] * from_ppm),
                "hno3":  clip(gas_minor[:, 4] * from_ppm),
                "oh":    clip(gas_minor[:, 5] * from_ppm),
                "hf":    clip(gas_minor[:, 6] * from_ppm),
                "hcl":   clip(gas_minor[:, 7] * from_ppm),
                "hbr":   clip(gas_minor[:, 8] * from_ppm),
                "clo":   clip(gas_minor[:,10] * from_ppm),
                "ocs":   clip(gas_minor[:,11] * from_ppm),
                "h2co":  clip(gas_minor[:,12] * from_ppm),
                "hocl":  clip(gas_minor[:,13] * from_ppm),
                "hcn":   clip(gas_minor[:,15] * from_ppm),
                "h2o2":  clip(gas_minor[:,17] * from_ppm),
                # Trace species
                "h2s":   clip(gas_trace[:, 2] * from_ppm),
        }

    @classmethod
    def midlat_from_datetime(cls, dt):
        if dt.month in range(4, 10):
            return cls("midlatitude_summer")
        else:
            return cls("midlatitude_winter")

    @property
    def pwv(self):
        column_density = self.column_density("h2o")
        return (column_density * MASS_WATER / RHO_WATER).to("mm")

    def column_density(self, specie):
        """
        Use the Ideal Gas Law to calculate air volume density and then use a
        trapezoidal sum to integrate the volume density as a function of
        altitude.
        """
        # c.k_B: Boltzmann's constant
        air_density = self.pressure / (c.k_B * self.temperature)
        mixing_ratio = self.mixing_ratio[specie]
        try:
            column_density = np.trapezoid(mixing_ratio * air_density, x=self.altitude)
        except AttributeError:
            column_density = np.trapz(mixing_ratio * air_density, x=self.altitude)
        return column_density.to("cm-2")

    def dobson_unit(self, specie):
        """
        The Dobson unit (DU) is defined as the thickness (in units of 10 μm) of
        that layer of pure gas which would be formed by the total column amount
        at standard conditions for temperature and pressure (STP).
        """
        return self.column_density(specie).to("m-2").value / 2.69e20


# Loading all six climatologies reads 18 data files, so the CLIMATOLOGIES
# dict is built lazily on first attribute access (PEP 562) and cached into
# the module globals. The full dict is built at once so that lookup, length,
# and iteration semantics match the previous eager definition exactly.
def __getattr__(name):
    if name == "CLIMATOLOGIES":
        globals()["CLIMATOLOGIES"] = climatologies = {
                n: Climatology(n) for n in Climatology.names}
        return climatologies
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | {"CLIMATOLOGIES"})
