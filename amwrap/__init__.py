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
# - Add options for including clouds through "lwp_abs_Rayleigh" and
#   "iwp_abs_Rayleigh" column types in layers.
# - Perform interpolations/extrapolation on vertical profiles
# - Have the PWV be configurable as a both a direct input in the mixing ratio as
#   well an exponential function with a scale height.
# - Add conversion between pressure and altitude (w/ Eotvos effect).

import os
import warnings
import subprocess
from io import BytesIO
from typing import Dict
from pathlib import Path
from datetime import datetime
from packaging import version

import numpy as np
import pandas as pd

from astropy import units as u
from astropy import constants as c


class AmExecutable:
    bin_dir = (Path(__file__).parent / "bin").absolute()

    def __init__(self, name):
        if name not in ("am", "am-serial"):
            raise ValueError(f"Invalid name: {name}")
        self.name = name
        try:
            # Test if the executable name is callable from the user's environment.
            result = subprocess.run([f"{name}", "-v"], capture_output=True)
            self.exec_name = name
            self.is_callable = True
            self.version = self._parse_version(result)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            # Fallback to the included executable.
            result = subprocess.run([f"{name}", "-v"], capture_output=True)
            self.exec_name = str(self.bin_dir / name)
            self.is_callable = True
            self.version = self.parse_version(result)

    @staticmethod
    def _parse_version(result):
        return version.parse(result.stdout.decode().split()[2])

# FIXME Use configuration system for "am"/"am-serial" executable names.
AM_PARALLEL = AmExecutable("am")
AM_SERIAL   = AmExecutable("am-serial")

NO_AM_CALLABLE = not AM_PARALLEL.is_callable and not AM_SERIAL.is_callable
BOTH_AM_CALLABLE = not NO_AM_CALLABLE
if NO_AM_CALLABLE:
    warnings.warn("No executable callable for AM.", UserWarning)
if BOTH_AM_CALLABLE and (AM_PARALLEL.version != AM_SERIAL.version):
    warnings.warn(
            "`am` and `am-serial` version mismatch: "
            f"{AM_PARALLEL.version} & {AM_SERIAL.version}",
            UserWarning,
    )

# Set environment variables for am
ENV = os.environ.copy()
CACHE_DIR = Path("/dev/shm/")
if CACHE_DIR.exists():
    ENV["AM_CACHE_PATH"] = CACHE_DIR

MOD_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
MASS_DRY_AIR = 28.9644  * u.u  # amu
MASS_WATER =   18.01528 * u.u
RHO_WATER =     0.9998395 * u.g / u.cm**3


def mixing_ratio_from_relative_humidity(pressure, temperature, relative_humidity):
    from metpy.units.pint import Quantity
    from metpy.calc import mixing_ratio_from_relative_humidity
    p  = pressure.to("hPa").value * Quantity("hPa")
    t  = temperature.to("deg_C", equivalencies=u.temperature()).value * Quantity("degC")
    rh = relative_humidity.value * Quantity("dimensionless")
    mr = mixing_ratio_from_relative_humidity(p, t, rh)
    # Returned value is mass mixing ratio, convert to volumetric mixing ratio using masses
    return mr.m * 1e6 * MASS_DRY_AIR / MASS_WATER * u.dimensionless_unscaled


class Climatology:
    names = (
            "midlatitude_summer",
            "midlatitude_winter",
            "subarctic_summer",
            "subarctic_winter",
            "tropical",
            "us_standard",
    )

    def __init__(self, name="midlatitude_winter"):
        if name not in self.names:
            raise ValueError(f"Invalid name: {name}")
        self.name = name
        climatology_path = MOD_DIR / "climatology"
        data = np.loadtxt(climatology_path / f"{name}.dat")
        gas_minor = np.loadtxt(climatology_path / "gas_minor.dat")
        gas_trace = np.loadtxt(climatology_path / "gas_trace.dat")
        # The slices will create views but multiplying by units will create copies.
        # Using the "<<" syntax to add units "in place" also creates copies in
        # this case, so just use the "*" for better readibility.
        self.altitude         = data[:, 0] * u.km
        self.pressure         = data[:, 1] * u.mbar
        self.density          = data[:, 2] * u.cm**-3
        self.temperature      = data[:, 3] * u.K
        # Technically volumetric mixing ratio has units of "mol/mol"
        # The units are in parts-per-million, normalize to unity as the "vmr"
        # parameter expects in AM.
        from_ppm = 1e-6 * u.dimensionless_unscaled
        self.mixing_ratio = {
                # Major species
                "h2o":   data[:, 4] * from_ppm,
                "co2":   data[:, 5] * from_ppm,
                "o3":    data[:, 6] * from_ppm,
                "n2o":   data[:, 7] * from_ppm,
                "co":    data[:, 8] * from_ppm,
                "ch4":   data[:, 9] * from_ppm,
                "o2":    data[:,10] * from_ppm,
                # Minor species
                "no":    gas_minor[:, 0] * from_ppm,
                "so2":   gas_minor[:, 1] * from_ppm,
                "no2":   gas_minor[:, 2] * from_ppm,
                "nh3":   gas_minor[:, 3] * from_ppm,
                "hno3":  gas_minor[:, 4] * from_ppm,
                "oh":    gas_minor[:, 5] * from_ppm,
                "hf":    gas_minor[:, 6] * from_ppm,
                "hcl":   gas_minor[:, 7] * from_ppm,
                "hbr":   gas_minor[:, 8] * from_ppm,
                "clo":   gas_minor[:,10] * from_ppm,
                "ocs":   gas_minor[:,11] * from_ppm,
                "h2co":  gas_minor[:,12] * from_ppm,
                "hocl":  gas_minor[:,13] * from_ppm,
                "hcn":   gas_minor[:,15] * from_ppm,
                "h2o2":  gas_minor[:,17] * from_ppm,
                # Trace species
                "h2s":   gas_trace[:, 2] * from_ppm,
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

CLIMATOLOGIES = {n: Climatology(n) for n in Climatology.names}


class Model:
    """
    AM model configuration.

    Attributes:
      background_temperature:
        Background brightness temperature. By default set to the temperature of
        the Cosmic Microwave Background.
      valid_output_descriptors:
        Map of valid AM output columns and default units.
    """
    background_temperature = 2.725 * u.K
    valid_output_descriptors = {
            # See Table B.3 of the AM manual.
            "frequency": "f GHz",
            "opacity": "tau neper",
            "transmittance": "tx none",
            "radiance": "I watt*cm-2*GHz-1*sr-1",
            "radiance difference": "I_diff watt*cm-2*GHz-1*sr-1",
            "brightness temperature": "Tb K",
            "Rayleigh-Jeans brightness temperature": "Trj K",
            "delay": "L mm",
            "absorption coefficient": "k auto",
    }
    valid_species = [
            "ch4", "12ch4", "13ch4", "12ch3d",
            "ch3cn", "12ch3_12c14n",
            "ch3oh", "12ch3_16oh",
            "co", "12c_16o", "13c_16o", "12c_18o", "12c_17o", "13c_18o", "13c_17o",
            "co2", "12c_16o2", "13c_16o2", "16o_12c_18o", "16o_12c_17o", "16o_13c_18o", "16o_13c_17o", "12c_18o2",
            "clo", "35cl_16o", "37cl_16o",
            "hbr", "h_79br", "h_81br",
            "hcn", "h_12c_14n", "h_13c_14n", "h_12c_15n",
            "h2co", "h2_12c_16o", "h2_13c_16o", "h2_12c_18o",
            "hcl", "h_35cl", "h_37cl",
            "hf", "h_19f",
            "hno3", "h_14n_16o3",
            "h2o", "h2o_lines", "h2_16o", "h2_18o", "h2_17o", "hd_16o", "hd_18o", "hd_17o",
            "h2o2", "h2_16o2",
            "ho2", "h_16o2",
            "hocl", "h_16o_35cl", "h_16o_37cl",
            "h2s", "h2_32s", "h2_34s", "h2_33s",
            "nh3", "14nh3", "15nh3",
            "n2o", "14n2_16o", "14n_15n_16o", "15n_14n_16o", "14n2_18o", "14n2_17o",
            "no", "14n_16o",
            "no2", "14n_16o2",
            "o", "16o",
            "o2", "o2_coupled", "16o2_coupled", "16o2", "16o_18o", "16o_17o",
            "o2_uncoupled", "16o2_uncoupled", "16o_18o_uncoupled", "16o_17o_uncoupled",
            "o3", "16o3", "16o2_18o", "16o_18o_16o", "16o2_17o", "16o_17o_16o",
            "ocs", "16o_12c_32s", "16o_12c_34s", "16o_13c_32s", "16o_12c_33s", "18o_12c_32s",
            "oh", "16oh", "18oh", "16od",
            "so2", "32s_16o2", "34s_16o2",
            "h2o_continuum",
            "h2o_air_continuum",
            "h2o_self_continuum",
            "n2n2",
            "n2air",
            "o2o2",
            "o2air",
            "lwp_abs_Rayleigh",
            "iwp_abs_Rayleigh",
    ]

    @u.quantity_input
    def __init__(self,
                pressure: u.Quantity["pressure"],  # noqa: F821
                temperature: u.Quantity["temperature"],  # noqa: F821
                mixing_ratio: Dict[str, u.Quantity["dimensionless"]|None]|None=None,  # noqa: F821
                zenith_angle: u.Quantity["angle"]=0*u.deg,  # noqa: F821
                freq_min: u.Quantity["frequency"]=18*u.GHz,  # noqa: F821
                freq_max: u.Quantity["frequency"]=26.5*u.GHz,  # noqa: F821
                freq_step: u.Quantity["frequency"]=10*u.MHz,  # noqa: F821
                troposphere_h2o_scaling=1.0,
                output_columns=(
                    "frequency",
                    "brightness temperature",
                    "opacity",
                    "delay",
                ),
                tolerance=1e-4,
                self_broadening_tolerance=0.003,
        ):
        """
        Initialize directly from array data of atmospheric properties. Mixing
        ratios are scaled to unity (not parts per million), and should range
        between 0 and 1.

        Args:
          pressure:
            Atmospheric pressure at the base of the layer.
          temperature:
            Atmospheric temperature at the base of the layer.
          mixing_ratio:
            Dictionary of volumetric mixing ratios indexed by a string for
            a specie name recognized by AM. If a specie is set to ``None``,
            then values are interpolated from the US Standard climatology.
          zenith_angle:
            Zenith angle of ray; 0 deg is towards zenith and 90 deg is towards
            the horizon.
          freq_min:
            Minimum frequency to compute output values over.
          freq_max:
            Maximum frequency to compute output values over.
          freq_step:
            Frequency step or "channel" size to compute output values over.
          troposphere_h2o_scaling:
            Scale the tropospheric water mixing ratio by this factor.
          output_columns:
            Values to calculate as part of the returned output. See available
            output columns in the `valid_output_descriptors` class attribute.
          tolerance:
            Numerical tolerance of line-by-line absorption coefficient
            computations. If set to zero, each line is computed over the entire
            frequency range.
          self_broadening_tolerance:
            Numerical tolerance for the self-broadening approximation described
            in S2.3.1 of the AM manual.
        """
        assert pressure.shape == temperature.shape
        assert pressure.shape[0] > 0
        self.pressure = pressure
        self.temperature = temperature
        assert 0 * u.deg <= zenith_angle <= 90 * u.deg
        self.zenith_angle = zenith_angle
        # Output frequency range.
        assert freq_max > freq_min
        assert freq_step > 0 * u.Hz
        assert (freq_max - freq_min) > freq_step
        self.freq_min  = freq_min
        self.freq_max  = freq_max
        self.freq_step = freq_step
        # Tropospheric water vapor scaling
        assert troposphere_h2o_scaling > 0
        self.troposphere_h2o_scaling = troposphere_h2o_scaling
        # Validate that the provided mixing ratios are available in AM. If `None`
        # is provided, then check that the specie is avialable in the standard
        # climatology and interpolate onto the pressure grid. Otherwise simply
        # check that the shapes match.
        if mixing_ratio is None:
            mixing_ratio = {}
        for specie, mr in mixing_ratio.items():
            if specie not in self.valid_species:
                raise ValueError(f"Column type unavailable in AM: {specie}")
            if mr is None:
                cl = CLIMATOLOGIES["us_standard"]
                if specie not in CLIMATOLOGIES["us_standard"].mixing_ratio:
                    raise ValueError(f"Default column type unavailable in climatology data: {specie}")
                # Not okay to mutate a dictionary while looping?
                mixing_ratio[specie] = np.interp(pressure, cl.pressure, cl.mixing_ratio[specie])
            else:
                assert np.all(0 <= mr) and np.all(mr <= 1)
                assert mr.shape == pressure.shape
        self.mixing_ratio = mixing_ratio
        # AM output column validation.
        for out in output_columns:
            if out not in self.valid_output_descriptors:
                raise ValueError(f"Invalid output column: {out}")
        self.output_columns = output_columns
        # Execution acceleration parameters.
        assert tolerance > 0
        self.tolerance = tolerance
        assert self_broadening_tolerance >= 0
        self.self_broadening_tolerance = self_broadening_tolerance

    @classmethod
    def from_climatology(cls, name, **kwargs):
        climatology = CLIMATOLOGIES[name]
        mixing_ratio = {
                k: v
                for k, v in climatology.mixing_ratio.items()
                if k in ("h2o", "o3")
        }
        return cls(
                climatology.pressure,
                climatology.temperature,
                mixing_ratio,
                **kwargs,
        )

    @property
    def output_units(self):
        return {
                n: self.valid_output_descriptors[n].split()[-1]
                for n in self.output_columns
        }

    @property
    def output_descriptor(self):
        # Use the NumPy binary output format "npy" as opposed to "text".
        return "output npy  " + "  ".join([
                self.valid_output_descriptors[n]
                for n in self.output_columns
        ])

    @property
    def increasing_pressure_order(self):
        return self.pressure[1] > self.pressure[0]

    @property
    def config_text(self):
        """
        AM configuration text. The Voigt-Kielkopf lineshape is used for
        pressures less than 1 mbar.

        Layers are labeled by pressure for:
        - thermosphere less than 0.01 mbar
        - mesosphere   from 0.01 to 0.1 mbar
        - stratosphere from 1 to 100 mbar
        - troposphere  greater than 100 mbar
        """
        # NOTE The astropy units have string representation with the units
        # appended, so will print as, e.g., "5.0 K". Even though explicit
        # conversions are not necessary, they are done here to make unit errors
        # raise an exception at this stage rather than choke in AM itself.

        ### Example header format
        # output f GHz  tau  Tb K
        # za <zenith_angle> deg
        # tol 1e-4
        # Nscale troposphere h2o <scaling_factor>
        # T0 2.7 K
        ###
        fmin  = self.freq_min.to("GHz").value
        fmax  = self.freq_max.to("GHz").value
        fstep = self.freq_step.to("MHz").value
        header = [
                "? File automatically generated by 'amwrap'.",
                f"f {fmin} GHz  {fmax} GHz  {fstep} MHz",
                self.output_descriptor,
                f"za {self.zenith_angle.to('deg').value} deg",
                f"tol {self.tolerance:1.4e}",
                f"selfbroad_vmr_tol {self.self_broadening_tolerance}",
                f"T0 {self.background_temperature.to('K').value} K",
                f"Nscale troposphere h2o {self.troposphere_h2o_scaling}",
        ]

        ### Example layer format:
        # layer mesosphere
        # Pbase 0.1 mbar
        # Tbase 222.3 K
        # lineshape Voigt-Kielkopf
        # column dry_air vmr
        # column h2o vmr 6.46e-06
        # column o3 vmr 1.77e-06
        ###
        layers = []
        for p, t in zip(
                self.pressure.to("mbar").value,
                self.temperature.to("K").value,
            ):
            if p < 0.01:  # mbar
                layer_type = "thermosphere"
            elif p < 1:
                layer_type = "mesosphere"
            elif p < 100:
                layer_type = "stratosphere"
            else:
                layer_type = "troposphere"
            layer = [
                    f"layer {layer_type}",
                    f"Pbase {p} mbar",
                    f"Tbase {t} K",
                    "column dry_air vmr",
            ]
            if layer_type == "mesosphere":
                layer.append("lineshape Voigt-Kielkopf")
            layers.append(layer)
        for specie, mr in self.mixing_ratio.items():
            for v, layer in zip(mr.to("").value, layers):
                layer.append(f"column {specie} vmr {v:1.4e}")
        # AM requires that pressure levels be specified from low- to
        # high-pressure. The inputs are typically ordered by increasing
        # altitude, so need to be reversed.
        if not self.increasing_pressure_order:
            layers.reverse()
        layers.insert(0, header)
        return "\n\n".join("\n".join(l) for l in layers)

    def write_config(self, outname="output", ext="amc"):
        if ext in (None, ""):
            filen = outname
        else:
            filen = f"{outname}.{ext}"
        with open(filen, "w") as f:
            f.write(self.config_text)

    def _parse_output(self, result):
        # The return code from AM will be 1 if any warnings were raised in
        # the calculations. A common warning is that lines are unresolved
        # at low pressures/high altitudes.
        warnings_returned = bool(result.returncode)
        # The array output is stored in the NumPy binary save-file ".npy"
        # format equivalent to what is generated by the `np.save` function.
        output_data = np.load(BytesIO(result.stdout))
        columns = [s.replace(" ", "_") for s in self.output_columns]
        df = pd.DataFrame(output_data, columns=columns)
        df.attrs["units"] = self.output_units
        df.attrs["stderr"] = result.stderr.decode()
        df.attrs["warnings?"] = warnings_returned
        df.attrs["species"] = list(self.mixing_ratio.keys())
        return df

    def run(self, parallel=False):
        """
        Execute AM for the given configuration.

        Args:
          parallel:
            Whether to call an AM executable that has been compiled for serial
            computations or multi-threaded parallel computation using OpenMP. 

        Returns:
          `pandas.DataFrame` of the output with other metadata set in the
          `attrs` attribute. Metadata includes output from `STDERR`, whether
          warnings were raised, and the units for each output column.
        """
        am = AM_PARALLEL if parallel else AM_SERIAL
        if not am.is_callable:
            raise RuntimeError(f"{am.name} is not callable: {am.exec_name}")
        # Call with subprocess and capture outputs to 'stdout' and 'stderr'.
        try:
            result = subprocess.run(
                    [am.exec_name, "-"],
                    env=ENV,
                    input=self.config_text.encode(),
                    capture_output=True,
            )
        except FileNotFoundError as e:
            raise RuntimeError(f"Could not call `{am.exec_name}`: {e}")
        return self._parse_output(result)


