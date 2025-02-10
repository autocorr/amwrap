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
#   "iwp_abs_Rayleigh" column types.
# - Parse diagnostic data from the STDERR output into the DataFrame.
# - Perform interpolations/extrapolation on vertical profiles
# - Have the PWV be configurable as a both a direct input in the mixing ratio as
#   well an exponential function with a scale height.

import os
import warnings
import subprocess
from io import BytesIO
from pathlib import Path
from datetime import datetime
from packaging import version

import numpy as np
import pandas as pd

from astropy import units as u
from astropy import constants as c


# Test whether the `am` and `am-serial` executables are callable.
# FIXME Use configuration system for "am" executable name.
def _test_callable(name):
    try:
        result = subprocess.run([f"{name}", "-v"], capture_output=True)
        v = version.parse(result.stdout.decode().split()[2])
        return True, v
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        warnings.warn(f"`{name}` executable not callable.", UserWarning)
        return False, None

CALLABLE, AM_VERSION = _test_callable("am")
CALLABLE_SERIAL, AM_SERIAL_VERSION = _test_callable("am-serial")

NO_AM_CALLABLE = not CALLABLE and not CALLABLE_SERIAL
BOTH_AM_CALLABLE = not NO_AM_CALLABLE
if NO_AM_CALLABLE:
    warnings.warn("No executable callable for AM.", UserWarning)
if BOTH_AM_CALLABLE and (AM_VERSION != AM_SERIAL_VERSION):
    raise RuntimeError(f"`am` and `am-serial` version mismatch: {AM_VERSION} & {AM_SERIAL_VERSION}")

# Set environment variables for am
ENV = os.environ.copy()
CACHE_DIR = Path("/dev/shm/")
if CACHE_DIR.exists():
    ENV["AM_CACHE_PATH"] = CACHE_DIR

MOD_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
MASS_DRY_AIR = 28.9644
MASS_WATER =   18.01528


def mixing_ratio_from_relative_humidity(cls, pressure, temperature, relative_humidity):
    from metpy.units.pint import Quantity
    from metpy.calc import mixing_ratio_from_relative_humidity
    p  = pressure.to("hPa").value * Quantity("hPa")
    t  = temperature.to("deg_C", equivalencies=u.temperature()).value * Quantity("degC")
    rh = relative_humidity.value * Quantity("dimensionless")
    mr = relative_humdity_from_mixing_ratio(p, t, rh)
    # Returned value is mass mixing ratio, convert to volumetric mixing ratio using masses
    return mr.m * 1e6 * MASS_AIR / MASS_WATER * u.dimensionless_unscaled


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
        data = np.loadtxt(MOD_DIR / "climatology" / f"{name}.dat")
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
        self.h2o_mixing_ratio = data[:, 4] * from_ppm
        self.co2_mixing_ratio = data[:, 5] * from_ppm
        self.o3_mixing_ratio  = data[:, 6] * from_ppm
        self.n2o_mixing_ratio = data[:, 7] * from_ppm
        self.co_mixing_ratio  = data[:, 8] * from_ppm
        self.ch4_mixing_ratio = data[:, 9] * from_ppm
        self.o2_mixing_ratio  = data[:,10] * from_ppm
        # TODO Integrate water vapor to PWV, useful for the tropospheric
        # relative scaling in AM calculations.

    @classmethod
    def midlat_from_datetime(cls, dt):
        if dt.month in range(4, 10):
            return cls("midlatitude_summer")
        else:
            return cls("midlatitude_winter")

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

    @u.quantity_input
    def __init__(self,
                pressure: u.Quantity["pressure"],
                temperature: u.Quantity["temperature"],
                h2o_mixing_ratio: u.Quantity["dimensionless"],
                zenith_angle: u.Quantity["angle"]=0*u.deg,
                freq_min: u.Quantity["frequency"]=18*u.GHz,
                freq_max: u.Quantity["frequency"]=26.5*u.GHz,
                freq_step: u.Quantity["frequency"]=10*u.MHz,
                troposphere_h2o_scaling=1.0,
                do_ozone=True,
                o3_mixing_ratio: u.Quantity["dimensionless"]|None=None,
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
          h2o_mixing_ratio:
            Volumetric mixing ratio of water.
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
          do_ozone:
            Whether to calculate ozone lines.
          o3_mixing_ratio:
            Volumetric mixing ratio of ozone.
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
        # FIXME AM wants layers listed in low-to-high pressure order, mark for
        # ascending or descending and then feed to iterator in config generation.
        assert pressure.shape == temperature.shape == h2o_mixing_ratio.shape
        assert pressure.shape[0] > 1
        self.pressure = pressure
        self.temperature = temperature
        self.h2o_mixing_ratio = h2o_mixing_ratio
        assert 0 * u.deg <= zenith_angle <= 90 * u.deg
        self.zenith_angle = zenith_angle
        # Output frequency range.
        assert freq_max > freq_min
        assert freq_step > 0 * u.Hz
        self.freq_min  = freq_min
        self.freq_max  = freq_max
        self.freq_step = freq_step
        # Tropospheric water vapor scaling
        assert troposphere_h2o_scaling > 0
        self.troposphere_h2o_scaling = troposphere_h2o_scaling
        # If ozone is requested but not supplied, use the US standard
        # atmosphere, interpolating to the same pressure axis if needed.
        self.do_ozone = do_ozone
        if do_ozone and o3_mixing_ratio is None:
            climatology = CLIMATOLOGIES["us_standard"]
            o3_mixing_ratio = climatology.o3_mixing_ratio
            if not np.allclose(pressure, climatology.pressure):
                # Dissimilar axes, interpolate.
                self.o3_mixing_ratio = np.interp(pressure, climatology.pressure, o3_mixing_ratio)
            else:
                # Same pressure axis, use directly.
                self.o3_mixing_ratio = o3_mixing_ratio
        else:
            assert h2o_mixing_ratio.shape == o3_mixing_ratio.shape
            self.o3_mixing_ratio = o3_mixing_ratio
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
        return cls(
                climatology.pressure,
                climatology.temperature,
                climatology.h2o_mixing_ratio,
                o3_mixing_ratio=climatology.o3_mixing_ratio,
                **kwargs,
        )

    @property
    def output_units(self):
        return [
                self.valid_output_descriptors[n].split()[-1]
                for n in self.output_columns
        ]

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
        for p, t, w_mr in zip(
                self.pressure.to("mbar").value,
                self.temperature.to("K").value,
                self.h2o_mixing_ratio.to("").value,
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
                    f"column h2o vmr {w_mr:1.4e}",
            ]
            if layer_type == "mesosphere":
                layer.append("lineshape Voigt-Kielkopf")
            layers.append(layer)
        if self.do_ozone:
            for item, o_mr in zip(layers, self.o3_mixing_ratio.to("").value):
                item.append(f"column o3 vmr {o_mr:1.4e}")
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
        df.attrs["ozone?"] = self.do_ozone
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
        # NOTE Sean found that using `subprocess` was slower than writing and
        # reading files. Maybe an encoding thing? Try with `subprocess` first,
        # and if slow then move to file IO on /dev/shm.
        executable_name = "am" if parallel else "am-serial"
        if parallel and not CALLABLE:
            raise RuntimeError("`am` is not callable.")
        if not parallel and not CALLABLE_SERIAL:
            raise RuntimeError("`am-serial` is not callable")
        # Call with subprocess and capture outputs to 'stdout' and 'stderr'.
        try:
            result = subprocess.run(
                    [executable_name, "-"],
                    env=ENV,
                    input=self.config_text.encode(),
                    capture_output=True,
            )
        except FileNotFoundError as e:
            raise RuntimeError(f"Could not call `{executable_name}`: {e}")
        return self._parse_output(result)


