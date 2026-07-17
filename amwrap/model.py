# -*- coding: utf-8 -*-
"""
AM model configuration, execution, and output parsing.
"""

from io import BytesIO
from typing import Dict

import numpy as np
import pandas as pd

from astropy import units as u

from . import driver
# Access CLIMATOLOGIES as a module attribute (not a from-import) so that the
# lazy loading in `climatology.py` is not triggered at import time. Aliased
# because `from_climatology` uses `climatology` as a local variable name.
from . import climatology as _climatology
from .climatology import Climatology
from .thermo import altitude_from_pressure


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
    water_cloud_type = "lwp_abs_Rayleigh"
    ice_cloud_type = "iwp_abs_Rayleigh"

    @u.quantity_input
    def __init__(
                self,
                pressure: u.Quantity["pressure"],  # noqa: F821
                temperature: u.Quantity[u.deg_C] | u.Quantity["temperature"],  # noqa: F821
                mixing_ratio: Dict[str, u.Quantity["dimensionless"]|None]|None=None,  # noqa: F821
                water_cloud: u.Quantity["surface_mass_density"]|None=None,  # noqa: F821
                ice_cloud: u.Quantity["surface_mass_density"]|None=None,  # noqa: F821
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
          water_cloud:
            Liquid water cloud mass surface density at the base of the layer.
          ice_cloud:
            Ice water cloud mass surface density at the base of the layer.
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
        if pressure.shape != temperature.shape:
            raise ValueError(
                    f"Shape mismatch: {pressure.shape=} != {temperature.shape=}")
        if pressure.shape[0] < 1:
            raise ValueError("At least one layer is required.")
        self.pressure = pressure
        self.temperature = temperature
        if not (0 * u.deg <= zenith_angle <= 90 * u.deg):
            raise ValueError(f"Zenith angle out of [0, 90] deg: {zenith_angle}")
        self.zenith_angle = zenith_angle
        # Output frequency range.
        if not freq_max > freq_min:
            raise ValueError(f"Require freq_max > freq_min: {freq_max=}, {freq_min=}")
        if not freq_step > 0 * u.Hz:
            raise ValueError(f"Require freq_step > 0: {freq_step=}")
        if not (freq_max - freq_min) > freq_step:
            raise ValueError(f"Require freq_max - freq_min > freq_step: {freq_step=}")
        self.freq_min  = freq_min
        self.freq_max  = freq_max
        self.freq_step = freq_step
        # Tropospheric water vapor scaling
        if not troposphere_h2o_scaling > 0:
            raise ValueError(f"Require troposphere_h2o_scaling > 0: {troposphere_h2o_scaling=}")
        self.troposphere_h2o_scaling = troposphere_h2o_scaling
        # Validate that the provided mixing ratios are available in AM. If `None`
        # is provided, then check that the specie is avialable in the standard
        # climatology and interpolate onto the pressure grid. Otherwise simply
        # check that the shapes match.
        if mixing_ratio is None:
            mixing_ratio = {}
        # Build a fresh dict with copied arrays so neither the caller's input nor
        # the shared CLIMATOLOGIES data can be mutated through this Model.
        resolved_mixing_ratio = {}
        for specie, mr in mixing_ratio.items():
            if specie not in self.valid_species:
                raise ValueError(f"Column type unavailable in AM: {specie}")
            if mr is None:
                cl = _climatology.CLIMATOLOGIES["us_standard"]
                if specie not in cl.mixing_ratio:
                    raise ValueError(f"Default column type unavailable in climatology data: {specie}")
                xp = cl.pressure[::-1].to(cl.pressure.unit).value
                fp = cl.mixing_ratio[specie][::-1].to("").value
                p  = pressure.to(cl.pressure.unit).value
                resolved_mixing_ratio[specie] = np.interp(p, xp, fp) * u.dimensionless_unscaled
            else:
                if not (np.all(0 <= mr) and np.all(mr <= 1)):
                    raise ValueError(f"Mixing ratio outside [0, 1] for {specie!r}")
                if mr.shape != pressure.shape:
                    raise ValueError(
                            f"Shape mismatch for {specie!r}: {mr.shape} != {pressure.shape}")
                resolved_mixing_ratio[specie] = mr.copy()
        self.mixing_ratio = resolved_mixing_ratio
        # Liquid and ice water clouds.
        for label, cloud in (("water_cloud", water_cloud), ("ice_cloud", ice_cloud)):
            if cloud is None:
                continue
            if cloud.shape != pressure.shape:
                raise ValueError(
                        f"Shape mismatch for {label}: {cloud.shape} != {pressure.shape}")
            if np.any(cloud < 0 * cloud.unit):
                raise ValueError(f"Negative {label} value.")
        self.water_cloud = water_cloud
        self.ice_cloud = ice_cloud
        # AM output column validation.
        for out in output_columns:
            if out not in self.valid_output_descriptors:
                raise ValueError(f"Invalid output column: {out}")
        self.output_columns = output_columns
        # Execution acceleration parameters.
        if not tolerance > 0:
            raise ValueError(f"Require tolerance > 0: {tolerance=}")
        self.tolerance = tolerance
        if not self_broadening_tolerance >= 0:
            raise ValueError(f"Require self_broadening_tolerance >= 0: {self_broadening_tolerance=}")
        self.self_broadening_tolerance = self_broadening_tolerance

    @classmethod
    def from_climatology(
                cls,
                clim: str | Climatology,
                **kwargs
        ):
        match clim:
            case str(name):
                climatology = _climatology.CLIMATOLOGIES[name]
            case Climatology():
                climatology = clim
            case _:
                raise ValueError(f"Invalid {clim=}")
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
    def altitude(self):
        return altitude_from_pressure(self.pressure)

    @property
    def output_descriptor(self):
        # Use the NumPy binary output format "npy" as opposed to "text".
        return "output npy  " + "  ".join([
                self.valid_output_descriptors[n]
                for n in self.output_columns
        ])

    @property
    def increasing_pressure_order(self):
        return len(self.pressure) > 1 and self.pressure[1] > self.pressure[0]

    @property
    def config_text(self):
        """
        AM configuration text. The Voigt-Kielkopf lineshape is used for
        pressures less than 1 mbar.

        Layers are labeled by pressure for:
        - thermosphere less than 0.01 mbar
        - mesosphere   from 0.01 to 1 mbar
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
                self.temperature.to("K", equivalencies=u.temperature()).value,
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
            # The Voigt-Kielkopf lineshape is used for pressures less than 1
            # mbar, covering both the mesosphere and thermosphere.
            if p < 1:  # mbar
                layer.append("lineshape Voigt-Kielkopf")
            layers.append(layer)
        for specie, mr in self.mixing_ratio.items():
            for v, layer in zip(mr.to("").value, layers):
                layer.append(f"column {specie} vmr {v:1.4e}")
        if (wc := self.water_cloud) is not None:
            for v, layer in zip(wc.to("kg m-2").value, layers):
                if v > 0:
                    layer.append(f"column lwp_abs_Rayleigh {v} kg*m^-2")
        if (ic := self.ice_cloud) is not None:
            for v, layer in zip(ic.to("kg m-2").value, layers):
                if v > 0:
                    layer.append(f"column iwp_abs_Rayleigh {v} kg*m^-2")
        # AM requires that pressure levels be specified from low- to
        # high-pressure. The inputs are typically ordered by increasing
        # altitude, so need to be reversed.
        if not self.increasing_pressure_order:
            layers.reverse()
        layers.insert(0, header)
        return "\n\n".join("\n".join(block) for block in layers)

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
        # at low pressures/high altitudes. Higher codes (or empty stdout)
        # indicate a hard error, with the cause on stderr.
        if not result.stdout or result.returncode > 1:
            raise RuntimeError(f"AM failed (exit {result.returncode}): "
                               f"{result.stderr.decode()}")
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

    def run(self, parallel=False, cache_dir=None):
        """
        Execute AM for the given configuration.

        Args:
          parallel:
            Whether to call an AM executable that has been compiled for serial
            computations or multi-threaded parallel computation using OpenMP.
          cache_dir:
            Override the AM disk cache directory for this run, setting
            ``AM_CACHE_PATH`` to the given path in the subprocess environment.
            Useful for isolating cache directories across parallel worker
            processes to avoid file-locking contention.

        Returns:
          `pandas.DataFrame` of the output with other metadata set in the
          `attrs` attribute. Metadata includes output from `STDERR`, whether
          warnings were raised, and the units for each output column.
        """
        am = driver.get_executable(parallel)
        result = driver.run_am(self.config_text, am, cache_dir=cache_dir)
        return self._parse_output(result)
