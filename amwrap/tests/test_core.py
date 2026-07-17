
import types
import pytest

from datetime import datetime

import numpy as np
from astropy import units as u

import amwrap
from amwrap import (
        Climatology,
        Model,
        AmExecutable,
        interp_by_pressure,
        altitude_from_pressure,
        mixing_ratio_from_relative_humidity,
        precipitable_water,
)


# Tests that shell out to the compiled AM binary. Everything else in this file
# is pure Python and runs without a compiled AM.
needs_am = pytest.mark.skipif(
        not amwrap.BOTH_AM_CALLABLE,
        reason="AM binary not callable",
)


@pytest.fixture
def f_cl():
    return Climatology("us_standard")


@pytest.fixture
def f_m():
    return Model.from_climatology("us_standard")


@pytest.fixture
def f_layers():
    """A four-layer model spanning all four pressure regimes."""
    return Model(
            pressure=[0.005, 0.5, 50, 500] * u.mbar,
            temperature=[250, 260, 270, 290] * u.K,
            mixing_ratio={"h2o": [1e-6, 2e-6, 3e-6, 4e-6] * u.dimensionless_unscaled},
    )


class TestClimatology:
    def test_read(self):
        assert Climatology("us_standard")
        assert Climatology("midlatitude_winter")
        assert Climatology("midlatitude_summer")
        assert Climatology.midlat_from_datetime(datetime(2025, 1, 1))

    def test_mixing_ratios(self, f_cl):
        assert len(f_cl.mixing_ratio) == 23
        assert np.isclose(f_cl.column_density("h2o").value, 4.80543e22)
        assert np.isclose(f_cl.column_density("o3").value,  9.28091e18)

    def test_pwv(self, f_cl):
        assert np.isclose(f_cl.pwv.value, 14.377778)

    def test_dobson_unit(self, f_cl):
        assert np.isclose(f_cl.dobson_unit("o3"), 345.01522)


class TestModel:
    def test_from_climatology(self):
        assert Model.from_climatology("us_standard")

    def test_init(self, f_cl):
        assert Model(
                f_cl.pressure,
                f_cl.temperature,
                f_cl.mixing_ratio,
        )

    def test_freq(self):
        with pytest.raises(ValueError):
            Model.from_climatology("us_standard", freq_min=10*u.GHz, freq_max=8*u.GHz)
        with pytest.raises(ValueError):
            Model.from_climatology("us_standard", freq_min=1*u.GHz, freq_max=2*u.GHz, freq_step=10*u.GHz)
        with pytest.raises(u.UnitsError):
            Model.from_climatology("us_standard", freq_min=5*u.m)

    @needs_am
    def test_run(self, f_m):
        assert f_m.run().shape
        assert f_m.run(parallel=True).shape
        df = f_m.run()
        assert df.shape == (851, 4)
        assert np.isclose(df.loc[0, "frequency"], 18.0)
        assert np.isclose(df.loc[0, "brightness_temperature"], 11.579323)
        assert np.isclose(df.loc[0, "opacity"], 0.0335504)
        assert np.isclose(df.loc[0, "delay"], 2301.044259)
        assert len(df.attrs["units"]) == 4
        assert df.attrs["stderr"]
        assert df.attrs["warnings?"]
        assert len(df.attrs["species"]) == 2

    @needs_am
    def test_run_cache_dir(self, f_m, tmp_path):
        import amwrap
        original_cache_path = amwrap.ENV.get("AM_CACHE_PATH")
        df = f_m.run(cache_dir=tmp_path)
        assert df.shape == (851, 4)
        assert any(tmp_path.glob("am_*")), "AM wrote no cache files to cache_dir"
        assert amwrap.ENV.get("AM_CACHE_PATH") == original_cache_path

    @needs_am
    def test_run_single_layer(self):
        df = Model([500] * u.mbar, [290] * u.K).run()
        assert df.shape[0] > 0

    @needs_am
    def test_run_cloud(self, f_cl):
        unit = u.kg / u.m**2
        water_cloud = np.zeros_like(f_cl.pressure).value * unit
        ice_cloud = water_cloud.copy()
        water_cloud[2] += 0.1 * unit  # ~2 km
        ice_cloud[10]  += 0.1 * unit  # ~10 km
        m = Model(f_cl.pressure, f_cl.temperature, water_cloud=water_cloud, ice_cloud=ice_cloud)
        df = m.run()
        assert np.isclose(df.loc[0, "brightness_temperature"], 7.985342)


class TestConfigText:
    """Direct unit tests for the AM configuration rendering (no binary)."""

    def test_layer_classification(self, f_layers):
        text = f_layers.config_text
        assert "layer thermosphere" in text
        assert "layer mesosphere" in text
        assert "layer stratosphere" in text
        assert "layer troposphere" in text

    def test_layer_boundaries(self):
        # Pressures placed just inside each regime boundary, as a single model
        # so each layer's own block can be checked. (config_text needs >1 layer.)
        cases = {
                0.001: "thermosphere",  # < 0.01
                0.009: "thermosphere",
                0.01:  "mesosphere",    # [0.01, 1)
                0.5:   "mesosphere",
                1.0:   "stratosphere",  # [1, 100)
                50.0:  "stratosphere",
                100.0: "troposphere",   # >= 100
                900.0: "troposphere",
        }
        pressures = list(cases.keys())
        m = Model(pressures * u.mbar, [270] * len(pressures) * u.K)
        blocks = m.config_text.split("\n\n")
        for p, expected in cases.items():
            block = next(b for b in blocks if f"Pbase {p} mbar" in b)
            assert f"layer {expected}" in block, f"{p} mbar -> {expected}"

    def test_voigt_kielkopf_below_one_mbar(self, f_layers):
        # The lineshape line should follow both the thermosphere and mesosphere
        # layers (p < 1 mbar) and no others.
        blocks = f_layers.config_text.split("\n\n")
        for block in blocks:
            has_ls = "lineshape Voigt-Kielkopf" in block
            if "layer thermosphere" in block or "layer mesosphere" in block:
                assert has_ls, f"missing lineshape:\n{block}"
            elif "layer stratosphere" in block or "layer troposphere" in block:
                assert not has_ls, f"unexpected lineshape:\n{block}"

    def test_pressure_order_reversal(self):
        # Input ordered by increasing altitude (decreasing pressure) must be
        # emitted low- to high-pressure for AM.
        m = Model([500, 100, 10, 1] * u.mbar, [290, 270, 260, 250] * u.K)
        assert not m.increasing_pressure_order
        text = m.config_text
        positions = [text.index(f"Pbase {p} mbar") for p in (1.0, 10.0, 100.0, 500.0)]
        assert positions == sorted(positions), "layers not ordered low->high pressure"

    def test_pressure_order_already_increasing(self):
        m = Model([1, 10, 100, 500] * u.mbar, [250, 260, 270, 290] * u.K)
        assert m.increasing_pressure_order
        text = m.config_text
        positions = [text.index(f"Pbase {p} mbar") for p in (1.0, 10.0, 100.0, 500.0)]
        assert positions == sorted(positions)

    def test_cloud_column_gating(self):
        unit = u.kg / u.m**2
        p = [10, 100, 500] * u.mbar
        t = [250, 270, 290] * u.K
        water = [0, 0.1, 0] * unit
        ice = [0.2, 0, 0] * unit
        text = Model(p, t, water_cloud=water, ice_cloud=ice).config_text
        # Only the seeded layers carry a cloud column.
        assert text.count("lwp_abs_Rayleigh") == 1
        assert text.count("iwp_abs_Rayleigh") == 1
        blocks = text.split("\n\n")
        water_block = next(b for b in blocks if "Pbase 100.0" in b)
        ice_block = next(b for b in blocks if "Pbase 10.0" in b)
        assert "lwp_abs_Rayleigh" in water_block
        assert "iwp_abs_Rayleigh" in ice_block

    def test_header_fields(self):
        m = Model(
                [400, 500] * u.mbar,
                [280, 290] * u.K,
                zenith_angle=30 * u.deg,
                freq_min=20 * u.GHz,
                freq_max=30 * u.GHz,
                freq_step=5 * u.MHz,
                troposphere_h2o_scaling=0.8,
                tolerance=1e-5,
                self_broadening_tolerance=0.01,
        )
        text = m.config_text
        assert "f 20.0 GHz  30.0 GHz  5.0 MHz" in text
        assert "za 30.0 deg" in text
        assert "tol 1.0000e-05" in text
        assert "selfbroad_vmr_tol 0.01" in text
        assert "T0 2.725 K" in text
        assert "Nscale troposphere h2o 0.8" in text

    def test_mixing_ratio_format(self):
        m = Model(
                [400, 500] * u.mbar,
                [280, 290] * u.K,
                mixing_ratio={"h2o": [1.23456e-6, 2e-6] * u.dimensionless_unscaled},
        )
        assert "column h2o vmr 1.2346e-06" in m.config_text

    def test_output_descriptor(self, f_layers):
        assert f_layers.output_descriptor.startswith("output npy")
        assert "Tb K" in f_layers.output_descriptor

    def test_single_layer_config(self):
        m = Model([500] * u.mbar, [290] * u.K)
        assert m.increasing_pressure_order is False
        text = m.config_text
        assert "Pbase 500.0 mbar" in text
        assert "layer troposphere" in text


class TestUtilities:
    """Module-level pure functions."""

    def test_altitude_from_pressure_reference(self):
        import amwrap as A
        assert np.isclose(altitude_from_pressure(A.STD_PRESSURE).to("km").value, 0.0)

    def test_altitude_from_pressure_monotonic(self):
        alt = altitude_from_pressure([1000, 500, 100] * u.hPa)
        assert alt.unit.is_equivalent(u.km)
        assert alt[0] < alt[1] < alt[2]

    def test_interp_identity_without_base(self):
        vals = np.array([1.0, 2.0, 3.0])
        out = interp_by_pressure(vals)
        assert np.array_equal(out, vals)

    def test_interp_shape_mismatch(self):
        with pytest.raises(ValueError):
            interp_by_pressure(np.array([1.0, 2.0]), pressure=[1, 2, 3] * u.mbar)

    def test_interp_base_without_pressure(self):
        with pytest.raises(ValueError):
            interp_by_pressure(np.array([1.0, 2.0]), pressure_base=5 * u.mbar)

    def test_interp_base_out_of_bounds(self):
        vals = np.array([1.0, 2.0, 3.0])
        p = [100, 200, 300] * u.mbar
        with pytest.raises(ValueError):
            interp_by_pressure(vals, pressure=p, pressure_base=500 * u.mbar)
        with pytest.raises(ValueError):
            interp_by_pressure(vals, pressure=p, pressure_base=50 * u.mbar)

    def test_interp_clips_to_base(self):
        # Profile ordered by decreasing pressure (increasing altitude).
        vals = np.array([10.0, 20.0, 30.0, 40.0])
        p = [1000, 750, 500, 250] * u.mbar
        out = interp_by_pressure(vals, pressure=p, pressure_base=600 * u.mbar)
        # Base falls between 750 and 500 mbar. The retained slice starts at the
        # base value (interpolated at 600 mbar) and keeps all lower-pressure
        # (higher-altitude) layers: [interp@600, 30, 40].
        assert len(out) == 3
        expected_base = np.interp(600, [250, 500, 750, 1000], [40, 30, 20, 10])
        assert np.isclose(out[0], expected_base)
        assert np.isclose(out[-1], 40.0)

    def test_mixing_ratio_from_relative_humidity(self):
        mr = mixing_ratio_from_relative_humidity(
                1000 * u.hPa, 300 * u.K, 0.5 * u.dimensionless_unscaled)
        assert mr.unit.is_equivalent(u.dimensionless_unscaled)
        assert 0 < mr.value < 0.1
        # Higher RH yields a larger mixing ratio.
        mr_hi = mixing_ratio_from_relative_humidity(
                1000 * u.hPa, 300 * u.K, 0.9 * u.dimensionless_unscaled)
        assert mr_hi > mr

    def test_precipitable_water(self):
        p = [1000, 900, 800] * u.hPa
        t = [300, 295, 290] * u.K
        rh = [0.5, 0.5, 0.5] * u.dimensionless_unscaled
        pwv = precipitable_water(p, t, rh)
        assert pwv.unit.is_equivalent(u.mm)
        assert pwv.to("mm").value > 0


class TestModelValidation:
    """Construction-time validation branches."""

    def test_invalid_species(self):
        with pytest.raises(ValueError):
            Model([500] * u.mbar, [290] * u.K,
                  mixing_ratio={"bogus": [1e-6] * u.dimensionless_unscaled})

    def test_none_species_interpolates(self, f_cl):
        m = Model(f_cl.pressure, f_cl.temperature, mixing_ratio={"h2o": None})
        assert m.mixing_ratio["h2o"].shape == f_cl.pressure.shape
        # On the native climatology grid the interpolation is the identity.
        assert np.allclose(
                m.mixing_ratio["h2o"].to("").value,
                f_cl.mixing_ratio["h2o"].to("").value,
        )

    def test_none_species_interpolation_values(self):
        # Interpolate onto a shifted pressure grid and compare to np.interp on
        # the reversed (increasing) climatology grid.
        cl = amwrap.CLIMATOLOGIES["us_standard"]
        target = [900, 500, 100] * u.mbar
        m = Model(target, [270, 260, 250] * u.K, mixing_ratio={"h2o": None})
        xp = cl.pressure[::-1].to("mbar").value
        fp = cl.mixing_ratio["h2o"][::-1].to("").value
        expected = np.interp(target.to("mbar").value, xp, fp)
        assert np.allclose(m.mixing_ratio["h2o"].to("").value, expected)
        # A decreasing-xp (buggy) interpolation would collapse to an endpoint.
        assert not np.allclose(expected, expected[0])

    def test_none_species_not_in_climatology(self):
        with pytest.raises(ValueError):
            # "co" is a valid AM species and in the climatology, but "ho2" is
            # a valid AM species absent from the climatology defaults.
            Model([500] * u.mbar, [290] * u.K, mixing_ratio={"ho2": None})

    def test_mixing_ratio_out_of_range(self):
        with pytest.raises(ValueError):
            Model([500] * u.mbar, [290] * u.K,
                  mixing_ratio={"h2o": [2.0] * u.dimensionless_unscaled})

    def test_zenith_angle_bounds(self):
        with pytest.raises(ValueError):
            Model([500] * u.mbar, [290] * u.K, zenith_angle=95 * u.deg)
        with pytest.raises(ValueError):
            Model([500] * u.mbar, [290] * u.K, zenith_angle=-5 * u.deg)

    def test_invalid_output_column(self):
        with pytest.raises(ValueError):
            Model([500] * u.mbar, [290] * u.K, output_columns=("frequency", "bogus"))

    def test_tolerance_bounds(self):
        with pytest.raises(ValueError):
            Model([500] * u.mbar, [290] * u.K, tolerance=0)
        with pytest.raises(ValueError):
            Model([500] * u.mbar, [290] * u.K, self_broadening_tolerance=-0.1)

    def test_from_climatology_instance(self, f_cl):
        assert Model.from_climatology(f_cl)

    def test_from_climatology_invalid_type(self):
        with pytest.raises(ValueError):
            Model.from_climatology(42)

    def test_output_units_and_altitude(self, f_layers):
        assert f_layers.output_units["frequency"] == "GHz"
        assert f_layers.altitude.unit.is_equivalent(u.km)

    def test_mixing_ratio_not_aliased(self):
        # Mutating a Model's mixing ratio must not corrupt the shared
        # climatology it was resolved from.
        cl = amwrap.CLIMATOLOGIES["us_standard"]
        before = cl.mixing_ratio["h2o"].copy()
        m = Model(cl.pressure, cl.temperature, mixing_ratio={"h2o": None})
        m.mixing_ratio["h2o"] *= 0.5
        assert np.array_equal(cl.mixing_ratio["h2o"].to("").value,
                              before.to("").value)

    def test_input_dict_not_mutated(self):
        provided = {"h2o": [1e-6, 2e-6] * u.dimensionless_unscaled}
        m = Model([400, 500] * u.mbar, [280, 290] * u.K, mixing_ratio=provided)
        m.mixing_ratio["h2o"] *= 0.5
        assert np.allclose(provided["h2o"].to("").value, [1e-6, 2e-6])

    def test_cloud_shape_mismatch(self):
        with pytest.raises(ValueError):
            Model([400, 500] * u.mbar, [280, 290] * u.K,
                  water_cloud=[0.1] * (u.kg / u.m**2))

    def test_cloud_negative(self):
        with pytest.raises(ValueError):
            Model([400, 500] * u.mbar, [280, 290] * u.K,
                  ice_cloud=[0.1, -0.1] * (u.kg / u.m**2))


class TestClimatologyEdge:
    def test_invalid_name(self):
        with pytest.raises(ValueError):
            Climatology("not_a_climatology")

    def test_pressure_base_clip(self):
        full = Climatology("us_standard")
        clipped = Climatology("us_standard", pressure_base=500 * u.mbar)
        assert len(clipped.pressure) < len(full.pressure)
        assert np.isclose(clipped.pressure[0].to("mbar").value, 500.0)
        assert clipped.mixing_ratio["h2o"].shape == clipped.pressure.shape

    def test_midlat_from_datetime_summer(self):
        assert Climatology.midlat_from_datetime(datetime(2025, 7, 1)).name \
                == "midlatitude_summer"

    def test_midlat_from_datetime_winter(self):
        assert Climatology.midlat_from_datetime(datetime(2025, 1, 1)).name \
                == "midlatitude_winter"

    def test_am_executable_invalid_name(self):
        with pytest.raises(ValueError):
            AmExecutable("bogus")


class TestParseOutput:
    """`_parse_output` in isolation, feeding a fake subprocess result."""

    def _fake_result(self, returncode=0):
        from io import BytesIO
        data = np.arange(12, dtype=float).reshape(3, 4)
        buf = BytesIO()
        np.save(buf, data)
        return types.SimpleNamespace(
                returncode=returncode,
                stdout=buf.getvalue(),
                stderr=b"some stderr text",
        )

    def test_parse_output_dataframe(self):
        m = Model([500] * u.mbar, [290] * u.K,
                  mixing_ratio={"h2o": [1e-6] * u.dimensionless_unscaled})
        df = m._parse_output(self._fake_result(returncode=0))
        assert df.shape == (3, 4)
        assert list(df.columns) == ["frequency", "brightness_temperature", "opacity", "delay"]
        assert df.attrs["units"]["frequency"] == "GHz"
        assert df.attrs["stderr"] == "some stderr text"
        assert df.attrs["warnings?"] is False
        assert df.attrs["species"] == ["h2o"]

    def test_parse_output_warnings_flag(self):
        m = Model([500] * u.mbar, [290] * u.K)
        df = m._parse_output(self._fake_result(returncode=1))
        assert df.attrs["warnings?"] is True

    def test_parse_output_hard_error(self):
        m = Model([500] * u.mbar, [290] * u.K)
        result = types.SimpleNamespace(
                returncode=2, stdout=b"", stderr=b"am: config error on line 3")
        with pytest.raises(RuntimeError, match="config error on line 3"):
            m._parse_output(result)

    def test_parse_output_empty_stdout(self):
        m = Model([500] * u.mbar, [290] * u.K)
        result = types.SimpleNamespace(
                returncode=0, stdout=b"", stderr=b"unexpected failure")
        with pytest.raises(RuntimeError, match="unexpected failure"):
            m._parse_output(result)
