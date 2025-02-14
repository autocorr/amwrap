
import pytest

from datetime import datetime

import numpy as np
from astropy import units as u

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).absolute().parent.parent.parent))
from amwrap import (Climatology, Model)


@pytest.fixture
def f_cl():
    return Climatology("us_standard")


@pytest.fixture
def f_m():
    return Model.from_climatology("us_standard")


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
        with pytest.raises(AssertionError):
            Model.from_climatology("us_standard", freq_min=10*u.GHz, freq_max=8*u.GHz)
        with pytest.raises(AssertionError):
            Model.from_climatology("us_standard", freq_min=1*u.GHz, freq_max=2*u.GHz, freq_step=10*u.GHz)
        with pytest.raises(u.UnitsError):
            Model.from_climatology("us_standard", freq_min=5*u.m)

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

