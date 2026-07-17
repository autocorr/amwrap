
import io
import sys
import stat
import types
import subprocess
import textwrap

import pytest

import numpy as np
from astropy import units as u

import amwrap
from amwrap import Model, driver


class TestPublicApi:
    def test_all_names_resolve(self):
        for name in amwrap.__all__:
            assert getattr(amwrap, name) is not None

    def test_dir_includes_lazy_names(self):
        for name in ("AM_PARALLEL", "AM_SERIAL", "BOTH_AM_CALLABLE",
                     "NO_AM_CALLABLE", "CLIMATOLOGIES"):
            assert name in dir(amwrap)

    def test_unknown_attribute_raises(self):
        with pytest.raises(AttributeError):
            amwrap.does_not_exist


class TestLazyImport:
    def test_import_is_lazy(self):
        # Run in a fresh interpreter so attribute access from other tests
        # cannot have triggered the probe or the climatology load already.
        script = textwrap.dedent("""
            import amwrap, amwrap.driver, amwrap.climatology
            assert "AM_PARALLEL" not in vars(amwrap.driver)
            assert "CLIMATOLOGIES" not in vars(amwrap.climatology)
            amwrap.BOTH_AM_CALLABLE
            assert "AM_PARALLEL" in vars(amwrap.driver)
            amwrap.CLIMATOLOGIES
            assert "CLIMATOLOGIES" in vars(amwrap.climatology)
        """)
        result = subprocess.run(
                [sys.executable, "-c", script],
                capture_output=True,
        )
        assert result.returncode == 0, result.stderr.decode()


class TestRunAm:
    def make_stub(self, tmp_path, payload, returncode=0, stderr=""):
        """Create a fake AM executable that emits canned bytes on stdout."""
        data_path = tmp_path / "canned.npy"
        data_path.write_bytes(payload)
        stub_path = tmp_path / "am-stub"
        stub_path.write_text(
                "#!/bin/sh\n"
                f"cat {data_path}\n"
                f"echo -n '{stderr}' >&2\n"
                f"exit {returncode}\n"
        )
        stub_path.chmod(stub_path.stat().st_mode | stat.S_IXUSR)
        return types.SimpleNamespace(
                name="am-stub", exec_name=str(stub_path), is_callable=True)

    def test_run_am_injected_executable(self, tmp_path):
        m = Model(
                pressure=[50, 500] * u.mbar,
                temperature=[260, 290] * u.K,
        )
        values = np.arange(20, dtype=float).reshape(5, 4)
        buf = io.BytesIO()
        np.save(buf, values)
        stub = self.make_stub(tmp_path, buf.getvalue())
        result = driver.run_am(m.config_text, stub)
        df = m._parse_output(result)
        assert df.shape == (5, 4)
        assert list(df.columns) == [
                "frequency", "brightness_temperature", "opacity", "delay"]
        assert not df.attrs["warnings?"]

    def test_run_am_missing_executable(self, tmp_path):
        stub = types.SimpleNamespace(
                name="am-stub",
                exec_name=str(tmp_path / "no-such-binary"),
                is_callable=True,
        )
        with pytest.raises(RuntimeError, match="Could not call"):
            driver.run_am("", stub)

    def test_run_am_cache_dir_does_not_mutate_env(self, tmp_path):
        buf = io.BytesIO()
        np.save(buf, np.zeros((1, 4)))
        stub = self.make_stub(tmp_path, buf.getvalue())
        env_before = dict(amwrap.ENV)
        driver.run_am("", stub, cache_dir=tmp_path)
        assert dict(amwrap.ENV) == env_before
