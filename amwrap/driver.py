# -*- coding: utf-8 -*-
"""
Execution boundary for the AM binary: executable discovery, subprocess
environment, and invocation.
"""

import os
import warnings
import subprocess
from pathlib import Path
from packaging import version


class AmExecutable:
    bin_dir = (Path(__file__).parent / "bin").absolute()

    def __init__(self, name):
        if name not in ("am", "am-serial"):
            raise ValueError(f"Invalid name: {name}")
        self.name = name
        # Try the PATH executable first, then the bundled binary. If neither is
        # callable, mark not callable rather than raising at import time.
        for exec_name in (name, str(self.bin_dir / name)):
            am_version = self._probe_version(exec_name)
            if am_version is not None:
                self.exec_name = exec_name
                self.is_callable = True
                self.version = am_version
                break
        else:
            self.exec_name = str(self.bin_dir / name)
            self.is_callable = False
            self.version = None

    @staticmethod
    def _probe_version(exec_name):
        try:
            result = subprocess.run([exec_name, "-v"], capture_output=True)
        except (OSError, subprocess.SubprocessError):
            return None
        try:
            return version.parse(result.stdout.decode().split()[2])
        except (IndexError, UnicodeDecodeError, version.InvalidVersion):
            return None

# The executable singletons are probed lazily on first attribute access (PEP
# 562) rather than at import time, so importing the package does not spawn
# subprocesses. Results are cached into module globals, after which
# `__getattr__` is no longer consulted.
_LAZY_NAMES = ("AM_PARALLEL", "AM_SERIAL", "NO_AM_CALLABLE", "BOTH_AM_CALLABLE")


def _probe():
    g = globals()
    if "AM_PARALLEL" in g:
        return
    # FIXME Use configuration system for "am"/"am-serial" executable names.
    g["AM_PARALLEL"] = am_parallel = AmExecutable("am")
    g["AM_SERIAL"]   = am_serial   = AmExecutable("am-serial")
    g["NO_AM_CALLABLE"] = no_am = (
            not am_parallel.is_callable and not am_serial.is_callable)
    g["BOTH_AM_CALLABLE"] = both_am = (
            am_parallel.is_callable and am_serial.is_callable)
    if no_am:
        warnings.warn("No executable callable for AM.", UserWarning)
    if both_am and (am_parallel.version != am_serial.version):
        warnings.warn(
                "`am` and `am-serial` version mismatch: "
                f"{am_parallel.version} & {am_serial.version}",
                UserWarning,
        )


def __getattr__(name):
    if name in _LAZY_NAMES:
        _probe()
        return globals()[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(_LAZY_NAMES))


# Set environment variables for am
ENV = os.environ.copy()
CACHE_DIR = Path("/dev/shm/")
if CACHE_DIR.exists():
    ENV["AM_CACHE_PATH"] = str(CACHE_DIR)


def get_executable(parallel=False):
    """
    Select the AM executable for serial or OpenMP-parallel computation,
    raising `RuntimeError` if it is not callable.
    """
    _probe()
    am = AM_PARALLEL if parallel else AM_SERIAL
    if not am.is_callable:
        raise RuntimeError(f"{am.name} is not callable: {am.exec_name}")
    return am


def run_am(config_text, executable, env=None, cache_dir=None):
    """
    Run an AM executable on the given configuration text.

    Args:
      config_text:
        AM configuration file contents, passed to AM on `STDIN`.
      executable:
        `AmExecutable` instance (or any object with an ``exec_name``
        attribute) to invoke.
      env:
        Environment mapping for the subprocess; defaults to the frozen
        module-level `ENV`.
      cache_dir:
        Override the AM disk cache directory for this run, setting
        ``AM_CACHE_PATH`` to the given path in the subprocess environment.

    Returns:
      `subprocess.CompletedProcess` with captured ``stdout`` and ``stderr``.
    """
    if env is None:
        env = ENV
    if cache_dir is not None:
        env = {**env, "AM_CACHE_PATH": str(cache_dir)}
    # Call with subprocess and capture outputs to 'stdout' and 'stderr'.
    try:
        return subprocess.run(
                [executable.exec_name, "-"],
                env=env,
                input=config_text.encode(),
                capture_output=True,
        )
    except FileNotFoundError as e:
        raise RuntimeError(f"Could not call `{executable.exec_name}`: {e}")
