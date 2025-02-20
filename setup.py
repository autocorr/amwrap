import os
import subprocess
from pathlib import Path

import setuptools


def build_am():
    """
    This is a major hack to build AM by running commands at the module-scope
    but I wasn't able to get the `setuptools.Command` scheme to work nor
    successfully adapt the solution from this forum post:
      https://discuss.python.org/t/need-to-run-make-when-building-wheel-from-sdist-or-installing-the-latter-how/64417/11
    """
    root_dir = Path.cwd()
    src_dir  = root_dir / "src" / "am-14.0"
    bin_dir  = root_dir / "amwrap" / "bin"
    os.chdir(src_dir)
    # compile serial version
    subprocess.check_call(["make", "-d", "-j", "serial"])
    (src_dir / "am").rename(bin_dir / "am-serial")
    subprocess.check_call(["make", "clean"])
    print("-- Compiled serial version of AM")
    # compile standard/multithreaded version
    subprocess.check_call(["make", "-d", "-j", "am"])
    (src_dir / "am").rename(bin_dir / "am")
    subprocess.check_call(["make", "clean"])
    print("-- Compiled parallel version of AM")
    os.chdir(root_dir)


build_am()
setuptools.setup()

