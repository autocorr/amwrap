import subprocess
from pathlib import Path

from hatchling.builders.hooks.plugin.interface import BuildHookInterface


class CustomBuildHook(BuildHookInterface):
    PLUGIN_NAME = "custom"

    def initialize(self, version, build_data):
        root    = Path(self.root)
        src_dir = root / "src" / "am-14.0"
        bin_dir = root / "amwrap" / "bin"
        bin_dir.mkdir(exist_ok=True)

        # The wheel bundles compiled binaries, so it is platform-specific.
        build_data["pure_python"] = False
        build_data["infer_tag"] = True

        # Clean any stale objects (e.g. OpenMP-compiled .o from a prior parallel
        # build) so the serial link doesn't pull in unresolved GOMP_* symbols.
        subprocess.run(["make", "clean"], cwd=src_dir, check=False)

        # Serial build — required
        subprocess.check_call(["make", "-j", "serial"], cwd=src_dir)
        (src_dir / "am").rename(bin_dir / "am-serial")
        subprocess.check_call(["make", "clean"], cwd=src_dir)

        # Parallel (OpenMP) build — optional
        try:
            subprocess.check_call(["make", "-j", "am"], cwd=src_dir)
            (src_dir / "am").rename(bin_dir / "am")
        except subprocess.CalledProcessError:
            pass
        subprocess.run(["make", "clean"], cwd=src_dir, check=False)

        # Include compiled binaries in the build (not VCS-tracked)
        build_data["artifacts"].append("amwrap/bin/am-serial")
        if (bin_dir / "am").exists():
            build_data["artifacts"].append("amwrap/bin/am")
