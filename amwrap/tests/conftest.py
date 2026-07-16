import sys
from pathlib import Path

# Make the repo root importable so `import amwrap` works when running the tests
# in-place from a source checkout without an editable install.
sys.path.insert(0, str(Path(__file__).absolute().parent.parent.parent))
