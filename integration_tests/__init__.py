import os
from pathlib import Path

import pytest


if "CNS_EXEC" in os.environ:
    CNS_EXEC = os.environ["CNS_EXEC"]
else:
    CNS_EXEC = None


DATA_DIR = Path(Path(__file__).parent.parent / "examples")

has_cns = pytest.mark.skipif(
    CNS_EXEC is None, reason="CNS_EXEC environment variable not set"
)

tests_path = Path(__file__).resolve().parents[0]
GOLDEN_DATA = Path(tests_path, "golden_data")

# Placeholder decorator to check if ATTRACT is installed
# TODO: Implement
ATTRACT_IS_INSTALLED = False
has_attract = pytest.mark.skipif(
    not ATTRACT_IS_INSTALLED, reason="ATTRACT not installed"
)
