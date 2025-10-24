import os
import platform
from pathlib import Path

import pytest

from haddock.libs.libgrid import ping_dirac

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

try:
    import mpi4py

    MPI_ENABLED = True
except ImportError:
    MPI_ENABLED = False

has_mpi = pytest.mark.skipif(not MPI_ENABLED, reason="MPI is not enabled")
has_grid = pytest.mark.skipif(not ping_dirac(), reason="Dirac not reachable")
is_linux_x86_64 = pytest.mark.skipif(
    platform.system().lower() != "linux" or platform.machine().lower() != "x86_64",
    reason="Only runs on x86_64 Linux systems",
)
is_not_linux_x86_64 = pytest.mark.skipif(
    platform.system().lower() == "linux" and platform.machine().lower() == "x86_64",
    reason="Only runs on non-x86_64-linux systems",
)
