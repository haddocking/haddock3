import os
from pathlib import Path

from haddock.libs.libgrid import ping_dirac
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

try:
    import mpi4py

    MPI_ENABLED = True
except ImportError:
    MPI_ENABLED = False

has_mpi = pytest.mark.skipif(not MPI_ENABLED, reason="MPI is not enabled")
has_grid = pytest.mark.skipif(not ping_dirac(), reason="Dirac not reachable")
