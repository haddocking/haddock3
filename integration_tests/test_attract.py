import os
import shutil
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.sampling.attract import \
    DEFAULT_CONFIG as DEFAULT_ATTRACT_CONFIG
from haddock.modules.sampling.attract import HaddockModule as AttractModule

from . import GOLDEN_DATA, has_attract


class MockPreviousIO:
    """Mock the previous IO module."""

    def __init__(self, path):
        self.path = path

    def retrieve_models(self):
        """Mock the retrieval of some models"""
        shutil.copy(
            Path(GOLDEN_DATA, "prot.pdb"),
            Path(self.path, "prot.pdb"),
        )

        shutil.copy(
            Path(GOLDEN_DATA, "rna.pdb"),
            Path(self.path, "rna.pdb"),
        )

        return [
            PDBFile(file_name="prot.pdb", path=self.path),
            PDBFile(file_name="rna.pdb", path=self.path),
        ]

    def output(self) -> None:
        """Mock the output"""
        return None


@pytest.fixture(name="attract_module")
def fixture_attract_module():
    """Initialize the attract module"""
    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)
        yield AttractModule(
            order=0, path=Path("."), initial_params=DEFAULT_ATTRACT_CONFIG
        )


@has_attract
@pytest.mark.skip(reason="work-in-progress")
def test_attract(attract_module):
    """Integration test for the attract module"""

    attract_module.previous_io = MockPreviousIO(path=attract_module.path)

    # attract_module.attract_tools = ???
    # attract_module.nalib = ???

    attract_module.run()

    assert len(attract_module.output_models) > 0
