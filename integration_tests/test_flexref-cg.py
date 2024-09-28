import shutil
import tempfile
import gzip
from pathlib import Path

import pytest

from haddock.libs.libontology import Format, PDBFile, Persistent
from haddock.modules.refinement.flexref import (
    DEFAULT_CONFIG as DEFAULT_FLEXREF_CONFIG,
)
from haddock.modules.refinement.flexref import HaddockModule as FlexrefModule

from integration_tests import GOLDEN_DATA


@pytest.fixture
def flexref_module():
    with tempfile.TemporaryDirectory() as tmpdir:
        flexref = FlexrefModule(
            order=0, path=Path(tmpdir), initial_params=DEFAULT_FLEXREF_CONFIG
        )
        yield flexref


class MockPreviousIO:
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, crossdock: bool = False):
        shutil.copy(
            Path(GOLDEN_DATA, "e2a-hpr-rigid-cg.pdb"),
            Path(self.path, "e2a-hpr-rigid-cg.pdb"),
        )
        shutil.copy(
            Path(GOLDEN_DATA, "e2a_haddock_cg.psf"),
            Path(self.path, "e2a_haddock_cg.psf"),
        )

        shutil.copy(
            Path(GOLDEN_DATA, "hpr_haddock_cg.psf"),
            Path(self.path, "hpr_haddock_cg.psf"),
        )

        model = PDBFile(
            file_name="e2a-hpr-rigid-cg.pdb",
            path=self.path,
            topology=(
                Persistent(
                    file_name="e2a_haddock_cg.psf",
                    path=self.path,
                    file_type=Format.TOPOLOGY,
                ),
                Persistent(
                    file_name="hpr_haddock_cg.psf",
                    path=self.path,
                    file_type=Format.TOPOLOGY,
                ),
            ),
        )

        model.seed = 42  # type: ignore

        return [model]

    def output(self):
        return None


def test_flexref_defaults(flexref_module, calc_fnat):

    flexref_module.previous_io = MockPreviousIO(path=flexref_module.path)
    flexref_module.params["epsilon"] = 10
    flexref_module.params["debug"] = True

    flexref_module.run()

    expected_inp = Path(flexref_module.path, "flexref_1.inp")
    expected_pdb = Path(flexref_module.path, "flexref_1.pdb")
    expected_gz = Path(flexref_module.path, "flexref_1.out.gz")
    
    assert expected_inp.exists(), f"{expected_inp} does not exist"
    assert expected_gz.exists(), f"{expected_gz} does not exist"
    assert expected_pdb.exists(), f"{expected_pdb} does not exist"

    #check for correct definition of non-bonded options for CG
    file_content = gzip.open(Path(flexref_module.path, expected_gz), 'rt').read()
    assert 'CTONNB=  12.000' in file_content

