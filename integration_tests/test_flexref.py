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
            Path(GOLDEN_DATA, "2oob.pdb"),
            Path(self.path, "2oob.pdb"),
        )
        shutil.copy(
            Path(GOLDEN_DATA, "2oob_A.psf"),
            Path(self.path, "2oob_A.psf"),
        )

        shutil.copy(
            Path(GOLDEN_DATA, "2oob_B.psf"),
            Path(self.path, "2oob_B.psf"),
        )

        model = PDBFile(
            file_name="2oob.pdb",
            path=self.path,
            topology=(
                Persistent(
                    file_name="2oob_A.psf",
                    path=self.path,
                    file_type=Format.TOPOLOGY,
                ),
                Persistent(
                    file_name="2oob_B.psf",
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

    flexref_module.run()

    assert Path(flexref_module.path, "flexref_1.pdb").exists()

    fnat = calc_fnat(
        model=Path(flexref_module.path, "flexref_1.pdb"),
        native=Path(GOLDEN_DATA, "2oob.pdb"),
    )

    assert fnat == pytest.approx(0.9, abs=0.1)


def test_flexref_fle(flexref_module, calc_fnat):

    flexref_module.previous_io = MockPreviousIO(path=flexref_module.path)

    flexref_module.params["nfle "] = 1

    flexref_module.params["fle_sta_1 "] = 66
    flexref_module.params["fle_end_1 "] = 77
    flexref_module.params["fle_seg_1  "] = "B"
    flexref_module.params["log_level"] = "verbose"
    flexref_module.params["debug"] = True

    flexref_module.run()

    assert Path(flexref_module.path, "flexref_1.pdb").exists()
    assert Path(flexref_module.path, "flexref_1.out.gz").exists()
    file_content = gzip.open(Path(flexref_module.path, "flexref_1.out.gz"), 'rt').read()
    assert '$SAPROTOCOL.TADFACTOR set to    4.00000' in file_content

    fnat = calc_fnat(
        model=Path(flexref_module.path, "flexref_1.pdb"),
        native=Path(GOLDEN_DATA, "2oob.pdb"),
    )

    assert fnat == pytest.approx(0.9, abs=0.1)


def test_flexref_mutliple_fle(flexref_module, calc_fnat):

    flexref_module.previous_io = MockPreviousIO(path=flexref_module.path)

    flexref_module.params["nfle "] = 2

    flexref_module.params["fle_sta_1 "] = 66
    flexref_module.params["fle_end_1 "] = 77
    flexref_module.params["fle_seg_1  "] = "B"

    flexref_module.params["fle_sta_2 "] = 41
    flexref_module.params["fle_end_2 "] = 47
    flexref_module.params["fle_seg_2  "] = "B"

    flexref_module.params["debug"] = True

    flexref_module.run()

    assert Path(flexref_module.path, "flexref_1.pdb").exists()
    assert Path(flexref_module.path, "flexref_1.out.gz").exists()

    fnat = calc_fnat(
        model=Path(flexref_module.path, "flexref_1.pdb"),
        native=Path(GOLDEN_DATA, "2oob.pdb"),
    )

    assert fnat == pytest.approx(0.9, abs=0.1)
