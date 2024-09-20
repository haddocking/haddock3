import shutil
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import Format, PDBFile, Persistent
from haddock.modules.refinement.mdref import DEFAULT_CONFIG as DEFAULT_MDREF_CONFIG
from haddock.modules.refinement.mdref import HaddockModule as FlexrefModule

from integration_tests import GOLDEN_DATA


@pytest.fixture
def mdref_module():
    with tempfile.TemporaryDirectory() as tmpdir:
        mdref = FlexrefModule(
            order=0, path=Path(tmpdir), initial_params=DEFAULT_MDREF_CONFIG
        )
        yield mdref


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


def test_mdref_defaults(mdref_module, calc_fnat):

    mdref_module.previous_io = MockPreviousIO(path=mdref_module.path)
    mdref_module.params["debug"] = True
    mdref_module.run()

    assert Path(mdref_module.path, "mdref_1.pdb").exists()
    assert Path(mdref_module.path, "mdref_1.out.gz").exists()

    fnat = calc_fnat(
        model=Path(mdref_module.path, "mdref_1.pdb"),
        native=Path(GOLDEN_DATA, "2oob.pdb"),
    )

    assert fnat == pytest.approx(0.9, abs=0.1)


def test_mdref_fle(mdref_module, calc_fnat):

    mdref_module.previous_io = MockPreviousIO(path=mdref_module.path)

    mdref_module.params["nfle "] = 1

    mdref_module.params["fle_sta_1 "] = 66
    mdref_module.params["fle_end_1 "] = 77
    mdref_module.params["fle_seg_1  "] = "B"

    mdref_module.params["debug"] = True

    mdref_module.run()

    assert Path(mdref_module.path, "mdref_1.pdb").exists()
    assert Path(mdref_module.path, "mdref_1.out.gz").exists()

    fnat = calc_fnat(
        model=Path(mdref_module.path, "mdref_1.pdb"),
        native=Path(GOLDEN_DATA, "2oob.pdb"),
    )

    assert fnat == pytest.approx(0.8, abs=0.1)


def test_mdref_mutliple_fle(mdref_module, calc_fnat):

    mdref_module.previous_io = MockPreviousIO(path=mdref_module.path)

    mdref_module.params["nfle "] = 2

    mdref_module.params["fle_sta_1 "] = 66
    mdref_module.params["fle_end_1 "] = 77
    mdref_module.params["fle_seg_1  "] = "B"

    mdref_module.params["fle_sta_2 "] = 41
    mdref_module.params["fle_end_2 "] = 47
    mdref_module.params["fle_seg_2  "] = "B"

    mdref_module.params["debug"] = True

    mdref_module.run()

    assert Path(mdref_module.path, "mdref_1.pdb").exists()
    assert Path(mdref_module.path, "mdref_1.out.gz").exists()

    fnat = calc_fnat(
        model=Path(mdref_module.path, "mdref_1.pdb"),
        native=Path(GOLDEN_DATA, "2oob.pdb"),
    )

    assert fnat == pytest.approx(0.8, abs=0.1)
