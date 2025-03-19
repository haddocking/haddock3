import shutil
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import Format, PDBFile, Persistent
from haddock.modules.refinement.emref import DEFAULT_CONFIG as DEFAULT_EMREF_CONFIG
from haddock.modules.refinement.emref import HaddockModule as FlexrefModule

from integration_tests import GOLDEN_DATA


@pytest.fixture
def emref_module():
    with tempfile.TemporaryDirectory() as tmpdir:
        emref = FlexrefModule(
            order=0, path=Path(tmpdir), initial_params=DEFAULT_EMREF_CONFIG
        )
        yield emref


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


def test_emref_defaults(emref_module, calc_fnat):

    emref_module.previous_io = MockPreviousIO(path=emref_module.path)
    emref_module.params["debug"] = True

    emref_module.run()

    assert Path(emref_module.path, "emref_1.pdb").exists()
    assert Path(emref_module.path, "emref_1.out.gz").exists()

    fnat = calc_fnat(
        model=Path(emref_module.path, "emref_1.pdb"),
        native=Path(GOLDEN_DATA, "2oob.pdb"),
    )

    assert fnat == pytest.approx(0.95, abs=0.05)


def test_emref_fle(emref_module, calc_fnat):

    emref_module.previous_io = MockPreviousIO(path=emref_module.path)

    emref_module.params["nfle "] = 1

    emref_module.params["fle_sta_1 "] = 66
    emref_module.params["fle_end_1 "] = 77
    emref_module.params["fle_seg_1  "] = "B"

    emref_module.params["debug"] = True

    emref_module.run()

    assert Path(emref_module.path, "emref_1.pdb").exists()
    assert Path(emref_module.path, "emref_1.out.gz").exists()

    fnat = calc_fnat(
        model=Path(emref_module.path, "emref_1.pdb"),
        native=Path(GOLDEN_DATA, "2oob.pdb"),
    )

    assert fnat == pytest.approx(0.95, abs=0.05)


def test_emref_mutliple_fle(emref_module, calc_fnat):

    emref_module.previous_io = MockPreviousIO(path=emref_module.path)

    emref_module.params["nfle "] = 2

    emref_module.params["fle_sta_1 "] = 66
    emref_module.params["fle_end_1 "] = 77
    emref_module.params["fle_seg_1  "] = "B"

    emref_module.params["fle_sta_2 "] = 41
    emref_module.params["fle_end_2 "] = 47
    emref_module.params["fle_seg_2  "] = "B"

    emref_module.params["debug"] = True

    emref_module.run()

    assert Path(emref_module.path, "emref_1.pdb").exists()
    assert Path(emref_module.path, "emref_1.out.gz").exists()

    fnat = calc_fnat(
        model=Path(emref_module.path, "emref_1.pdb"),
        native=Path(GOLDEN_DATA, "2oob.pdb"),
    )

    assert fnat == pytest.approx(0.95, abs=0.05)
