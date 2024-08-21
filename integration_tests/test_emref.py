import shutil
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import Format, PDBFile, Persistent
from haddock.modules.refinement.emref import DEFAULT_CONFIG as DEFAULT_EMREF_CONFIG
from haddock.modules.refinement.emref import HaddockModule as FlexrefModule

from . import golden_data


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
            Path(golden_data, "protprot_complex.pdb"),
            Path(".", "protprot_complex.pdb"),
        )
        shutil.copy(
            Path(golden_data, "e2aP_1F3G_haddock.psf"),
            Path(".", "e2aP_1F3G_haddock.psf"),
        )

        shutil.copy(
            Path(golden_data, "hpr_ensemble_1_haddock.psf"),
            Path(".", "hpr_ensemble_1_haddock.psf"),
        )

        model_list = [
            PDBFile(
                file_name="protprot_complex.pdb",
                path=".",
                topology=(
                    Persistent(
                        file_name="hpr_ensemble_1_haddock.psf",
                        path=".",
                        file_type=Format.TOPOLOGY,
                    ),
                    Persistent(
                        file_name="e2aP_1F3G_haddock.psf",
                        path=".",
                        file_type=Format.TOPOLOGY,
                    ),
                ),
            ),
        ]

        return model_list

    def output(self):
        return None


def test_emref_defaults(emref_module, calc_fnat):

    emref_module.previous_io = MockPreviousIO(Path(golden_data))

    emref_module.run()

    assert Path(emref_module.path, "emref_1.pdb").exists()
    assert Path(emref_module.path, "emref_1.out.gz").exists()

    fnat = calc_fnat(
        model=Path(emref_module.path, "emref_1.pdb"),
        native=Path(golden_data, "protprot_complex.pdb"),
    )

    assert fnat == pytest.approx(0.90, abs=0.05)


def test_emref_fle(emref_module, calc_fnat):

    emref_module.previous_io = MockPreviousIO(Path(golden_data))

    emref_module.params["nfle "] = 1
    emref_module.params["fle_sta_1 "] = 141
    emref_module.params["fle_end_1 "] = 146
    emref_module.params["fle_seg_1  "] = "A"

    emref_module.run()

    assert Path(emref_module.path, "emref_1.pdb").exists()
    assert Path(emref_module.path, "emref_1.out.gz").exists()

    fnat = calc_fnat(
        model=Path(emref_module.path, "emref_1.pdb"),
        native=Path(golden_data, "protprot_complex.pdb"),
    )

    assert fnat == pytest.approx(0.90, abs=0.05)


def test_emref_mutliple_fle(emref_module, calc_fnat):

    emref_module.previous_io = MockPreviousIO(Path(golden_data))

    emref_module.params["nfle "] = 2
    emref_module.params["fle_sta_1 "] = 141
    emref_module.params["fle_end_1 "] = 146
    emref_module.params["fle_seg_1  "] = "A"

    emref_module.params["fle_sta_2 "] = 66
    emref_module.params["fle_end_2 "] = 71
    emref_module.params["fle_seg_2  "] = "B"

    emref_module.run()

    assert Path(emref_module.path, "emref_1.pdb").exists()
    assert Path(emref_module.path, "emref_1.out.gz").exists()

    fnat = calc_fnat(
        model=Path(emref_module.path, "emref_1.pdb"),
        native=Path(golden_data, "protprot_complex.pdb"),
    )

    assert fnat == pytest.approx(0.90, abs=0.05)
