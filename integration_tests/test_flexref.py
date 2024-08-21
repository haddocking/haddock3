import shutil
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import Format, PDBFile, Persistent
from haddock.modules.refinement.flexref import (
    DEFAULT_CONFIG as DEFAULT_FLEXREF_CONFIG,
)
from haddock.modules.refinement.flexref import HaddockModule as FlexrefModule

from . import golden_data


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
            Path(golden_data, "protprot_complex.pdb"),
            Path(self.path, "protprot_complex.pdb"),
        )
        shutil.copy(
            Path(golden_data, "e2aP_1F3G_haddock.psf"),
            Path(self.path, "e2aP_1F3G_haddock.psf"),
        )

        shutil.copy(
            Path(golden_data, "hpr_ensemble_1_haddock.psf"),
            Path(self.path, "hpr_ensemble_1_haddock.psf"),
        )

        model_list = [
            PDBFile(
                file_name="protprot_complex.pdb",
                path=self.path,
                topology=(
                    Persistent(
                        file_name="hpr_ensemble_1_haddock.psf",
                        path=self.path,
                        file_type=Format.TOPOLOGY,
                    ),
                    Persistent(
                        file_name="e2aP_1F3G_haddock.psf",
                        path=self.path,
                        file_type=Format.TOPOLOGY,
                    ),
                ),
            ),
        ]

        return model_list

    def output(self):
        return None


def test_flexref_defaults(flexref_module, calc_fnat):

    flexref_module.previous_io = MockPreviousIO(path=flexref_module.path)

    flexref_module.run()

    assert Path(flexref_module.path, "flexref_1.pdb").exists()
    assert Path(flexref_module.path, "flexref_1.out.gz").exists()

    fnat = calc_fnat(
        model=Path(flexref_module.path, "flexref_1.pdb"),
        native=Path(golden_data, "protprot_complex.pdb"),
    )

    assert fnat == pytest.approx(0.9, abs=0.1)


def test_flexref_fle(flexref_module, calc_fnat):

    flexref_module.previous_io = MockPreviousIO(path=flexref_module.path)

    flexref_module.params["nfle "] = 1
    flexref_module.params["fle_sta_1 "] = 141
    flexref_module.params["fle_end_1 "] = 146
    flexref_module.params["fle_seg_1  "] = "A"

    flexref_module.run()

    assert Path(flexref_module.path, "flexref_1.pdb").exists()
    assert Path(flexref_module.path, "flexref_1.out.gz").exists()

    fnat = calc_fnat(
        model=Path(flexref_module.path, "flexref_1.pdb"),
        native=Path(golden_data, "protprot_complex.pdb"),
    )

    assert fnat == pytest.approx(0.9, abs=0.1)


def test_flexref_mutliple_fle(flexref_module, calc_fnat):

    flexref_module.previous_io = MockPreviousIO(path=flexref_module.path)

    flexref_module.params["nfle "] = 2
    flexref_module.params["fle_sta_1 "] = 141
    flexref_module.params["fle_end_1 "] = 146
    flexref_module.params["fle_seg_1  "] = "A"

    flexref_module.params["fle_sta_2 "] = 66
    flexref_module.params["fle_end_2 "] = 71
    flexref_module.params["fle_seg_2  "] = "B"

    flexref_module.run()

    assert Path(flexref_module.path, "flexref_1.pdb").exists()
    assert Path(flexref_module.path, "flexref_1.out.gz").exists()

    fnat = calc_fnat(
        model=Path(flexref_module.path, "flexref_1.pdb"),
        native=Path(golden_data, "protprot_complex.pdb"),
    )

    assert fnat == pytest.approx(0.0, abs=0.1)
