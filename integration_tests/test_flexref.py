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


def test_flexref_defaults(flexref_module, run_dockq_analysis):

    flexref_module.previous_io = MockPreviousIO(Path(golden_data))

    flexref_module.run()

    assert Path(flexref_module.path, "flexref_1.pdb").exists()
    assert Path(flexref_module.path, "flexref_1.out.gz").exists()

    dockq_analysis = run_dockq_analysis(
        model=Path(flexref_module.path, "flexref_1.pdb"),
        native=Path(golden_data, "protprot_complex.pdb"),
    )
    # assert dockq_analysis["DockQ"] == pytest.approx(0.65, abs=0.05)
    # assert dockq_analysis["iRMSD"] == pytest.approx(1.4, abs=0.5)
    assert dockq_analysis["fnat"] == pytest.approx(0.9, abs=0.1)
    # assert dockq_analysis["LRMSD"] == pytest.approx(7.1, abs=0.1)


def test_flexref_fle(flexref_module, run_dockq_analysis):

    flexref_module.previous_io = MockPreviousIO(Path(golden_data))

    flexref_module.params["nfle "] = 1
    flexref_module.params["fle_sta_1 "] = 141
    flexref_module.params["fle_end_1 "] = 146
    flexref_module.params["fle_seg_1  "] = "A"

    flexref_module.run()

    assert Path(flexref_module.path, "flexref_1.pdb").exists()
    assert Path(flexref_module.path, "flexref_1.out.gz").exists()

    dockq_analysis = run_dockq_analysis(
        model=Path(flexref_module.path, "flexref_1.pdb"),
        native=Path(golden_data, "protprot_complex.pdb"),
    )
    # assert dockq_analysis["DockQ"] == pytest.approx(0.65, abs=0.05)
    # assert dockq_analysis["iRMSD"] == pytest.approx(1.4, abs=0.5)
    assert dockq_analysis["fnat"] == pytest.approx(0.9, abs=0.1)
    # assert dockq_analysis["LRMSD"] == pytest.approx(6.0, abs=1.0)


def test_flexref_mutliple_fle(flexref_module, run_dockq_analysis):

    flexref_module.previous_io = MockPreviousIO(Path(golden_data))

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

    dockq_analysis = run_dockq_analysis(
        model=Path(flexref_module.path, "flexref_1.pdb"),
        native=Path(golden_data, "protprot_complex.pdb"),
    )
    # assert dockq_analysis["DockQ"] == pytest.approx(0.65, abs=0.05)
    # assert dockq_analysis["iRMSD"] == pytest.approx(1.4, abs=0.5)
    assert dockq_analysis["fnat"] == pytest.approx(0.9, abs=0.1)
    # assert dockq_analysis["LRMSD"] == pytest.approx(7.1, abs=0.1)
