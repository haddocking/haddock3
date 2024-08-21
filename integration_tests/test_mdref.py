import shutil
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import Format, PDBFile, Persistent
from haddock.modules.refinement.mdref import DEFAULT_CONFIG as DEFAULT_MDREF_CONFIG
from haddock.modules.refinement.mdref import HaddockModule as FlexrefModule

from . import golden_data


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


def test_mdref_defaults(mdref_module, run_dockq_analysis):

    mdref_module.previous_io = MockPreviousIO(Path(golden_data))

    mdref_module.run()

    assert Path(mdref_module.path, "mdref_1.pdb").exists()
    assert Path(mdref_module.path, "mdref_1.out.gz").exists()

    dockq_analysis = run_dockq_analysis(
        model=Path(mdref_module.path, "mdref_1.pdb"),
        native=Path(golden_data, "protprot_complex.pdb"),
    )
    # assert dockq_analysis["DockQ"] == pytest.approx(0.85, abs=0.05)
    # assert dockq_analysis["iRMSD"] == pytest.approx(0.75, abs=0.05)
    assert dockq_analysis["fnat"] == pytest.approx(0.9, abs=0.1)
    # assert dockq_analysis["LRMSD"] == pytest.approx(0.85, abs=0.05)


def test_mdref_fle(mdref_module, run_dockq_analysis):

    mdref_module.previous_io = MockPreviousIO(Path(golden_data))

    mdref_module.params["nfle "] = 1
    mdref_module.params["fle_sta_1 "] = 141
    mdref_module.params["fle_end_1 "] = 146
    mdref_module.params["fle_seg_1  "] = "A"

    mdref_module.run()

    assert Path(mdref_module.path, "mdref_1.pdb").exists()
    assert Path(mdref_module.path, "mdref_1.out.gz").exists()

    dockq_analysis = run_dockq_analysis(
        model=Path(mdref_module.path, "mdref_1.pdb"),
        native=Path(golden_data, "protprot_complex.pdb"),
    )
    # assert dockq_analysis["DockQ"] == pytest.approx(0.85, abs=0.05)
    # assert dockq_analysis["iRMSD"] == pytest.approx(0.75, abs=0.05)
    assert dockq_analysis["fnat"] == pytest.approx(0.75, abs=0.1)
    # assert dockq_analysis["LRMSD"] == pytest.approx(0.85, abs=0.05)


def test_mdref_mutliple_fle(mdref_module, run_dockq_analysis):

    mdref_module.previous_io = MockPreviousIO(Path(golden_data))

    mdref_module.params["nfle "] = 2
    mdref_module.params["fle_sta_1 "] = 141
    mdref_module.params["fle_end_1 "] = 146
    mdref_module.params["fle_seg_1  "] = "A"

    mdref_module.params["fle_sta_2 "] = 66
    mdref_module.params["fle_end_2 "] = 71
    mdref_module.params["fle_seg_2  "] = "B"

    mdref_module.run()

    assert Path(mdref_module.path, "mdref_1.pdb").exists()
    assert Path(mdref_module.path, "mdref_1.out.gz").exists()

    dockq_analysis = run_dockq_analysis(
        model=Path(mdref_module.path, "mdref_1.pdb"),
        native=Path(golden_data, "protprot_complex.pdb"),
    )
    # assert dockq_analysis["DockQ"] == pytest.approx(0.85, abs=0.05)
    # assert dockq_analysis["iRMSD"] == pytest.approx(0.75, abs=0.05)
    assert dockq_analysis["fnat"] == pytest.approx(0.9, abs=0.1)
    # assert dockq_analysis["LRMSD"] == pytest.approx(0.85, abs=0.05)
