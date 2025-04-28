from haddock.modules.analysis.saxsscore import HaddockModule as SaxsModule

import shutil
from pathlib import Path
import pytest
from haddock.libs.libontology import PDBFile
from tests import golden_data

from haddock.modules.analysis.saxsscore import (
    DEFAULT_CONFIG as DEFAULT_SAXSSCORE_CONFIG,
)
from haddock.modules.analysis.saxsscore.saxs import run_crysol


class MockPreviousIO:
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = False):
        shutil.copy(
            Path(golden_data, "e2aP_1F3G_haddock.pdb"),
            Path(".", "e2aP_1F3G_haddock.pdb"),
        )
        model_list = [
            PDBFile(
                file_name="e2aP_1F3G_haddock.pdb",
                path=".",
                score=0.0,
            ),
        ]

        return model_list

    def output(self):
        return None


@pytest.fixture
def saxsscore():
    return SaxsModule(
        order=0,
        path=Path("."),
        init_params=DEFAULT_SAXSSCORE_CONFIG,
    )


@pytest.fixture
def saxs_data():
    yield Path(golden_data, "saxs.dat")


ATSAS_PATH = "/home/rodrigo/software/ATSAS-3.2.1-1"


def test_saxsscore(saxsscore, saxs_data):

    saxsscore.previous_io = MockPreviousIO(path=saxsscore.path)
    saxsscore.params["saxs_fname"] = saxs_data
    saxsscore.params["atsas_path"] = ATSAS_PATH

    saxsscore.run()

    assert Path("saxsscore.tsv").exists()


def test_crysol_is_executable(saxs_data):
    pdb_f = Path(golden_data, "hpr_ensemble_1_haddock.pdb")
    run_crysol(
        atsas_path=ATSAS_PATH,
        input_f=pdb_f,
        saxs_data=saxs_data,
        lm=20,
        ns=20,
    )

    assert Path(f"{Path(pdb_f.stem)}_saxs.fit").exists()
