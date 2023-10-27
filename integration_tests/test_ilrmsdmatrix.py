import os
import tempfile
from pathlib import Path
import shutil
import pytest
import pytest_mock

from haddock.libs.libontology import PDBFile

from haddock.modules.analysis.ilrmsdmatrix import DEFAULT_CONFIG as DEFAULT_ILRMSD_CONFIG
from haddock.modules.analysis.ilrmsdmatrix import HaddockModule as IlrmsdmatrixModule
DATA_DIR = Path(Path(__file__).parent.parent / "tests" / "golden_data")


@pytest.fixture
def ilrmsdmatrix_module():
    with tempfile.TemporaryDirectory() as tmpdir:
        ilrmsdmatrix = IlrmsdmatrixModule(
            order=0, path=tmpdir, initial_params=DEFAULT_ILRMSD_CONFIG
        )
        ilrmsdmatrix.params["receptor_chain"] = "R"
        ilrmsdmatrix.params["ligand_chains"] = ["S"]
        yield ilrmsdmatrix


class MockPreviousIO():
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = False):
        shutil.copy(Path(DATA_DIR, "protprot_1bkd_1.pdb"), Path(".", "protprot_complex_1.pdb"))
        shutil.copy(Path(DATA_DIR, "protprot_1bkd_2.pdb"), Path(".", "protprot_complex_2.pdb"))
        
        model_list = [
            PDBFile(file_name="protprot_complex_1.pdb", path="."),
            PDBFile(file_name="protprot_complex_2.pdb", path="."),
        ]

        return model_list


def test_ilrmsdmatrix_default(ilrmsdmatrix_module, mocker):
    """Test the topoaa module."""
    ilrmsdmatrix_module.previous_io = MockPreviousIO(path=ilrmsdmatrix_module.path)
    mocker.patch("haddock.modules.BaseHaddockModule.export_io_models", return_value = None)
    ilrmsdmatrix_module.run()
    # expected paths
    exp_ilrmsd_matrix = Path(ilrmsdmatrix_module.path, "ilrmsd.matrix")
    exp_contacts_file = Path(ilrmsdmatrix_module.path, "receptor_contacts.con")
    assert exp_ilrmsd_matrix.exists(), "ilrmsd.matrix does not exist"
    assert exp_contacts_file.exists(), "receptor_contacts.con does not exist"
    # open files and check content
    with open(exp_ilrmsd_matrix) as f:
        assert f.readline() == "1 2 20.020\n"
    with open(exp_contacts_file) as f:
        lines = f.readlines()
        assert lines[0] == "R 12 13 15 17 21 25 30 31 32 33 34 35 37 38 39 40 41 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 73 95 102 103 105\n"
        assert lines[1] == "S 809 810 814 822 824 825 826 828 829 832 833 836 869 872 875 876 879 880 881 882 884 908 909 910 911 912 913 929 930 931 932 934 935 936 938 939 940 942 943 944 945 963 1002 1003 1006 1007 1010 1019\n"
