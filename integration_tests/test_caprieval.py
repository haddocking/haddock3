import math
import shutil
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.caprieval import (
    DEFAULT_CONFIG as DEFAULT_CAPRIEVAL_CONFIG,
)
from haddock.modules.analysis.caprieval import HaddockModule as CaprievalModule
from tests import golden_data


@pytest.fixture
def caprieval_module():
    with tempfile.TemporaryDirectory() as tmpdir:
        yield CaprievalModule(
            order=0,
            path=Path(tmpdir),
            init_params=DEFAULT_CAPRIEVAL_CONFIG,
        )


@pytest.fixture
def model_list():
    return [
        PDBFile(
            file_name="protprot_complex_1.pdb",
            path=".",
            unw_energies={"energy_term": 0.0},
        ),
        PDBFile(
            file_name="protprot_complex_2.pdb",
            path=".",
            unw_energies={"energy_term": 0.0},
        ),
    ]


class MockPreviousIO:
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = True):
        shutil.copy(
            Path(golden_data, "protprot_complex_1.pdb"),
            Path(".", "protprot_complex_1.pdb"),
        )
        shutil.copy(
            Path(golden_data, "protprot_complex_2.pdb"),
            Path(".", "protprot_complex_2.pdb"),
        )
        model_list = [
            PDBFile(
                file_name="protprot_complex_1.pdb",
                path=".",
                unw_energies={"energy_term": 0.0},
            ),
            PDBFile(
                file_name="protprot_complex_2.pdb",
                path=".",
                unw_energies={"energy_term": 0.0},
            ),
        ]

        return model_list

    def output(self):
        return None


def evaluate_caprieval_execution(module: CaprievalModule, model_list):
    """Helper function to check if `caprieval` executed properly."""

    # Check if the files were written
    expected_files = [
        "blosum62.izone",
        "blosum62_A.aln",
        "blosum62_B.aln",
        "capri_clt.tsv",
        "capri_ss.tsv",
    ]

    for file in expected_files:
        assert Path(module.path, file).exists(), f"{file} does not exist"
        assert Path(module.path, file).stat().st_size > 0, f"{file} is empty"

    # Check if the `output_models` are the same as the `input_models`
    #  they will change locations, so check the filenames
    assert module.output_models[0].file_name == model_list[0].file_name
    assert module.output_models[1].file_name == model_list[1].file_name

    # The models do not hold the capri metrics, so check the output files to see if they were written
    with open(Path(module.path, "capri_ss.tsv")) as f:
        lines = f.readlines()

    # Check the values
    assert len(lines) == 3, "There should be 3 lines in the capri_ss.tsv"

    # Check the header
    # model   md5     caprieval_rank  score   irmsd   fnat    lrmsd   ilrmsd  dockq   cluster_id      cluster_ranking model-cluster_ranking
    expected_colnames = [
        "model",
        "md5",
        "caprieval_rank",
        "score",
        "irmsd",
        "fnat",
        "lrmsd",
        "ilrmsd",
        "dockq",
        "cluster_ranking",
        "model-cluster_ranking",
        "energy_term",
    ]
    header = lines[0].strip().split("\t")
    for col_name in expected_colnames:
        assert col_name in header, f"{col_name} not found in the header"

    # Check the values
    data = lines[1:]
    expected_data = [
        {
            "model": "",
            "md5": "-",
            "caprieval_rank": 1,
            "score": float("nan"),
            "irmsd": 0.000,
            "fnat": 1.000,
            "lrmsd": 0.000,
            "ilrmsd": 0.000,
            "dockq": 1.000,
            "cluster_id": "-",
            "cluster_ranking": "-",
            "model-cluster_ranking": "-",
            "energy_term": 0.000,
        },
        {
            "model": "",
            "md5": "-",
            "caprieval_rank": 2,
            "score": float("nan"),
            "irmsd": 8.327,
            "fnat": 0.050,
            "lrmsd": 20.937,
            "ilrmsd": 18.248,
            "dockq": 0.074,
            "cluster_id": "-",
            "cluster_ranking": "-",
            "model-cluster_ranking": "-",
            "energy_term": 0.000,
        },
    ]
    oberseved_data = []
    for line in data:
        values = line.strip().split("\t")
        data_dict = {
            "model": "",  # don't check this, it's a path and will change
            "md5": str(values[1]),
            "caprieval_rank": int(values[2]),
            "score": float(values[3]),
            "irmsd": float(values[4]),
            "fnat": float(values[5]),
            "lrmsd": float(values[6]),
            "ilrmsd": float(values[7]),
            "dockq": float(values[8]),
            "cluster_id": str(values[9]),
            "cluster_ranking": str(values[10]),
            "model-cluster_ranking": str(values[11]),
            "energy_term": float(values[12]),
        }
        oberseved_data.append(data_dict)

    for k, v in zip(expected_data, oberseved_data):
        for key in k:
            v1 = k[key]
            v2 = v[key]

            # Check the type
            assert type(v1) == type(v2), f"Type mismatch for {key}"

            if isinstance(v1, float):
                if math.isnan(v1):
                    assert math.isnan(v2), f"Value mismatch for {key}"
                else:
                    assert (
                        pytest.approx(v1, rel=1e-3) == v2
                    ), f"Value mismatch for {key}"
            else:
                assert v1 == v2, f"Value mismatch for {key}"


def test_caprieval_default(caprieval_module, model_list):

    caprieval_module.previous_io = MockPreviousIO(path=caprieval_module.path)
    caprieval_module.run()

    evaluate_caprieval_execution(caprieval_module, model_list)


def test_caprieval_less_io(caprieval_module, model_list):
    caprieval_module.previous_io = MockPreviousIO(path=caprieval_module.path)
    caprieval_module.params["less_io"] = True

    caprieval_module.run()

    evaluate_caprieval_execution(caprieval_module, model_list)
