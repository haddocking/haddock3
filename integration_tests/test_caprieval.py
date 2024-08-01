import math
import shutil
import tempfile
from pathlib import Path
from typing import Union

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


@pytest.fixture
def expected_clt_data() -> list[dict[str, Union[int, str, float]]]:
    return [
        {
            "cluster_rank": "-",
            "cluster_id": "-",
            "n": 2,
            "under_eval": "yes",
            "score": float("nan"),
            "score_std": float("nan"),
            "irmsd": 4.163,
            "irmsd_std": 4.163,
            "fnat": 0.525,
            "fnat_std": 0.475,
            "lrmsd": 10.469,
            "lrmsd_std": 10.469,
            "dockq": 0.537,
            "dockq_std": 0.463,
            "ilrmsd": 9.124,
            "ilrmsd_std": 9.124,
            "air": float("nan"),
            "air_std": float("nan"),
            "bsa": float("nan"),
            "bsa_std": float("nan"),
            "desolv": float("nan"),
            "desolv_std": float("nan"),
            "elec": float("nan"),
            "elec_std": float("nan"),
            "total": float("nan"),
            "total_std": float("nan"),
            "vdw": float("nan"),
            "vdw_std": float("nan"),
            "caprieval_rank": 1,
        }
    ]


@pytest.fixture
def expected_ss_data() -> list[dict[str, Union[int, str, float]]]:
    return [
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


def _cast_float_str_int(v: Union[int, str, float]) -> Union[int, str, float]:
    """Helper function to cast a value string to a float, int or str."""
    try:
        return int(v)
    except ValueError:
        try:
            return float(v)
        except ValueError:
            return v


def _compare_polymorphic_data(
    expected_data: list[dict[str, Union[int, str, float]]],
    oberseved_data: list[dict[str, Union[int, str, float]]],
):
    """Helper function to compare a list of dictionaries with polymorphic values."""
    for k, v in zip(expected_data, oberseved_data):
        for key in k:
            v1 = k[key]
            v2 = v[key]

            # Check the type match
            assert type(v1) == type(v2), f"Type mismatch for {key}"

            # Check float
            if isinstance(v1, float):
                if math.isnan(v1):
                    assert isinstance(v2, float) and math.isnan(
                        v2
                    ), f"Value mismatch for {key}"
                else:
                    assert (
                        isinstance(v2, (int, float))
                        and pytest.approx(v1, rel=1e-3) == v2
                    ), f"Value mismatch for {key}"

            # Check int
            elif isinstance(v1, int):
                assert v1 == v2, f"Value mismatch for {key}"

            # Check str
            elif isinstance(v1, str):
                assert v1 == v2, f"Value mismatch for {key}"

            # Value is not float, int or str
            else:
                raise ValueError(f"Unknown type for {key}")


def _check_capri_ss_tsv(
    capri_file: str, expected_data: list[dict[str, Union[int, str, float]]]
):
    """Helper function to check the content of the capri_ss.tsv file."""
    with open(capri_file) as f:
        lines = f.readlines()

    # Check the header
    expected_header_cols = list(expected_data[0].keys())
    observed_header_cols = lines[0].strip().split("\t")

    # Check if all they have the same lenght
    assert len(observed_header_cols) == len(expected_header_cols), "Header mismatch"

    for col_name in expected_header_cols:
        assert col_name in observed_header_cols, f"{col_name} not found in the header"

    oberseved_data: list[dict[str, Union[int, str, float]]] = []
    data = lines[1:]
    for line in data:
        values = line.strip().split("\t")

        # Check there is one value for each column
        assert len(values) == len(expected_header_cols), "Values mismatch"

        data_dict = {}
        for h, v in zip(expected_header_cols, values):
            data_dict[h] = _cast_float_str_int(v)

        oberseved_data.append(data_dict)

    # Cannot compare the names of the models, since the observed will be a random string
    [d.pop("model") for d in expected_data], [d.pop("model") for d in oberseved_data]

    _compare_polymorphic_data(expected_data, oberseved_data)


def _check_capri_clt_tsv(
    capri_file: str, expected_data: list[dict[str, Union[int, str, float]]]
):
    """Helper function to check the content of the capri_clt.tsv file."""
    with open(capri_file) as f:
        lines = f.readlines()

    # There are several `#` lines in the file, these are comments and can be ignored
    lines = [line for line in lines if not line.startswith("#")]

    # Check header
    expected_header_cols = list(expected_data[0].keys())
    observed_header_cols = lines[0].strip().split("\t")

    # Check if all the columns are present
    assert len(observed_header_cols) == len(expected_header_cols), "Header mismatch"

    for col in expected_header_cols:
        assert col in observed_header_cols, f"{col} not found in the header"

    data = lines[1:]
    oberseved_data: list[dict[str, Union[int, str, float]]] = []
    for line in data:
        values = line.strip().split("\t")

        # Check if there is one value for each column
        assert len(values) == len(expected_header_cols), "Values mismatch"

        data_dic = {}
        for h, v in zip(expected_header_cols, values):
            data_dic[h] = _cast_float_str_int(v)

        oberseved_data.append(data_dic)

    assert len(oberseved_data) == len(expected_data), "Data mismatch"

    _compare_polymorphic_data(expected_data, oberseved_data)


def evaluate_caprieval_execution(
    module: CaprievalModule, model_list, ss_data, clt_data
):
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

    _check_capri_ss_tsv(
        capri_file=str(Path(module.path, "capri_ss.tsv")),
        expected_data=ss_data,
    )

    _check_capri_clt_tsv(
        capri_file=str(Path(module.path, "capri_clt.tsv")),
        expected_data=clt_data,
    )


def test_caprieval_default(
    caprieval_module, model_list, expected_ss_data, expected_clt_data
):

    caprieval_module.previous_io = MockPreviousIO(path=caprieval_module.path)
    caprieval_module.run()

    evaluate_caprieval_execution(
        caprieval_module, model_list, expected_ss_data, expected_clt_data
    )


def test_caprieval_less_io(
    caprieval_module, model_list, expected_ss_data, expected_clt_data
):
    caprieval_module.previous_io = MockPreviousIO(path=caprieval_module.path)
    caprieval_module.params["less_io"] = True

    caprieval_module.run()

    evaluate_caprieval_execution(
        caprieval_module, model_list, expected_ss_data, expected_clt_data
    )
