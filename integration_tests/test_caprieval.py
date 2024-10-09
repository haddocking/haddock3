import math
import shutil
import tempfile
import random
from pathlib import Path
from typing import Union

import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.caprieval import (
    DEFAULT_CONFIG as DEFAULT_CAPRIEVAL_CONFIG,
    capri,
)
from haddock.modules.analysis.caprieval import HaddockModule as CaprievalModule
from tests import golden_data as UNITTESTS_GOLDEN_DATA
from . import GOLDEN_DATA


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
            score=0.0,
        ),
        PDBFile(
            file_name="protprot_complex_2.pdb",
            path=".",
            unw_energies={"energy_term": 0.0},
            score=100.0,
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
            "score": 50.0,
            "score_std": 50.0,
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
            "rmsd": 4.311,
            "rmsd_std": 4.311,
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
            "score": 0.0,
            "irmsd": 0.000,
            "fnat": 1.000,
            "lrmsd": 0.000,
            "ilrmsd": 0.000,
            "dockq": 1.000,
            "rmsd": 0.000,
            "cluster_id": "-",
            "cluster_ranking": "-",
            "model-cluster_ranking": "-",
            "energy_term": 0.000,
        },
        {
            "model": "",
            "md5": "-",
            "caprieval_rank": 2,
            "score": 100.0,
            "irmsd": 8.327,
            "fnat": 0.050,
            "lrmsd": 20.937,
            "ilrmsd": 18.248,
            "dockq": 0.074,
            "rmsd": 8.623,
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
            Path(UNITTESTS_GOLDEN_DATA, "protprot_complex_1.pdb"),
            Path(".", "protprot_complex_1.pdb"),
        )
        shutil.copy(
            Path(UNITTESTS_GOLDEN_DATA, "protprot_complex_2.pdb"),
            Path(".", "protprot_complex_2.pdb"),
        )
        model_list = [
            PDBFile(
                file_name="protprot_complex_1.pdb",
                path=".",
                unw_energies={"energy_term": 0.0},
                score=0.0,
            ),
            PDBFile(
                file_name="protprot_complex_2.pdb",
                path=".",
                unw_energies={"energy_term": 0.0},
                score=100.0,
            ),
        ]

        return model_list

    def output(self):
        return None


class MockPreviousIO_with_models_to_be_clustered:
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = True):
        model_list = []
        for m in [f"rigidbody_{i}.pdb" for i in range(1, 11)]:
            src = Path(GOLDEN_DATA, "models_for_clustering", m)
            dst = Path(self.path, m)
            shutil.copy(src, dst)
            p = PDBFile(file_name=m, path=".", score=random.uniform(-100, 100))
            model_list.append(p)
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
    for expected, observed in zip(expected_data, oberseved_data):
        # Loop over keys
        for key in expected.keys():
            # Point corresponding values
            exp_val = expected[key]
            obs_val = observed[key]

            # Check the type match
            assert type(exp_val) == type(obs_val), f"Type mismatch for {key}"

            # Check float
            if isinstance(exp_val, float):
                if math.isnan(exp_val):
                    assert isinstance(obs_val, float) and math.isnan(
                        obs_val
                    ), f"Value mismatch for {key}"
                else:
                    assert (
                        isinstance(obs_val, (int, float))
                        and pytest.approx(exp_val, rel=1e-3) == obs_val
                    ), f"Value mismatch for {key}"

            # Check int
            elif isinstance(exp_val, int):
                assert exp_val == obs_val, f"Value mismatch for {key}"

            # Check str
            elif isinstance(exp_val, str):
                assert exp_val == obs_val, f"Value mismatch for {key}"

            # Value is not float, int or str
            else:
                raise ValueError(f"Unknown type for {key}")


def _check_capri_ss_tsv(
    capri_file: str, expected_data: list[dict[str, Union[int, str, float]]]
):
    """Helper function to check the content of the capri_ss.tsv file."""
    with open(capri_file) as f:
        lines = [_.strip().split("\t") for _ in f.readlines() if not _.startswith("#")]

    # Check the header
    expected_header_cols = list(expected_data[0].keys())
    observed_header_cols = lines[0]

    # Check if all of them have the same lenght
    assert len(observed_header_cols) == len(expected_header_cols), "Header mismatch"
    # Make sure column names are the same
    for col_name in expected_header_cols:
        assert col_name in observed_header_cols, f"{col_name} not found in the header"

    oberseved_data: list[dict[str, Union[int, str, float]]] = []
    data = lines[1:]
    for line in data:

        # Check there is one value for each column
        assert len(line) == len(expected_header_cols), "Values mismatch"

        data_dict = {}
        for h, v in zip(expected_header_cols, line):
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
        # There are several `#` lines in the file, these are comments and can be ignored
        lines = [_.strip().split("\t") for _ in f.readlines() if not _.startswith("#")]

    # Check header
    expected_header_cols = list(expected_data[0].keys())
    observed_header_cols = lines[0]

    # Check if all the columns are present
    assert len(observed_header_cols) == len(expected_header_cols), "Header mismatch"

    for col in expected_header_cols:
        assert col in observed_header_cols, f"{col} not found in the header"

    data = lines[1:]
    oberseved_data: list[dict[str, Union[int, str, float]]] = []
    for line in data:

        # Check if there is one value for each column
        assert len(line) == len(expected_header_cols), "Values mismatch"

        data_dic = {}
        for h, v in zip(expected_header_cols, line):
            data_dic[h] = _cast_float_str_int(v)

        oberseved_data.append(data_dic)

    assert len(oberseved_data) == len(expected_data), "Data mismatch"

    _compare_polymorphic_data(expected_data, oberseved_data)


def _check_means_match(
    capri_ss_f: Path,
    capri_clt_f: Path,
    target_metric: str,
    top_n: int,
):
    """Helper function to check if the means of `capri_ss` and `capri_clt` match"""

    # Read the `capri_ss.tsv` file
    assert capri_ss_f.exists()

    assert capri_clt_f.exists()

    # find the column that contains the target metric
    with open(capri_ss_f, "r") as fh:
        capri_ss_l = fh.readlines()

    capri_ss_header = capri_ss_l[0].strip().split("\t")
    capri_ss_data = capri_ss_l[1:]
    try:
        metric_idx = capri_ss_header.index(target_metric)
    except ValueError:
        # Metric not found
        return

    assert metric_idx is not None
    values: list[float] = []
    for entry in capri_ss_data:
        data = entry.strip().split(sep="\t")
        v = data[metric_idx]
        ranking = int(data[2])
        if ranking <= top_n:
            values.append(float(v) if v != "-" else float("nan"))

    # get the topX and calculate mean
    mean_ss_v = sum(values[:top_n]) / float(top_n)

    with open(capri_clt_f, "r") as fh:
        capri_clt_l = fh.readlines()

    # remove lines with comments
    capri_clt_l = [line for line in capri_clt_l if "#" not in line]

    # find the metric index
    capri_clt_header = capri_clt_l[0].strip().split("\t")
    capri_clt_data = capri_clt_l[1].strip().split("\t")
    for metric_idx, metric in enumerate(capri_clt_header):
        if metric == target_metric:
            break

    assert metric_idx is not None
    mean_clt_v = float(capri_clt_data[metric_idx])

    assert mean_clt_v == pytest.approx(mean_ss_v, 0.01), f"{target_metric} do not match"


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


def test_ss_clt_relation(caprieval_module):
    """Check if the values in the ss.tsv match the ones in clt.tsv"""

    caprieval_module.previous_io = MockPreviousIO_with_models_to_be_clustered(
        path=caprieval_module.path
    )

    caprieval_module.run()

    metrics_to_be_evaluated = [
        "score",
        "irmsd",
        "fnat",
        "lrmsd",
        "dockq",
        "ilrmsd",
        "rmsd",
    ]

    for metric in metrics_to_be_evaluated:
        _check_means_match(
            capri_ss_f=Path(caprieval_module.path, "capri_ss.tsv"),
            capri_clt_f=Path(caprieval_module.path, "capri_clt.tsv"),
            target_metric=metric,
            top_n=caprieval_module.params["clt_threshold"],
        )
