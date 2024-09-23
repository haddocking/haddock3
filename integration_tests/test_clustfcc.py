import os
import shutil
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.clustfcc import DEFAULT_CONFIG as clustfcc_pars
from haddock.modules.analysis.clustfcc import HaddockModule as ClustFCCModule

from integration_tests import GOLDEN_DATA


class MockPreviousIO:
    """MockPreviousIO injects models into the test class"""

    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = False):
        shutil.copy(
            Path(GOLDEN_DATA, "protprot_complex_1.pdb"),
            Path(self.path, "protprot_complex_1.pdb"),
        )

        shutil.copy(
            Path(GOLDEN_DATA, "protprot_complex_2.pdb"),
            Path(self.path, "protprot_complex_2.pdb"),
        )

        model_list = [
            PDBFile(
                file_name="protprot_complex_1.pdb",
                path=self.path,
            ),
            PDBFile(
                file_name="protprot_complex_2.pdb",
                path=self.path,
            ),
        ]
        return model_list

    def output(self) -> None:
        """Placeholder for the output method"""
        return None


@pytest.fixture(name="output_list")
def fixture_output_list():
    """Clustfcc output list."""
    return [
        "fcc.matrix",
        "cluster.out",
        "protprot_complex_1.con",
        "protprot_complex_2.con",
        "clustfcc.txt",
        "io.json",
        "clustfcc.tsv",
        "fcc_matrix.html",
    ]


@pytest.fixture(name="fcc_module")
def fixture_fcc_module():
    """Clustfcc module."""
    with tempfile.TemporaryDirectory() as tempdir:
        yield ClustFCCModule(order=1, path=Path(tempdir), initial_params=clustfcc_pars)


def test_clustfcc(fcc_module, output_list):
    """Test clustfcc output."""
    fcc_module.previous_io = MockPreviousIO(path=fcc_module.path)
    fcc_module.params["plot_matrix"] = True

    fcc_module.run()

    for _f in output_list:
        expected_file = Path(fcc_module.path, _f)

        # Check that the file exists
        assert expected_file.exists()

        # check if the file is not empty
        assert expected_file.stat().st_size > 0

    # Test the fcc matrix contents
    with open(Path(fcc_module.path, "fcc.matrix"), encoding="utf-8", mode="r") as f:
        observed_fcc_matrix = f.read()

    expected_fcc_output = "1 2 0.05 0.062" + os.linesep
    assert observed_fcc_matrix == expected_fcc_output

    # Check .con files
    expected_output_length = [100, 119]

    observed_contact_files = [
        Path(fcc_module.path, "protprot_complex_1.con"),
        Path(fcc_module.path, "protprot_complex_2.con"),
    ]

    for exp, obs_f in zip(expected_output_length, observed_contact_files):
        with open(obs_f, encoding="utf-8", mode="r") as f:
            obs = len(f.read().splitlines())
            assert obs == exp

    # Check cluster.out file
    expected_cluster_output = [
        "Cluster 1 -> 2 " + os.linesep,
        "Cluster 2 -> 1 " + os.linesep,
    ]
    with open(Path(fcc_module.path, "cluster.out"), encoding="utf-8", mode="r") as f:
        for expected, line in zip(expected_cluster_output, f.readlines()):
            assert line == expected
