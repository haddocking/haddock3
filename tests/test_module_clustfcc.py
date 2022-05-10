"""Test the clustfcc module."""
import os
from pathlib import Path

import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.clustfcc import DEFAULT_CONFIG as clustfcc_pars
from haddock.modules.analysis.clustfcc import HaddockModule

from . import golden_data


@pytest.fixture
def protprot_input_list():
    """Prot-prot input."""
    return [
        PDBFile(Path(golden_data, "protprot_complex_1.pdb"), path=golden_data),
        PDBFile(Path(golden_data, "protprot_complex_2.pdb"), path=golden_data)
        ]


@pytest.fixture
def output_list():
    """Clustfcc output list."""
    return [
        "fcc.matrix",
        "cluster.out",
        "protprot_complex_1.con",
        "protprot_complex_2.con",
        "clustfcc.txt",
        "io.json"
        ]


def test_clustfcc_output_existence(protprot_input_list, output_list):
    """Test clustfcc output."""
    fcc_module = HaddockModule(
        order=1,
        path=Path("1_emscoring"),
        initial_params=clustfcc_pars
        )
    fcc_module.previous_io.output = protprot_input_list

    fcc_module._run()

    ls = os.listdir()

    for el in output_list:
        assert el in ls


def test_matrix_output():
    """Check fcc.matrix file."""
    fcc_file = Path("fcc.matrix")

    observed_output = open(fcc_file).read()

    expected_output = "1 2 0.05 0.062" + os.linesep

    assert observed_output == expected_output


def test_contacts():
    """Check .con files."""
    expected_output_length = [100, 119]

    observed_output_lengths = []

    con1_len = len(open("protprot_complex_1.con").readlines())
    observed_output_lengths.append(con1_len)

    con2_len = len(open("protprot_complex_2.con").readlines())
    observed_output_lengths.append(con2_len)

    assert observed_output_lengths == expected_output_length


def remove_clustfcc_files(output_list):
    """Remove clustfcc files."""
    for f in output_list:
        path_f = Path(f)
        if path_f.exists():
            os.unlink(path_f)


def test_cluster_out(output_list):
    """Check cluster.out file."""
    expected_output = "Cluster 1 -> 2 " + os.linesep
    expected_output += "Cluster 2 -> 1 " + os.linesep

    observed_output = open("cluster.out").read()

    assert expected_output == observed_output

    remove_clustfcc_files(output_list)
