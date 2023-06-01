"""Test the libclust library."""
import os
from pathlib import Path

import pytest

from haddock.libs.libclust import write_structure_list
from haddock.libs.libontology import PDBFile

from . import golden_data


@pytest.fixture
def protprot_input_models():
    """Prot-prot input."""
    return [
        PDBFile(
            Path(golden_data, "protprot_complex_1.pdb"),
            path=golden_data,
            score=-42
            ),
        PDBFile(
            Path(golden_data, "protprot_complex_2.pdb"),
            path=golden_data,
            score=-17
            )
        ]


def test_write_structure_list(protprot_input_models):
    """Test write_structure_list function."""
    # fake clustering
    clustered_models = [protprot_input_models[0]]
    clustered_models[0].clt_id = 1
    cl_fname = "clustfcc.tsv"
    write_structure_list(protprot_input_models, clustered_models, cl_fname)
    observed_file_content = open(cl_fname, "r").read()
    expected_file_content = (
        f'rank\tmodel_name\tscore\tcluster_id{os.linesep}'
        f'1\tprotprot_complex_1.pdb\t-42.00\t1{os.linesep}'
        f'2\tprotprot_complex_2.pdb\t-17.00\t-{os.linesep}'
        f'{os.linesep}'
        )
    assert observed_file_content == expected_file_content
    os.unlink(cl_fname)
