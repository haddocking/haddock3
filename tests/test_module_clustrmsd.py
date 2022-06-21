"""Test the clustrmsd module."""
import os
from pathlib import Path

import numpy as np
import pytest

from haddock.libs.libontology import ModuleIO, PDBFile, RMSDFile
from haddock.modules.analysis.clustrmsd import DEFAULT_CONFIG as clustrmsd_pars
from haddock.modules.analysis.clustrmsd import HaddockModule
from haddock.modules.analysis.clustrmsd.clustrmsd import (
    get_clusters,
    get_dendrogram,
    read_matrix,
    )
from haddock.modules.analysis.rmsdmatrix import DEFAULT_CONFIG as rmsd_pars
from haddock.modules.analysis.rmsdmatrix import HaddockModule as HaddockRMSD

from . import golden_data


def write_rmsd_matrix(output_name, rmsd_vec):
    """Write RMSD matrix."""
    with open(output_name, "w") as wf:
        for data in rmsd_vec:
            data_str = f"{data[0]:.0f} {data[1]:.0f} {data[2]:.3f}"
            data_str += os.linesep
            wf.write(data_str)
    return


def write_malformed_rmsd_matrix(output_name, rmsd_vec):
    """Write RMSD matrix."""
    with open(output_name, "w") as wf:
        for data in rmsd_vec:
            data_str = ""
            for el in data:
                data_str += " "
                data_str += str(el)
            data_str += os.linesep
            wf.write(data_str)
    return


def save_rmsd_json(output_name, json_name, npairs):
    """Save the RMSD json."""
    matrix_io = ModuleIO()
    rmsd_matrix_file = RMSDFile(
        output_name,
        npairs=npairs
        )
    matrix_io.add(rmsd_matrix_file)
    matrix_io.save(filename=json_name)
    return


def read_rmsd_json(json_name):
    """Read the RMSD json."""
    matrix_json = ModuleIO()
    matrix_json.load(json_name)
    return matrix_json


@pytest.fixture
def correct_rmsd_vec():
    """Correct RMSD vector."""
    return [[1, 2, 1.234],
            [1, 3, 5.678],
            [2, 3, 4.567]]


@pytest.fixture
def short_rmsd_vec():
    """Short RMSD vector."""
    return [[1, 2, 1.234],
            [1, 3, 5.678]]


@pytest.fixture
def malformed_rmsd_vec():
    """Malformed RMSD vector."""
    return [[1, 2, 1.234],
            [1, 3],
            [2, 3, 4.567]]


@pytest.fixture
def correct_rmsd_array():
    """Correct RMSD array."""
    return np.array([1.0, 2.0, 2.5, 2.0, 3.0, 0.5])


@pytest.fixture
def input_protdna_models():
    """Prot-DNA models using for emscoring output."""
    return [
        PDBFile(
            Path(golden_data, "protdna_complex_1.pdb"),
            path=golden_data,
            score=42.0
            ),
        PDBFile(
            Path(golden_data, "protdna_complex_2.pdb"),
            path=golden_data,
            score=28.0
            )]


def test_read_rmsd_matrix(correct_rmsd_vec):
    """Check correct reading of rmsd matrix."""
    rmsd_vec = correct_rmsd_vec
    
    output_name = "fake_rmsd.matrix"
    json_name = "fake_rmsd.json"
    
    write_rmsd_matrix(output_name, rmsd_vec)

    save_rmsd_json(output_name, json_name, 3)

    matrix_json = read_rmsd_json(json_name)
    
    matrix = read_matrix(matrix_json.input[0])
    
    assert len(matrix) == len(rmsd_vec)

    for n in range(len(rmsd_vec)):
        assert rmsd_vec[n][2] == matrix[n]
    
    os.unlink(json_name)
    os.unlink(output_name)


def test_read_matrix_input(correct_rmsd_vec):
    """Test wrong input to read_matrix."""
    rmsd_vec = correct_rmsd_vec
    
    output_name = "fake_rmsd.matrix"

    write_rmsd_matrix(output_name, rmsd_vec)

    # testing input : has to be RMSDFile
    with pytest.raises(Exception):
        read_matrix(output_name)
    
    with pytest.raises(Exception):
        read_matrix(Path(output_name))

    os.unlink(output_name)


def test_read_matrix_binomial(short_rmsd_vec):
    """Test rmsd matrix with length != binomial coefficient."""
    rmsd_vec = short_rmsd_vec
    
    output_name = "fake_rmsd.matrix"
    json_name = "fake_rmsd.json"
    
    write_rmsd_matrix(output_name, rmsd_vec)

    save_rmsd_json(output_name, json_name, 3)

    matrix_json = read_rmsd_json(json_name)

    with pytest.raises(Exception):
        read_matrix(matrix_json)
    
    os.unlink(json_name)
    os.unlink(output_name)


def test_read_matrix_npairs(correct_rmsd_vec):
    """Test rmsd file with npairs != binomial coefficient."""
    rmsd_vec = correct_rmsd_vec
    
    output_name = "fake_rmsd.matrix"
    json_name = "fake_rmsd.json"
    
    write_rmsd_matrix(output_name, rmsd_vec)

    save_rmsd_json(output_name, json_name, 2)

    matrix_json = read_rmsd_json(json_name)

    with pytest.raises(Exception):
        read_matrix(matrix_json)
    
    os.unlink(json_name)
    os.unlink(output_name)


def test_read_matrix_malformed(malformed_rmsd_vec):
    """Test read malformed rmsd file."""
    rmsd_vec = malformed_rmsd_vec
    
    output_name = "fake_rmsd.matrix"
    json_name = "fake_rmsd.json"
    
    write_malformed_rmsd_matrix(output_name, rmsd_vec)

    save_rmsd_json(output_name, json_name, 3)

    matrix_json = read_rmsd_json(json_name)

    with pytest.raises(Exception):
        read_matrix(matrix_json)
    
    os.unlink(json_name)
    os.unlink(output_name)


def test_correct_clusters(correct_rmsd_array):
    """Test correct average-linkage hierarchical clustering."""
    rmsd_matrix = correct_rmsd_array

    observed_dendrogram = get_dendrogram(
        rmsd_matrix,
        linkage_type="average"
        )

    expected_dendrogram = np.array(
        [[2., 3., 0.5, 2.],
         [0., 1., 1.0, 2.],
         [4., 5., 2.375, 4.]]
        )

    assert (observed_dendrogram == expected_dendrogram).all()

    observed_clusters = get_clusters(
        observed_dendrogram,
        tolerance=2,
        criterion="maxclust"
        )
    
    expected_clusters = np.array([2, 2, 1, 1])

    assert (observed_clusters == expected_clusters).all()

# TODO: add tests for the other categories of clustering


def test_correct_output(input_protdna_models):
    """Test correct clustrmsd output."""
    rmsd_module = HaddockRMSD(
        order=2,
        path=Path(""),
        initial_params=rmsd_pars
        )
    rmsd_module.previous_io.output = input_protdna_models
    rmsd_module._run()

    clustrmsd_module = HaddockModule(
        order=3,
        path=Path(""),
        initial_params=clustrmsd_pars
        )
    clustrmsd_module._run()

    ls = os.listdir()

    expected_out_filename = "cluster.out"
    expected_txt_filename = "clustrmsd.txt"

    assert expected_out_filename in ls

    assert expected_txt_filename in ls

    expected_out_content = "Cluster 1 -> 1 2"
    expected_out_content += os.linesep

    observed_out_content = open(expected_out_filename).read()

    assert observed_out_content == expected_out_content

    # TODO: check the content of clustrmsd.txt

    os.unlink("io.json")
    os.unlink("rmsd.matrix")
    os.unlink("rmsd_matrix.json")
    os.unlink(expected_out_filename)
    os.unlink(expected_txt_filename)
