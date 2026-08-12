"""Test the libclust library."""

import os
import random
import tempfile
import numpy as np
from pathlib import Path
from typing import cast

import pytest

from haddock.libs.libclust import (
    MAX_NB_ENTRY_HTML_MATRIX,
    add_cluster_info,
    plot_cluster_matrix,
    write_structure_list,
)
from haddock.libs.libontology import PDBFile

from . import golden_data


@pytest.fixture
def protprot_input_models():
    """Prot-prot input."""
    return [
        PDBFile(
            Path(golden_data, "protprot_complex_1.pdb"), path=golden_data, score=-42
        ),
        PDBFile(
            Path(golden_data, "protprot_complex_2.pdb"), path=golden_data, score=-17
        ),
    ]


@pytest.fixture
def big_distance_matrix_data():
    """Big distance matrix data."""
    bigmatrix = []
    for i in range(1, MAX_NB_ENTRY_HTML_MATRIX + 1):
        for j in range(i + 1, MAX_NB_ENTRY_HTML_MATRIX + 2):
            bigmatrix.append(f"{i}\t{j}\t{random.random():.2f}")
    return bigmatrix


@pytest.fixture
def small_distance_matrix_data():
    """Distance matrix data."""
    small = 10
    matrix = []
    for i in range(1, small):
        for j in range(i + 1, small + 1):
            matrix.append(f"{i}\t{j}\t{random.random():.2f}")
    return matrix


def test_write_structure_list(protprot_input_models):
    """Test write_structure_list function."""
    # fake clustering
    clustered_models = [protprot_input_models[0]]
    clustered_models[0].clt_id = 1
    cl_fname = "clustfcc.tsv"
    write_structure_list(protprot_input_models, clustered_models, cl_fname)
    observed_file_content = open(cl_fname, "r").read()
    expected_file_content = (
        f"rank\tmodel_name\tscore\tcluster_id{os.linesep}"
        f"1\tprotprot_complex_1.pdb\t-42.00\t1{os.linesep}"
        f"2\tprotprot_complex_2.pdb\t-17.00\t-{os.linesep}"
        f"{os.linesep}"
    )
    assert observed_file_content == expected_file_content
    os.unlink(cl_fname)


def test_plot_cluster_matrix_big(big_distance_matrix_data):
    """Test test_plot_cluster_matrix_big function with big matrix."""
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        # Write matrix
        matrix_path = f"{tmpdir}bigmatrix.mat"
        with open(matrix_path, "w") as f:
            f.write("\n".join(big_distance_matrix_data))
        # Run function
        figure_path = plot_cluster_matrix(
            matrix_path,
            list(range(MAX_NB_ENTRY_HTML_MATRIX + 1)),
            [str(i + 1) for i in range(MAX_NB_ENTRY_HTML_MATRIX + 1)],
            dttype="random",
            diag_fill=0.5,
            color_scale="Blues",
            reverse=False,
            output_fname="clust_matrix_test_big",
        )
        # Check output
        assert figure_path is None
        assert not Path("clust_matrix_test_big.html").exists()
        Path(matrix_path).unlink(missing_ok=False)


def test_plot_cluster_matrix(small_distance_matrix_data):
    """Test test_plot_cluster_matrix_big function with small matrix."""
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        # Write matrix
        matrix_path = f"{tmpdir}smallmatrix.mat"
        with open(matrix_path, "w") as f:
            f.write("\n".join(small_distance_matrix_data))
        # Run function
        figure_path = plot_cluster_matrix(
            matrix_path,
            list(range(10)),
            [str(i + 1) for i in range(10)],
            dttype="random",
            diag_fill=0.5,
            color_scale="Blues",
            reverse=True,
            output_fname="clust_matrix_test_small",
        )
        # Check output
        assert os.path.exists(figure_path)
        assert Path(figure_path).stat().st_size != 0
        assert Path(figure_path).suffix == ".html"
        Path(figure_path).unlink(missing_ok=False)
        Path(matrix_path).unlink(missing_ok=False)


def test_add_cluster_info(protprot_input_models):

    # first 5 -> pdbs_cluster_2 (worse cluster), last 5 -> pdbs_cluster_1 (better cluster)
    scores = [+10.0, +20.0, +30.0, +40.0, +50.0, -410.0, -420.0, -430.0, -440.0, -450.0]
    input_pdbs = [
        PDBFile(file_name=f"{i}.pdb", score=score) for i, score in enumerate(scores)
    ]

    # split so we get two sets of models
    mid = len(input_pdbs) // 2
    pdbs_cluster_1 = input_pdbs[mid:]
    avg_1 = float(np.average(scores[mid:]))
    pdbs_cluster_2 = input_pdbs[:mid]
    avg_2 = np.average(scores[:mid])

    # create cluster information
    # [(cluster_id, average-score)]
    input_l = [(1, avg_1), (2, avg_2)]

    # Make sure the PDBFiles are not created with any cluster information
    for p in input_pdbs:
        assert p.clt_id is None
        assert p.clt_rank is None
        assert p.clt_model_rank is None

    observed = add_cluster_info(input_l, clt_dic={1: pdbs_cluster_1, 2: pdbs_cluster_2})

    # assert that the cluster information was added to the models
    for p in observed:
        assert p.clt_id is not None
        assert p.clt_rank is not None
        assert p.clt_model_rank is not None

    # assert the order has changed
    assert observed[0] is not input_pdbs[0]

    # check they were properly sorted
    assert cast(float, observed[0].score) < cast(float, observed[-1].score)
