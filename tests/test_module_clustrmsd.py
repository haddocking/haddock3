"""Test the clustrmsd module."""

import os
import tempfile
from pathlib import Path

import numpy as np
import pytest

from haddock.libs.libontology import ModuleIO, RMSDFile
from haddock.modules.analysis.clustrmsd import DEFAULT_CONFIG as clustrmsd_pars
from haddock.modules.analysis.clustrmsd import HaddockModule
from haddock.modules.analysis.clustrmsd.clustrmsd import (
    apply_min_population,
    cond_index,
    get_cluster_center,
    get_clusters,
    get_dendrogram,
    iterate_min_population,
    order_clusters,
    read_matrix,
    )
from haddock.modules.analysis.rmsdmatrix import DEFAULT_CONFIG as rmsd_pars
from haddock.modules.analysis.rmsdmatrix import HaddockModule as HaddockRMSD


@pytest.fixture(name="output_list")
def fixture_output_list():
    """Clustrmsd output list."""
    return [
        "rmsd.matrix",
        "rmsd_matrix.json",
        "cluster.out",
        "clustrmsd.txt",
        "clustrmsd.tsv",
        "io.json",
    ]


def remove_clustrmsd_files(output_list):
    """Remove clustrmsd files."""
    for f in output_list:
        path_f = Path(f)
        if path_f.exists():
            os.unlink(path_f)


def write_rmsd_matrix(output_name, rmsd_vec):
    """Write RMSD matrix."""
    with open(output_name, "w") as wf:
        for data in rmsd_vec:
            data_str = f"{data[0]:.0f} {data[1]:.0f} {data[2]:.3f}"
            data_str += os.linesep
            wf.write(data_str)


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


def save_rmsd_json(output_name, json_name, npairs):
    """Save the RMSD json."""
    matrix_io = ModuleIO()
    rmsd_matrix_file = RMSDFile(output_name, npairs=npairs)
    matrix_io.add(rmsd_matrix_file)
    matrix_io.save(filename=json_name)


def read_rmsd_json(json_name):
    """Read the RMSD json."""
    matrix_json = ModuleIO()
    matrix_json.load(json_name)
    return matrix_json


@pytest.fixture(name="correct_rmsd_vec")
def fixture_correct_rmsd_vec():
    """Correct RMSD vector."""
    return [[1, 2, 1.234], [1, 3, 5.678], [2, 3, 4.567]]


@pytest.fixture(name="short_rmsd_vec")
def fixture_short_rmsd_vec():
    """Short RMSD vector."""
    return [[1, 2, 1.234], [1, 3, 5.678]]


@pytest.fixture(name="malformed_rmsd_vec")
def fixture_malformed_rmsd_vec():
    """Malformed RMSD vector."""
    return [[1, 2, 1.234], [1, 3], [2, 3, 4.567]]


@pytest.fixture(name="correct_rmsd_array")
def fixture_correct_rmsd_array():
    """Correct RMSD array."""
    return np.array([1.0, 2.0, 2.5, 2.0, 3.0, 0.5])


def test_read_rmsd_matrix(correct_rmsd_vec):
    """Check correct reading of rmsd matrix."""
    rmsd_vec = correct_rmsd_vec

    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)

        output_name = Path(Path.cwd(), "fake_rmsd.matrix")
        json_name = Path(Path.cwd(), "fake_rmsd.json")

        write_rmsd_matrix(output_name, rmsd_vec)

        save_rmsd_json(output_name, json_name, 3)

        matrix_json = read_rmsd_json(json_name)

        matrix = read_matrix(matrix_json.input[0])

        assert len(matrix) == len(rmsd_vec)

        for n in range(len(rmsd_vec)):
            assert rmsd_vec[n][2] == matrix[n]


def test_read_matrix_input(correct_rmsd_vec):
    """Test wrong input to read_matrix."""
    rmsd_vec = correct_rmsd_vec

    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)
        output_name = Path(Path.cwd(), "fake_rmsd.matrix")

        write_rmsd_matrix(output_name, rmsd_vec)

        # testing input : has to be RMSDFile
        with pytest.raises(TypeError):
            read_matrix(output_name)

        with pytest.raises(TypeError):
            read_matrix(Path(output_name))


def test_read_matrix_binomial(short_rmsd_vec):
    """Test rmsd matrix with length != binomial coefficient."""
    rmsd_vec = short_rmsd_vec

    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)
        output_name = Path(Path.cwd(), "fake_rmsd.matrix")
        json_name = Path(Path.cwd(), "fake_rmsd.json")

        write_rmsd_matrix(output_name, rmsd_vec)

        save_rmsd_json(output_name, json_name, 3)

        matrix_json = read_rmsd_json(json_name)

        with pytest.raises(ValueError):
            read_matrix(matrix_json.input[0])


def test_read_matrix_npairs(correct_rmsd_vec):
    """Test rmsd file with npairs != binomial coefficient."""
    rmsd_vec = correct_rmsd_vec

    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)
        output_name = Path(Path.cwd(), "fake_rmsd.matrix")
        json_name = Path(Path.cwd(), "fake_rmsd.json")

        write_rmsd_matrix(output_name, rmsd_vec)

        save_rmsd_json(output_name, json_name, 2)

        matrix_json = read_rmsd_json(json_name)

        with pytest.raises(ValueError):
            read_matrix(matrix_json.input[0])


def test_read_matrix_malformed(malformed_rmsd_vec):
    """Test read malformed rmsd file."""
    rmsd_vec = malformed_rmsd_vec

    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)
        output_name = Path(Path.cwd(), "fake_rmsd.matrix")
        json_name = Path(Path.cwd(), "fake_rmsd.json")

        write_malformed_rmsd_matrix(output_name, rmsd_vec)

        save_rmsd_json(output_name, json_name, 3)

        matrix_json = read_rmsd_json(json_name)

        with pytest.raises(ValueError):
            read_matrix(matrix_json.input[0])


def test_correct_clusters(correct_rmsd_array):
    """Test correct average-linkage hierarchical clustering."""
    rmsd_matrix = correct_rmsd_array

    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)
        observed_dendrogram = get_dendrogram(rmsd_matrix, linkage_type="average")

        expected_dendrogram = np.array(
            [[2.0, 3.0, 0.5, 2.0], [0.0, 1.0, 1.0, 2.0], [4.0, 5.0, 2.375, 4.0]]
        )

        assert (observed_dendrogram == expected_dendrogram).all()

        observed_clusters = get_clusters(
            observed_dendrogram, tolerance=2, criterion="maxclust"
        )

        expected_clusters = np.array([2, 2, 1, 1])

        assert (observed_clusters == expected_clusters).all()


# TODO: add tests for the other categories of clustering


def test_correct_output(protdna_input_list, output_list):
    """Test correct clustrmsd output."""
    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)
        rmsd_module = HaddockRMSD(order=2, path=Path(""), initial_params=rmsd_pars)
        rmsd_module.previous_io.output = protdna_input_list
        rmsd_module._run()

        clustrmsd_module = HaddockModule(
            order=3, path=Path(""), initial_params=clustrmsd_pars
        )
        clustrmsd_module._run()

        ls = os.listdir()

        expected_out_filename = "cluster.out"
        expected_txt_filename = "clustrmsd.txt"

        assert expected_out_filename in ls

        assert expected_txt_filename in ls

        expected_out_content = f"Cluster 1 -> 1 2{os.linesep}"

        observed_out_content = open(expected_out_filename).read()

        assert observed_out_content == expected_out_content

        remove_clustrmsd_files(output_list)


def test_get_cluster_center(correct_rmsd_array):
    """Test get_cluster_center function."""
    obs_clt_center = get_cluster_center(
        npw=[0, 1, 2, 3], n_obs=4, rmsd_matrix=correct_rmsd_array
    )
    exp_clt_center = 2
    assert obs_clt_center == exp_clt_center
    # [1,2,3] cluster. 2 is still the center
    obs_clt_center = get_cluster_center(
        npw=[1, 2, 3], n_obs=4, rmsd_matrix=correct_rmsd_array
    )
    exp_clt_center = 2
    assert obs_clt_center == exp_clt_center


def test_cond_index():
    """Test cond_index function."""
    n_obs = 10
    idxs = [0, 0, 2, 8]
    jdxs = [1, 2, 3, 9]
    len_idxs = len(idxs)
    obs_c_idxs = [cond_index(idxs[n], jdxs[n], n_obs) for n in range(len_idxs)]
    exp_c_idxs = [0, 1, 17, 44]
    assert obs_c_idxs == exp_c_idxs


def test_apply_min_population():
    """Test apply_min_population function."""
    # defining cluster_arr
    cluster_arr = np.array([1, 1, 4, 1, 1, 2, 1, 1, 3, 1])
    # using a min_population of 2
    obs_cluster_arr = apply_min_population(cluster_arr, 2)
    exp_cluster_arr = np.array([1, 1, -1, 1, 1, -1, 1, 1, -1, 1])
    assert (obs_cluster_arr == exp_cluster_arr).all()
    # using a min_population of 1
    obs_cluster_arr = apply_min_population(cluster_arr, min_population=1)
    exp_cluster_arr = np.array([1, 1, 4, 1, 1, 2, 1, 1, 3, 1])
    assert (obs_cluster_arr == exp_cluster_arr).all()


def test_iterate_min_population():
    """Test iterate_min_population function."""
    cluster_arr = np.array([1, 1, 2, 3, 4])
    obs_cluster_arr, obs_min_population = iterate_min_population(
        cluster_arr,
        min_population=4,
    )
    exp_cluster_arr = np.array([1, 1, -1, -1, -1])
    assert obs_min_population == 2
    assert (obs_cluster_arr == exp_cluster_arr).all()


def test_order_clusters():
    """Test order_clusters function."""
    cluster_arr = np.array([1, 1, 2, 3, 4])
    obs_clusters, obs_cluster_arr = order_clusters(cluster_arr)
    exp_clusters = [1, 2, 3, 4]
    exp_cluster_arr = np.array([1, 1, 2, 3, 4])
    assert obs_clusters == exp_clusters
    assert (obs_cluster_arr == exp_cluster_arr).all()
    # now with a less trivial cluster_arr
    cluster_arr = np.array([3, 3, 2, 4, 3, 1, 3, 3, 1, 3, 4])
    obs_clusters, obs_cluster_arr = order_clusters(cluster_arr)
    exp_clusters = [1, 2, 3, 4]
    exp_cluster_arr = np.array([1, 1, 4, 3, 1, 2, 1, 1, 2, 1, 3])
    assert obs_clusters == exp_clusters
    assert (obs_cluster_arr == exp_cluster_arr).all()
    # yet another cluster_arr with unclustered structures
    cluster_arr = np.array([3, 3, 2, 4, -1, -1, 3, 3, -1, -1, 4, 1])
    obs_clusters, obs_cluster_arr = order_clusters(cluster_arr)
    exp_clusters = [1, 2, 3, 4]
    exp_cluster_arr = np.array([1, 1, 4, 2, -1, -1, 1, 1, -1, -1, 2, 3])
    assert obs_clusters == exp_clusters
    assert (obs_cluster_arr == exp_cluster_arr).all()
