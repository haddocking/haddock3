"""Test finding the best structures in libplots."""

import os

from pathlib import Path

import pandas as pd
import pytest
import tempfile

from haddock.libs.libplots import (
    box_plot_handler,   
    create_other_cluster,
    find_best_struct,
    make_alascan_plot,
    offline_js_manager,
    scatter_plot_handler,
    read_capri_table,
    )
from . import data_folder, golden_data


@pytest.fixture(name="example_capri_ss")
def fixture_example_capri_ss():
    """Provide example capri_ss.tsv filepath."""
    return Path(golden_data, "capri_ss_example.tsv")


@pytest.fixture(name="loaded_capri_ss")
def fixture_loaded_capri_ss(example_capri_ss):
    """Provide loaded example capri_ss.tsv."""
    return read_capri_table(example_capri_ss)


@pytest.fixture(name="example_capri_clt")
def fixture_example_capri_clt():
    """Provide example capri_clt.tsv filepath."""
    return Path(golden_data, "capri_ss_example.tsv")

@pytest.fixture(name="cluster_ranking")
def fixture_cluster_ranking(example_capri_clt):
    """Provide top 10 ranking from capri_clt.tsv."""
    from haddock.clis.cli_analyse import get_cluster_ranking
    _cluster_ranking = get_cluster_ranking(
        example_capri_clt,
        10,
        )
    return _cluster_ranking


def test_find_best_struct_best1(loaded_capri_ss):
    """Finds one best structure."""
    result = find_best_struct(loaded_capri_ss, 1)
    cluster1 = result[result["cluster_id"] == 1]

    expected = pd.DataFrame(
        {
            "cluster_id": {0: 1},
            "best1": {
                0: "../1_rigidbody/rigidbody_383.pdb",
            },
        }
    )
    pd.testing.assert_frame_equal(cluster1, expected)


def test_find_best_struct_best5_withmssing(loaded_capri_ss):
    """Finds 20 best structures if possible."""
    result = find_best_struct(loaded_capri_ss, 5)
    # cluster 37 has only 4 structures
    cluster37 = result[result["cluster_id"] == 37].reset_index(drop=True)
    expected = pd.DataFrame(
        {
            "cluster_id": {0: 37},
            "best1": {0: "../1_rigidbody/rigidbody_130.pdb"},
            "best2": {0: "../1_rigidbody/rigidbody_788.pdb"},
            "best3": {0: "../1_rigidbody/rigidbody_485.pdb"},
            "best4": {0: "../1_rigidbody/rigidbody_496.pdb"},
            "best5": {0: ""},
        }
    )
    pd.testing.assert_frame_equal(cluster37, expected)


def test_find_best_struct_bestdefault(loaded_capri_ss):
    """Finds 4 best structures if possible."""
    result = find_best_struct(loaded_capri_ss)
    cluster1 = result[result["cluster_id"] == 1]

    expected = pd.DataFrame(
        {
            "cluster_id": {0: 1},
            "best1": {0: "../1_rigidbody/rigidbody_383.pdb"},
            "best2": {0: "../1_rigidbody/rigidbody_265.pdb"},
            "best3": {0: "../1_rigidbody/rigidbody_231.pdb"},
            "best4": {0: "../1_rigidbody/rigidbody_218.pdb"},
        }
    )
    pd.testing.assert_frame_equal(cluster1, expected)


def test_find_best_struct_empty_column():
    ss_df = pd.DataFrame(
        {
            "cluster_id": [1, 1, 1],
            "model": [
                "../1_rigidbody/rigidbody_383.pdb",
                "../1_rigidbody/rigidbody_265.pdb",
                "../1_rigidbody/rigidbody_231.pdb",
            ],
            "model-cluster_ranking": [1, 2, 3],
            # Dont care about the rest of the columns, they are ignored by the function
        }
    )
    result = find_best_struct(ss_df, 1)
    expected = pd.DataFrame(
        {
            "cluster_id": {0: 1},
            "best1": {
                0: "../1_rigidbody/rigidbody_383.pdb",
            },
            # best2, best3 have been removed as they are empty
        }
    )
    pd.testing.assert_frame_equal(result, expected)


def test_create_other_cluster_nochange():
    clusters_df = pd.DataFrame(
        {
            "cluster_id": [1],
            "cluster_rank": [1],
            "n": [1],
            "caprieval_rank": [1],
            # Dont care about the rest of the columns, they are ignored by the function
        }
    )
    structs_df = pd.DataFrame(
        {
            "caprieval_rank": [1],
            "cluster_ranking": [1],
            "cluster_id": [1],
            "model-cluster_ranking": [1],
            # Dont care about the rest of the columns, they are ignored by the function
        }
    )
    new_clusters_df, new_structs_df = create_other_cluster(clusters_df, structs_df, 2)
    pd.testing.assert_frame_equal(new_structs_df, structs_df)
    pd.testing.assert_frame_equal(new_clusters_df, clusters_df)


def test_create_other_cluster_maxsamelen():
    clusters_df = pd.DataFrame(
        {
            "cluster_id": [1, 2],
            "cluster_rank": [1, 2],
            "n": [1, 1],
            "caprieval_rank": [1, 2],
            # Dont care about the rest of the columns, they are ignored by the function
        }
    )
    structs_df = pd.DataFrame(
        {
            "caprieval_rank": [1, 2],
            "cluster_ranking": [1, 2],
            "cluster_id": [1, 2],
            "model-cluster_ranking": [1, 2],
            # Dont care about the rest of the columns, they are ignored by the function
        }
    )
    new_clusters_df, new_structs_df = create_other_cluster(clusters_df, structs_df, 2)
    pd.testing.assert_frame_equal(new_structs_df, structs_df)
    pd.testing.assert_frame_equal(new_clusters_df, clusters_df)


def test_create_other_cluster_maxover():
    clusters_df = pd.DataFrame(
        {
            "cluster_id": [1, 2, 3],
            "cluster_rank": [1, 2, 3],
            "n": [1, 1, 2],
            "caprieval_rank": [1, 2, 3],
            "score": [1.0, 2.0, 3.0],
            "score_std": [0.0, 0.0, 0.0],
            # Dont care about the rest of the columns, they are ignored by the function
        }
    )
    structs_df = pd.DataFrame(
        {
            "caprieval_rank": [1, 2, 3, 4],
            "cluster_ranking": [1, 2, 3, 3],
            "cluster_id": [1, 2, 3, 3],
            "model-cluster_ranking": [1, 1, 1, 2],
            "score": [1.0, 2.0, 3.0, 4.0],
            # Dont care about the rest of the columns, they are ignored by the function
        }
    )
    new_clusters_df, new_structs_df = create_other_cluster(clusters_df, structs_df, 2)
    expected_clusters_df = pd.DataFrame(
        {
            "cluster_id": [1, "Other"],
            "cluster_rank": [1, 2],
            "n": [1, 3],
            "caprieval_rank": [1, 2],
            "score": [1.0, 3.0],
            "score_std": [0.0, 1.0],
            # Dont care about the rest of the columns, they are ignored by the function
        }
    )
    expected_structs_df = pd.DataFrame(
        {
            "caprieval_rank": [1, 2, 3, 4],
            "cluster_ranking": [1, 2, 2, 2],
            "cluster_id": [1, "Other", "Other", "Other"],
            "model-cluster_ranking": [1, 1, 2, 3],
            "score": [1.0, 2.0, 3.0, 4.0],
            # Dont care about the rest of the columns, they are ignored by the function
        }
    )
    pd.testing.assert_frame_equal(new_structs_df, expected_structs_df)
    pd.testing.assert_frame_equal(new_clusters_df, expected_clusters_df)


@pytest.fixture
def example_capri_ss_dashcluster():
    """Provide example capri_ss.tsv filename."""
    return Path(data_folder, "capri_ss_-cluster.tsv")

# TODO find a way to not duplicate this fixture from test_module_alascan.py
# maybe move to conftest.py file?
@pytest.fixture
def example_df_scan_clt():
    """Return example alascan clt DataFrame."""
    example_clt_data = [
        ["A", 38, "ASP", "A-38-ASP", -2.0, -1.0, -0.4, -2.3, -0.5, -7.2, 1.0],
        ["A", 69, "LYS", "A-38-ASP", -0.0, 1.0, -0.4, 0.8, 0.5, -7.2, 1.0],
    ]
    columns = [
        "chain",
        "resnum",
        "resname",
        "full_resname",
        "score",
        "delta_score",
        "delta_vdw",
        "delta_elec",
        "delta_desolv",
        "delta_bsa",
        "frac_pres",
    ]

    example_df_scan_clt = pd.DataFrame(example_clt_data, columns=columns)
    yield example_df_scan_clt


def test_make_alascan_plot(example_df_scan_clt):
    """Test make_alascan_plot."""
    make_alascan_plot(example_df_scan_clt, clt_id="-")
    # assert existence of plot
    assert Path("scan_clt_-.html").exists()


def test_plotly_cdn_url():
    """Test to obtain plotly cdn url function."""
    from plotly.io._utils import plotly_cdn_url
    plotly_cdn_full_url = plotly_cdn_url()
    assert type(plotly_cdn_full_url) == str


def test_plotly_offline():
    """Test to obtain plotly self contained javascript."""
    from plotly.offline.offline import get_plotlyjs
    plotly_self_contained = get_plotlyjs()
    assert type(plotly_self_contained) == str


def test_offline_js_manager():
    """Test offline manager behavior."""
    # Offline == True
    with tempfile.TemporaryDirectory(".") as tmpdir:
        figurepath = Path(tmpdir, "figure.html")
        offline_plotly_js = offline_js_manager(figurepath, offline=True)
        assert "plotly_bundle.js" in offline_plotly_js
        assert Path(tmpdir, "plotly_bundle.js").exists()
    # Offline == False
    from plotly.io._utils import plotly_cdn_url
    plotly_cdn_full_url = plotly_cdn_url()
    offline_plotly_js = offline_js_manager(tmpdir, offline=False)
    assert plotly_cdn_full_url in offline_plotly_js


def test_box_plot_handler(example_capri_ss, cluster_ranking):
    """Test box plot generation without format definition."""
    initdir = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        box_plot_handler(
            example_capri_ss.resolve(),
            cluster_ranking,
            None,  # When not set (default), goes for html file generation
            1.0,
            )
        assert len(list(Path(".").glob("*.html"))) > 0
        assert len(list(Path(".").glob("*.png"))) == 0
    os.chdir(initdir)


def test_box_plot_handler_format(example_capri_ss, cluster_ranking):
    """Test box plot generation with png format definition."""
    initdir = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        box_plot_handler(
            example_capri_ss.resolve(),
            cluster_ranking,
            "png",
            1.0,
            )
        assert len(list(Path(".").glob("*.html"))) > 0
        assert len(list(Path(".").glob("*.png"))) > 0
    os.chdir(initdir)


def test_scatter_plot_handler(example_capri_ss, cluster_ranking):
    """Test scatter plot generation without format definition."""
    initdir = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        scatter_plot_handler(
            example_capri_ss.resolve(),
            cluster_ranking,
            None,  # When not set (default), goes for html file generation
            1.0,
            )
        assert len(list(Path(".").glob("*.html"))) > 0
        assert len(list(Path(".").glob("*.png"))) == 0
    os.chdir(initdir)


def test_scatter_plot_handler_format(example_capri_ss, cluster_ranking):
    """Test scatter plot generation with png format definition."""
    initdir = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        scatter_plot_handler(
            example_capri_ss.resolve(),
            cluster_ranking,
            "png",
            1.0,
            )
        assert len(list(Path(".").glob("*.html"))) > 0
        assert len(list(Path(".").glob("*.png"))) > 0
    os.chdir(initdir)
