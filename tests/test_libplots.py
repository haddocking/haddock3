"""Test finding the best structures in libplots."""

from pathlib import Path
import numpy as np

import pandas as pd
import pytest

from haddock.libs.libplots import find_best_struct, make_alascan_plot, read_capri_table

from . import data_folder, golden_data


@pytest.fixture
def example_capri_ss():
    """Provide example capri_ss.tsv df."""
    return read_capri_table(Path(golden_data, "capri_ss_example.tsv"))


def test_find_best_struct_best1(example_capri_ss):
    """Finds one best structure."""
    result = find_best_struct(example_capri_ss, 1)
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


def test_find_best_struct_best5_withmssing(example_capri_ss):
    """Finds 20 best structures if possible."""
    result = find_best_struct(example_capri_ss, 5)
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


def test_find_best_struct_bestdefault(example_capri_ss):
    """Finds 4 best structures if possible."""
    result = find_best_struct(example_capri_ss)
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
            "cluster-id": [1, 1, 1],
            "model": [
                "../1_rigidbody/rigidbody_383.pdb",
                "../1_rigidbody/rigidbody_265.pdb",
                "../1_rigidbody/rigidbody_231.pdb",
            ],
            "model-cluster-ranking": [1, 2, 3],
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


@pytest.fixture
def example_capri_ss_dashcluster():
    """Provide example capri_ss.tsv filename."""
    return Path(data_folder, "capri_ss_-cluster.tsv")


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
