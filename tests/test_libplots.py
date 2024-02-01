"""Test finding the best structures in libplots."""
from pathlib import Path

import pandas as pd
import pytest

from haddock.libs.libplots import find_best_struct, make_alascan_plot

from . import data_folder, golden_data
from .test_module_alascan import example_df_scan_clt


@pytest.fixture
def example_capri_ss():
    """Provide example capri_ss.tsv filename."""
    return Path(golden_data, "capri_ss_example.tsv")


def test_find_best_struct_best1(example_capri_ss):
    """Finds one best structure."""
    result = find_best_struct(example_capri_ss, 1)
    cluster1 = result.loc[result["Cluster ID"] == 1]

    expected = pd.DataFrame(
        {
            "Cluster ID": {0: 1},
            "Cluster Rank": {0: 2},
            "Nr 01 best structure": {
                0: "../../1_rigidbody/rigidbody_383.pdb",
                },
            }
        )
    expected.columns.names = ["Structure"]
    pd.testing.assert_frame_equal(cluster1, expected)


def test_find_best_struct_best20(example_capri_ss):
    """Finds 20 best structures if possible."""
    result = find_best_struct(example_capri_ss, 20)
    cluster1 = result.loc[result["Cluster ID"] == 1]

    expected = pd.DataFrame(
        {
            "Cluster ID": {0: 1},
            "Cluster Rank": {0: 2},
            "Nr 01 best structure": {0: "../../1_rigidbody/rigidbody_383.pdb"},
            "Nr 02 best structure": {0: "../../1_rigidbody/rigidbody_265.pdb"},
            "Nr 03 best structure": {0: "../../1_rigidbody/rigidbody_231.pdb"},
            "Nr 04 best structure": {0: "../../1_rigidbody/rigidbody_218.pdb"},
            }
        )
    expected.columns.names = ["Structure"]
    pd.testing.assert_frame_equal(cluster1, expected)


def test_find_best_struct_best4(example_capri_ss):
    """Finds 4 best structures if possible."""
    result = find_best_struct(example_capri_ss, 4)
    cluster1 = result.loc[result["Cluster ID"] == 1]

    expected = pd.DataFrame(
        {
            "Cluster ID": {0: 1},
            "Cluster Rank": {0: 2},
            "Nr 01 best structure": {0: "../../1_rigidbody/rigidbody_383.pdb"},
            "Nr 02 best structure": {0: "../../1_rigidbody/rigidbody_265.pdb"},
            "Nr 03 best structure": {0: "../../1_rigidbody/rigidbody_231.pdb"},
            "Nr 04 best structure": {0: "../../1_rigidbody/rigidbody_218.pdb"},
            }
        )
    expected.columns.names = ["Structure"]
    pd.testing.assert_frame_equal(cluster1, expected)


@pytest.fixture
def example_capri_ss_dashcluster():
    """Provide example capri_ss.tsv filename."""
    return Path(data_folder, "capri_ss_-cluster.tsv")


def test_find_best_struct_unclustered(example_capri_ss_dashcluster):
    """Finds 4 best structures when there unclustered."""
    result = find_best_struct(example_capri_ss_dashcluster, 4)

    struct_dir = "../../../01_rigidbody"
    expected = pd.DataFrame(
        {
            "Cluster ID": {0: "-"},
            "Cluster Rank": {0: "Unclustered"},
            "Nr 01 best structure": {0: f"{struct_dir}/rigidbody_6.pdb"},
            "Nr 02 best structure": {0: f"{struct_dir}/rigidbody_16.pdb"},
            "Nr 03 best structure": {0: f"{struct_dir}/rigidbody_20.pdb"},
            "Nr 04 best structure": {0: f"{struct_dir}/rigidbody_14.pdb"},
            }
        )
    expected.columns.names = ["Structure"]
    pd.testing.assert_frame_equal(result, expected)


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