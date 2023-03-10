from pathlib import Path

import pytest
import pandas as pd

from haddock.libs.libplots import find_best_struct
from . import golden_data, data_folder


@pytest.fixture
def example_capri_ss():
    """Provide example capri_ss.tsv filename."""
    return Path(golden_data, "capri_ss_example.tsv")


def test_find_best_struct_best1(example_capri_ss):
    result = find_best_struct(example_capri_ss, 1)
    cluster1 = result.loc[result["Cluster ID"] == 1]

    expected = pd.DataFrame(
        {
            "Cluster ID": {0: 1},
            "Nr 01 best structure": {
                0: "../1_rigidbody/rigidbody_383.pdb",
            },
        }
    )
    expected.columns.names = ["Structure"]
    pd.testing.assert_frame_equal(cluster1, expected)


def test_find_best_struct_best20(example_capri_ss):
    result = find_best_struct(example_capri_ss, 20)
    cluster1 = result.loc[result["Cluster ID"] == 1]

    expected = pd.DataFrame(
        {
            "Cluster ID": {0: 1},
            "Nr 01 best structure": {0: "../1_rigidbody/rigidbody_383.pdb"},
            "Nr 02 best structure": {0: "../1_rigidbody/rigidbody_265.pdb"},
            "Nr 03 best structure": {0: "../1_rigidbody/rigidbody_231.pdb"},
            "Nr 04 best structure": {0: "../1_rigidbody/rigidbody_218.pdb"},
        }
    )
    expected.columns.names = ["Structure"]
    pd.testing.assert_frame_equal(cluster1, expected)


def test_find_best_struct_best10(example_capri_ss):
    result = find_best_struct(example_capri_ss, 4)
    cluster1 = result.loc[result["Cluster ID"] == 1]

    expected = pd.DataFrame(
        {
            "Cluster ID": {0: 1},
            "Nr 01 best structure": {0: "../1_rigidbody/rigidbody_383.pdb"},
            "Nr 02 best structure": {0: "../1_rigidbody/rigidbody_265.pdb"},
            "Nr 03 best structure": {0: "../1_rigidbody/rigidbody_231.pdb"},
            "Nr 04 best structure": {0: "../1_rigidbody/rigidbody_218.pdb"},
        }
    )
    expected.columns.names = ["Structure"]
    pd.testing.assert_frame_equal(cluster1, expected)


@pytest.fixture
def example_capri_ss_dashcluster():
    """Provide example capri_ss.tsv filename."""
    return Path(data_folder, "capri_ss_-cluster.tsv")


def test_find_best_struct_unclustered(example_capri_ss_dashcluster):
    result = find_best_struct(example_capri_ss_dashcluster, 4)

    expected = pd.DataFrame(
        {
            "Cluster ID": {0: "-"},
            "Nr 01 best structure": {0: "../../01_rigidbody/rigidbody_6.pdb"},
            "Nr 02 best structure": {0: "../../01_rigidbody/rigidbody_16.pdb"},
            "Nr 03 best structure": {0: "../../01_rigidbody/rigidbody_20.pdb"},
            "Nr 04 best structure": {0: "../../01_rigidbody/rigidbody_14.pdb"},
        }
    )
    expected.columns.names = ["Structure"]
    pd.testing.assert_frame_equal(result, expected)
