"""Test the Alascan module."""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.alascan.scan import (
    Scan,
    add_delta_to_bfactor,
    add_zscores,
    alascan_cluster_analysis,
    mutate,
    )

from . import golden_data


@pytest.fixture
def protprot_model_list():
    """Prot-prot input."""
    return [
        PDBFile(Path(golden_data, "protprot_complex_1.pdb"), path=golden_data),
        ]


@pytest.fixture
def params():
    return {
        "int_cutoff": 3.0, "plot": False, "scan_residue": "ALA"
        }


@pytest.fixture
def alascan_module(protprot_model_list, params):
    """Return example alascan module."""
    scan = Scan(
        model_list=protprot_model_list,
        output_name="alascan",
        path=Path("."),
        core=0,
        params=params
        )

    yield scan


@pytest.fixture
def example_df_scan():
    """Return example alascan DataFrame."""
    example_scan_data = [
        [
            "A", 19, "THR", "ALA", -108.540, -43.234,
            -275.271, -10.252, 1589.260, -0.243, 5.401,
            0.119, 28.482, -0.414, 10.470, -0.329
            ],
        [
            "A", 20, "ILE", "ALA", -18.540, -41.234,
            -270.271, -11.252, 589.260, +0.243, 4.401,
            0.109, 27.482, -0.3, 0.470, +0.329
            ],  # this one is truly random
        ]
    columns = ["chain", "res", "ori_resname", "end_resname", "score",
               "vdw", "elec", "desolv", "bsa", "delta_ori_score",
               "delta_score", "delta_vdw", "delta_elec", "delta_desolv",
               "delta_bsa", "z_score"]
    example_df_scan = pd.DataFrame(
        example_scan_data,
        columns=columns
        )

    yield example_df_scan


@pytest.fixture
def example_df_scan_clt():
    """Return example alascan clt DataFrame."""
    example_clt_data = [
        ["A", 38, "ASP", "A-38-ASP", -2.0, -1.0, -0.4, -2.3, -0.5, -7.2, 1.0],
        ["A", 69, "LYS", "A-38-ASP", -0.0, 1.0, -0.4, 0.8, 0.5, -7.2, 1.0],
        ]
    columns = ["chain", "resnum", "resname", "full_resname", "score",
               "delta_score", "delta_vdw", "delta_elec", "delta_desolv",
               "delta_bsa", "frac_pres"]

    example_df_scan_clt = pd.DataFrame(
        example_clt_data,
        columns=columns
        )
    yield example_df_scan_clt


def test_alascan_module(alascan_module, protprot_model_list):
    """Test the alascan module."""
    assert alascan_module.int_cutoff == 3.0
    assert alascan_module.scan_res == "ALA"

    # this requires CNS to be installed
    # alascan_module.run()
    # assert alascan_module.df_scan.empty

    # now test the alascan_cluster_analysis function
    clt_scan = alascan_cluster_analysis(protprot_model_list)
    assert clt_scan == {"-": {}}
    os.remove("scan_protprot_complex_1.pdb.csv")
    os.remove("scan_clt_-.csv")


def test_mutate(protprot_model_list):
    """Test the mutate function."""
    mut_fname = Path(golden_data, protprot_model_list[0].file_name)
    mut_pdb_fname = mutate(mut_fname, "A", 19, "ALA")
    assert mut_pdb_fname == Path("protprot_complex_1-A_T19A.pdb")
    assert os.path.exists(mut_pdb_fname)

    first_line = open(mut_pdb_fname).readline()
    assert first_line[17:20] == "ALA"
    # clean up
    os.remove(mut_pdb_fname)


def test_add_delta_to_bfactor(protprot_model_list, example_df_scan):
    """Test the add_delta_to_bfactor function."""
    mut_fname = Path(golden_data, protprot_model_list[0].file_name)
    mut_pdb_fname = mutate(mut_fname, "A", 19, "ALA")
    add_delta_to_bfactor(str(mut_pdb_fname), example_df_scan)
    # ALA 19 should have beta = 100.0
    assert np.isclose(100.0, float(open(mut_pdb_fname).readline()[60:66]))

    # clean up
    os.remove(mut_pdb_fname)


def test_add_zscores(example_df_scan_clt):
    """Test the add_zscores function."""
    obs_df_scan_clt = add_zscores(example_df_scan_clt)
    assert np.isclose(obs_df_scan_clt["z_score"].values[0], -1.0)
    assert np.isclose(obs_df_scan_clt["z_score"].values[1], 1.0)
