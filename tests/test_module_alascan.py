"""Test the Alascan module."""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import tempfile
import shutil
from unittest.mock import patch
from haddock.modules.analysis.alascan import HaddockModule as AlascanModule
from haddock.modules.analysis.alascan import DEFAULT_CONFIG

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.alascan.scan import (
    Scan,
    ScanJob,
    add_delta_to_bfactor,
    add_zscores,
    check_alascan_jobs,
    mutate
)

from . import golden_data


@pytest.fixture
def complex_pdb():
    """Return example complex pdb."""
    return Path(golden_data, "protprot_complex_1.pdb")


@pytest.fixture
def protprot_model_list(complex_pdb):
    """Prot-prot input."""
    return [
        PDBFile(file_name=complex_pdb, path=str(golden_data)),
    ]


@pytest.fixture
def params():
    return {"int_cutoff": 3.0, "plot": False, "scan_residue": "ALA"}


@pytest.fixture
def scan_obj(protprot_model_list, params):
    """Return example alascan module."""
    scan_obj = Scan(
        model_list=protprot_model_list,
        output_name="alascan",
        path=Path("."),
        core=0,
        params=params,
    )

    yield scan_obj


@pytest.fixture
def scanjob_obj(scan_obj):
    """Return example alascan module."""
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        yield ScanJob(
        scan_obj=scan_obj,
        output=Path("alascan"),
        params=params,
    )

    #scanjob_obj = ScanJob(
    #    scan_obj=scan_obj,
    #    output=Path("alascan"),
    #    params=params,
    #)

    #yield scanjob_obj


@pytest.fixture
def alascan():
    """Return alascan module."""

    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        yield AlascanModule(
            order=1,
            path=Path(tmpdir),
            initial_params=DEFAULT_CONFIG,
        )


@pytest.fixture
def example_df_scan():
    """Return example alascan DataFrame."""
    example_scan_data = [
        [
            "A",
            19,
            "THR",
            "ALA",
            -108.540,
            -43.234,
            -275.271,
            -10.252,
            1589.260,
            -0.243,
            5.401,
            0.119,
            28.482,
            -0.414,
            10.470,
            -0.329,
        ],
        [
            "A",
            20,
            "ILE",
            "ALA",
            -18.540,
            -41.234,
            -270.271,
            -11.252,
            589.260,
            +0.243,
            4.401,
            0.109,
            27.482,
            -0.3,
            0.470,
            +0.329,
        ],  # this one is truly random
    ]
    columns = [
        "chain",
        "res",
        "ori_resname",
        "end_resname",
        "score",
        "vdw",
        "elec",
        "desolv",
        "bsa",
        "delta_ori_score",
        "delta_score",
        "delta_vdw",
        "delta_elec",
        "delta_desolv",
        "delta_bsa",
        "z_score",
    ]
    example_df_scan = pd.DataFrame(example_scan_data, columns=columns)

    yield example_df_scan


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


def test_scan_obj(scan_obj, protprot_model_list):
    """Test the alascan module."""
    assert scan_obj.int_cutoff == 3.0
    assert scan_obj.scan_res == "ALA"


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


@pytest.mark.skip("Not implemented yet")
def test_get_index_list():
    assert False


def test_confirm_installation(alascan):
    assert alascan.confirm_installation() is None


def test_init(alascan):
    alascan.__init__(
        order=42,
        path=Path("0_anything"),
        initial_params=DEFAULT_CONFIG,
    )

    # Once a module is initialized, it should have the following attributes
    assert alascan.path == Path("0_anything")
    assert alascan._origignal_config_file == DEFAULT_CONFIG
    assert type(alascan.params) == dict
    assert len(alascan.params) != 0


@patch("haddock.modules.analysis.alascan.Scheduler")
@patch("haddock.modules.analysis.alascan.alascan_cluster_analysis")
@patch("haddock.modules.analysis.alascan.HaddockModule.export_io_models")
def test_run(
    MockScheduler,
    mock_alascan_cluster_analysis,
    mock_export_io_models,
    alascan,
):
    # Only add files that are used within the body of the `run` method!
    Path(alascan.path, "alascan_0.scan").touch()
    ala_path = Path(golden_data, "scan_protprot_complex_1.csv")
    shutil.copy(ala_path, Path(alascan.path, "scan_protprot_complex_1.csv"))
    
    alascan.params["plot"] = True
    alascan.params["output"] = True

    alascan.previous_io = MockPreviousIO()
    alascan.run()


class MockPreviousIO:
    # In the mocked method, add the arguments that are called by the original method
    #  that is being tested
    def retrieve_models(self, individualize: bool = False):
        return [
            PDBFile(file_name="protprot_complex_1.pdb", path=str(golden_data)),
        ]


def mock_alascan_cluster_analysis():
    # Return whatever is necessary for the `run` method to work
    return {
        0: {
            "X-X-X": {
                "delta_score": 0.0,
                "delta_vdw": 0.0,
                "delta_elec": 0.0,
                "delta_desolv": 0.0,
                "delta_bsa": 0.0,
                "frac_pr": 0.0,
            }
        },
    }


def test_check_alascan_jobs(scanjob_obj):
    # calling check_alascan_jobs should raise an Exception as the job 
    # has not run yet
    with pytest.raises(Exception):
        check_alascan_jobs([scanjob_obj])
    
    
