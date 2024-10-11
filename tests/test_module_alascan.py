"""Test the Alascan module."""

import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.alascan import DEFAULT_CONFIG
from haddock.modules.analysis.alascan import HaddockModule as AlascanModule
from haddock.modules.analysis.alascan.scan import (
    Scan,
    ScanJob,
    add_delta_to_bfactor,
    add_zscores,
    alascan_cluster_analysis,
    calc_score,
    create_alascan_plots,
    generate_alascan_output,
    mutate,
    )

from . import golden_data


@pytest.fixture(name="protprot_input_list")
def fixture_protprot_input_list():
    """Prot-prot input."""
    return [
        PDBFile(Path(golden_data, "protprot_complex_1.pdb"), path=golden_data),
        PDBFile(Path(golden_data, "protprot_complex_2.pdb"), path=golden_data),
    ]


@pytest.fixture(name="complex_pdb")
def fixture_complex_pdb():
    """Return example complex pdb."""
    return Path(golden_data, "protprot_complex_1.pdb")


@pytest.fixture(name="protprot_model_list")
def fixture_protprot_model_list(complex_pdb):
    """Prot-prot input."""
    return [
        PDBFile(file_name=complex_pdb, path=str(golden_data)),
    ]


@pytest.fixture(name="params")
def fixture_params():
    """Parameter fixture to be used in the tests"""
    return {
        "int_cutoff": 3.0,
        "plot": False,
        "scan_residue": "ALA",
        "resdic_A": [19, 20],
    }


@pytest.fixture(name="scan_obj")
def fixture_scan_obj(protprot_model_list, params):
    """Return example alascan module."""
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        yield Scan(
            model_list=protprot_model_list,
            output_name="alascan",
            path=Path("."),
            core=0,
            params=params,
        )


@pytest.fixture(name="scanjob_obj")
def fixture_scanjob_obj(scan_obj, params):
    """Return example alascan module."""
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        yield ScanJob(
            scan_obj=scan_obj,
            output=Path("."),
            params=params,
        )


@pytest.fixture(name="alascan")
def fixture_alascan():
    """Return alascan module."""
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        yield AlascanModule(
            order=1,
            path=Path("."),
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


@pytest.fixture(name="scan_file")
def fixture_scan_file():
    """Return example alascan file."""
    yield Path(golden_data, "scan_protprot_complex_1.csv")


def test_scan_obj(scan_obj, protprot_model_list):
    """Test the alascan module."""
    assert scan_obj.int_cutoff == 3.0
    assert scan_obj.scan_res == "ALA"


def test_mutate(protprot_model_list):
    """Test the mutate function."""
    with tempfile.TemporaryDirectory() as tmpdir:
        mut_fname = Path(golden_data, protprot_model_list[0].file_name)
        os.chdir(path=tmpdir)
        mut_pdb_fname = mutate(mut_fname, "A", 19, "ALA")

        assert mut_pdb_fname == Path("protprot_complex_1-A_T19A.pdb")
        assert os.path.exists(mut_pdb_fname)

        first_line = open(mut_pdb_fname).readline()
        assert first_line[17:20] == "ALA"

        with pytest.raises(KeyError):
            mutate(mut_fname, "A", 19, "HOH")


def test_add_delta_to_bfactor(protprot_model_list, example_df_scan):
    """Test the add_delta_to_bfactor function."""
    with tempfile.TemporaryDirectory() as tmpdir:
        mut_fname = Path(golden_data, protprot_model_list[0].file_name)
        os.chdir(path=tmpdir)
        mut_pdb_fname = mutate(mut_fname, "A", 19, "ALA")
        add_delta_to_bfactor(str(mut_pdb_fname), example_df_scan)
        # ALA 19 should have beta = 100.0
        assert np.isclose(100.0, float(open(mut_pdb_fname).readline()[60:66]))


def test_add_zscores(example_df_scan_clt):
    """Test the add_zscores function."""
    obs_df_scan_clt = add_zscores(example_df_scan_clt)
    assert np.isclose(obs_df_scan_clt["z_score"].values[0], -1.0)
    assert np.isclose(obs_df_scan_clt["z_score"].values[1], 1.0)


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


def test_run(
    mocker,
    alascan,
):
    # Only add files that are used within the body of the `run` method!

    Path(alascan.path, "alascan_0.scan").touch()
    ala_path = Path(golden_data, "scan_protprot_complex_1.csv")
    shutil.copy(ala_path, Path(alascan.path, "scan_protprot_complex_1.csv"))

    alascan.params["plot"] = True

    alascan.previous_io = MockPreviousIO()
    mocker.patch("haddock.libs.libparallel.Scheduler.run", return_value=None)
    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models", return_value=None
    )

    alascan.run()
    assert Path(alascan.path, "scan_clt_-.csv").exists()
    assert Path(alascan.path, "scan_clt_-.html").exists()

    alascan.params["output"] = True
    alascan.run()


def test_scanjob_run(scanjob_obj, mocker):
    mocker.patch("haddock.modules.analysis.alascan.scan.Scan.run", return_value=None)
    mocker.patch("haddock.modules.analysis.alascan.scan.Scan.output", return_value=None)
    scanjob_obj.run()


def test_scan_run_output(mocker, scan_obj):
    """Test Scan run and output method."""
    mocker.patch(
        "haddock.modules.analysis.alascan.scan.calc_score",
        return_value=(-106.7, -29, -316, -13, 1494),
    )
    scan_obj.run()
    assert Path(scan_obj.path, "scan_protprot_complex_1.csv").exists()
    assert scan_obj.df_scan.shape[0] == 2

    scan_obj.output()
    assert Path(scan_obj.path, "alascan").exists()


def test_scan_run_interface(mocker, scan_obj):
    """Test Scan run with empty filter_resdic."""
    scan_obj.filter_resdic = {"_": []}
    scan_obj.scan_res = "ASP"
    mocker.patch(
        "haddock.modules.analysis.alascan.scan.calc_score",
        return_value=(-106.7, -29, -316, -13, 1494),
    )
    scan_obj.run()

    assert Path(scan_obj.path, "scan_protprot_complex_1.csv").exists()
    assert scan_obj.df_scan.shape[0] == 5


def test_calc_score(mocker):
    """Test the run_scan method."""
    mocker.patch(
        "haddock.modules.analysis.alascan.scan.get_score_string",
        return_value=[
            "> starting calculations...",
            "> HADDOCK-score = (1.0 * vdw) + (0.2 * elec) + (1.0 * desolv) + (0.0 * air) + (0.0 * bsa)",
            "> HADDOCK-score (emscoring) = -106.7376",
            "> vdw=-29.5808,elec=-316.542,desolv=-13.8484,air=0.0,bsa=1494.73",
        ],
    )
    # sub with tmpdir
    scores = calc_score(Path(golden_data, "protprot_complex_1.pdb"), run_dir="tmp")

    assert scores == (-106.7376, -29.5808, -316.542, -13.8484, 1494.73)


def test_calc_score_wrong(mocker):
    mocker.patch(
        "haddock.modules.analysis.alascan.scan.get_score_string",
        return_value=["> could not calculate score"],
    )
    # now calc_score should raise an Exception
    with pytest.raises(Exception):
        calc_score(Path(golden_data, "protprot_complex_1.pdb"), run_dir="tmp")


def test_generate_alascan_output(mocker, protprot_model_list, scan_file):
    """Test the generate_alascan_output method."""
    with tempfile.TemporaryDirectory() as tmpdir:
        shutil.copy(scan_file, Path(tmpdir, "scan_protprot_complex_1.csv"))
        os.chdir(tmpdir)
        models_to_export = generate_alascan_output(protprot_model_list, path=tmpdir)
        assert len(models_to_export) == 1
        assert models_to_export[0].ori_name == "protprot_complex_1.pdb"
        assert models_to_export[0].file_name == "protprot_complex_1_alascan.pdb"


class MockPreviousIO:
    """Mock PreviousIO class."""

    # In the mocked method, add the arguments that are called by the original method
    #  that is being tested
    def retrieve_models(self, individualize: bool = False):
        return [
            PDBFile(file_name="protprot_complex_1.pdb", path=str(golden_data)),
        ]


def mock_alascan_cluster_analysis():
    """Mock alascan_cluster_analysis."""
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


def test_alascan_cluster_analysis(protprot_input_list, scan_file):
    """Test alascan_cluster_analysis."""
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        shutil.copy(scan_file, Path("scan_protprot_complex_1.csv"))
        shutil.copy(scan_file, Path("scan_protprot_complex_2.csv"))
        alascan_cluster_analysis(protprot_input_list)

        assert Path("scan_clt_-.csv").exists()

        protprot_input_list[1].clt_id = 1
        alascan_cluster_analysis(protprot_input_list)

        assert Path("scan_clt_1.csv").exists()
        assert Path("scan_clt_-.csv").exists()


def test_create_alascan_plots(mocker, caplog):
    """Test create_alascan_plots."""
    mocker.patch("pandas.read_csv", return_value=pd.DataFrame())
    create_alascan_plots({"-": []}, scan_residue="ALA")

    for record in caplog.records:
        assert record.levelname == "WARNING"

    # now assuming existing file but wrong plot creation
    mocker.patch("pandas.read_csv", return_value=[])
    mocker.patch("os.path.exists", return_value=True)

    create_alascan_plots({"-": []}, scan_residue="ALA")
    for record in caplog.records:
        assert record.levelname in ["INFO", "WARNING"]
