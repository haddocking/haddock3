"""Test the Alascan module."""

import os
import time
import signal
import shutil
import tempfile
import multiprocessing
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.alascan import DEFAULT_CONFIG
from haddock.libs.libparallel import Scheduler, GenericTask
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
    _SHUTDOWN_REQUESTED, 
    _handle_shutdown_signal,
    calculate_core_allocation,
    process_residue_task,
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
        "chains": [],
        "output_mutants": False,
    }


@pytest.fixture(name="scan_obj")
def fixture_scan_obj(protprot_model_list, params, monkeypatch):
    """Return example alascan module."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        yield Scan(
            model_list=protprot_model_list,
            output_name="alascan",
            path=Path("."),
            core=0,
            residue_ncores=1,
            params=params,
        )


@pytest.fixture(name="scanjob_obj")
def fixture_scanjob_obj(scan_obj, params, monkeypatch):
    """Return example alascan module."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        yield ScanJob(
            scan_obj=scan_obj,
            params=params,
        )


@pytest.fixture(name="alascan")
def fixture_alascan(monkeypatch):
    """Return alascan module."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
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
    yield Path(golden_data, "scan_protprot_complex_1.tsv")


def test_scan_obj(scan_obj, protprot_model_list):
    """Test the alascan module."""
    assert scan_obj.int_cutoff == 3.0
    assert scan_obj.scan_res == "ALA"


def test_mutate(protprot_model_list, monkeypatch):
    """Test the mutate function."""
    with tempfile.TemporaryDirectory() as tmpdir:
        mut_fname = Path(golden_data, protprot_model_list[0].file_name)
        monkeypatch.chdir(path=tmpdir)
        mut_pdb_fname = mutate(mut_fname, "A", 19, "ALA")

        assert mut_pdb_fname == Path("protprot_complex_1-A_T19A.pdb")
        assert os.path.exists(mut_pdb_fname)

        first_line = open(mut_pdb_fname).readline()
        assert first_line[17:20] == "ALA"

        with pytest.raises(KeyError):
            mutate(mut_fname, "A", 19, "HOH")


def test_add_delta_to_bfactor(protprot_model_list, example_df_scan, monkeypatch):
    """Test the add_delta_to_bfactor function."""
    with tempfile.TemporaryDirectory() as tmpdir:
        mut_fname = Path(golden_data, protprot_model_list[0].file_name)
        monkeypatch.chdir(path=tmpdir)
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
    ala_path = Path(golden_data, "scan_protprot_complex_1.tsv")
    shutil.copy(ala_path, Path(alascan.path, "scan_protprot_complex_1.tsv"))

    alascan.params["plot"] = True

    alascan.previous_io = MockPreviousIO()
    mocker.patch("haddock.libs.libparallel.Scheduler.run", return_value=None)
    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models", return_value=None
    )

    alascan.run()
    assert Path(alascan.path, "scan_clt_-.tsv").exists()
    assert Path(alascan.path, "scan_clt_-.html").exists()


def test_scanjob_run(scanjob_obj, mocker):
    mocker.patch("haddock.modules.analysis.alascan.scan.Scan.run", return_value=None)
    scanjob_obj.run()


def test_scan_run_output(mocker, scan_obj):
    """Test Scan run and output method."""
    mocker.patch(
        "haddock.modules.analysis.alascan.scan.calc_score",
        return_value=(-106.7, -29, -316, -13, 1494),
    )
    scan_obj.run()
    assert Path(scan_obj.path, "scan_protprot_complex_1.tsv").exists()
    assert scan_obj.df_scan.shape[0] == 2


def test_scan_run_interface(mocker, scan_obj):
    """Test Scan run with empty filter_resdic."""
    scan_obj.filter_resdic = {"_": []}
    scan_obj.scan_res = "ASP"
    mocker.patch(
        "haddock.modules.analysis.alascan.scan.calc_score",
        return_value=(-106.7, -29, -316, -13, 1494),
    )
    scan_obj.run()

    assert Path(scan_obj.path, "scan_protprot_complex_1.tsv").exists()
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


def test_generate_alascan_output(mocker, protprot_model_list, scan_file, monkeypatch):
    """Test the generate_alascan_output method."""
    with tempfile.TemporaryDirectory() as tmpdir:
        shutil.copy(scan_file, Path(tmpdir, "scan_protprot_complex_1.tsv"))
        monkeypatch.chdir(tmpdir)
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


def test_alascan_cluster_analysis(protprot_input_list, scan_file, monkeypatch):
    """Test alascan_cluster_analysis."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        shutil.copy(scan_file, Path("scan_protprot_complex_1.tsv"))
        shutil.copy(scan_file, Path("scan_protprot_complex_2.tsv"))
        alascan_cluster_analysis(protprot_input_list)

        assert Path("scan_clt_-.tsv").exists()

        protprot_input_list[1].clt_id = 1
        alascan_cluster_analysis(protprot_input_list)

        assert Path("scan_clt_1.tsv").exists()
        assert Path("scan_clt_-.tsv").exists()


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


# tests related to per-residue parallelization within alascan
@pytest.fixture(name="scan_obj_parallel")
def fixture_scan_obj_parallel(protprot_model_list, params, monkeypatch):
    """Return example alascan module configured for parallel processing."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        # Configure for parallel processing
        params_parallel = params.copy()
        yield Scan(
            model_list=protprot_model_list,
            path=Path("."),
            core=0,
            residue_ncores=4,  # Multiple cores for parallel processing
            params=params_parallel,
        )


def test_calculate_core_allocation_more_models_than_cores():
    """Test core allocation when we have more models than cores."""
    model_cores, residue_cores = calculate_core_allocation(nmodels=10, total_cores=4)
    assert model_cores == 4
    assert residue_cores == 1


def test_calculate_core_allocation_less_models_than_cores():
    """Test core allocation when we have fewer models than cores."""
    model_cores, residue_cores = calculate_core_allocation(nmodels=2, total_cores=8)
    assert model_cores == 2
    assert residue_cores == 4


def test_calculate_core_allocation_equal_models_and_cores():
    """Test core allocation when models equal cores."""
    model_cores, residue_cores = calculate_core_allocation(nmodels=4, total_cores=4)
    assert model_cores == 4
    assert residue_cores == 1


def test_calculate_core_allocation_single_model():
    """Test core allocation for single model case."""
    model_cores, residue_cores = calculate_core_allocation(nmodels=1, total_cores=4)
    assert model_cores == 1
    assert residue_cores == 4


def test_process_residue_task(mocker, monkeypatch):
    """Test normal execution of process_residue_task."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        mocker.patch(
            "haddock.modules.analysis.alascan.scan.mutate",
            return_value=Path("mutated.pdb")
        )
        mocker.patch(
            "haddock.modules.analysis.alascan.scan.calc_score",
            return_value=(-100.0, -30.0, -300.0, -10.0, 1000.0)
        )
        mocker.patch("os.path.exists", return_value=True)
        mocker.patch("os.remove", return_value=None)
        args = (
            "A", 19, {"A-19": "THR"}, Path("native.pdb"), "ALA",
            -105.0, -30.5, -310.0, -10.5, 1100.0,
            "score_dir", False
        )
        result = process_residue_task(args)
        expected = [
            "A", 19, "THR", "ALA",
            -100.0, -30.0, -300.0, -10.0, 1000.0,
            -5.0, -0.5, -10.0, -0.5, 100.0
        ]
        assert result == expected


def test_process_residue_task_shutdown(mocker):
    """Test process_residue_task handles respond shutdown."""
    mocker.patch('haddock.modules.analysis.alascan.scan._SHUTDOWN_REQUESTED', True)
    args = ("A", 19, {"A-19": "THR"}, Path("native.pdb"), "ALA",
        -105.0, -30.5, -310.0, -10.5, 1100.0,
        "score_dir", False)
    result = process_residue_task(args)
    assert result is None


def test_process_residue_task_exception(mocker):
    mock_log = mocker.patch("haddock.modules.analysis.alascan.scan.log")
    mocker.patch(
        "haddock.modules.analysis.alascan.scan.mutate",
        side_effect=Exception("Test exception"))
    args = (
        "A", 19, {"A-19": "THR"}, Path("native.pdb"), "ALA",
        -105.0, -30.5, -310.0, -10.5, 1100.0,
        "score_dir", False)
    result = process_residue_task(args)
    assert result is None
    mock_log.warning.assert_called_with(mocker.ANY)
    assert "Error processing A-19" in str(mock_log.warning.call_args)


def test_process_residue_task_keyboard_interrupt(mocker):
    """Test process_residue_task handles KeyboardInterrupt directly."""
    mocker.patch(
        "haddock.modules.analysis.alascan.scan.mutate",
        side_effect=KeyboardInterrupt("Direct interrupt")
    )
    args = (
        "A", 19, {"A-19": "THR"}, Path("native.pdb"), "ALA",
        -105.0, -30.5, -310.0, -10.5, 1100.0,
        "score_dir", False
    )
    result = process_residue_task(args)
    assert result is None


def test_run_sequential_residues(mocker, scan_obj_parallel):
    """Test sequential residue processing."""
    mock_process = mocker.patch(
        "haddock.modules.analysis.alascan.scan.process_residue_task",
        side_effect=[
            ["A", 19, "THR", "ALA", -100.0, -30.0, -300.0, -10.0, 1000.0, -5.0, -0.5, -10.0, -0.5, 100.0],
            None,  # failed task
            ["A", 21, "VAL", "ALA", -102.0, -32.0, -302.0, -11.0, 1002.0, -3.0, 1.5, -8.0, 0.5, 98.0]
        ])
    tasks = [
        ("A", 19, {"A-19": "THR"}, Path("native.pdb"), "ALA", -105.0, -30.5, -310.0, -10.5, 1100.0, "score_dir", False),
        ("A", 20, {"A-20": "ILE"}, Path("native.pdb"), "ALA", -101.0, -28.5, -306.0, -9.5, 1098.0, "score_dir", False),
        ("A", 21, {"A-21": "VAL"}, Path("native.pdb"), "ALA", -105.0, -33.5, -310.0, -11.5, 1100.0, "score_dir", False)]
    results = scan_obj_parallel._run_sequential_residues(tasks)
    assert len(results) == 2  # One failed task should be filtered out
    assert mock_process.call_count == 3


def test_run_parallel_residues(mocker, scan_obj_parallel):
    """Test successful parallel residue processing."""
    mock_process_results = [
        ["A", 19, "THR", "ALA", -100.0, -30.0, -300.0, -10.0, 1000.0, -5.0, -0.5, -10.0, -0.5, 100.0],
        ["A", 20, "ILE", "ALA", -98.0, -28.0, -298.0, -9.0, 998.0, -3.0, 0.5, -8.0, 0.5, 98.0]
    ]
    mocker.patch("haddock.modules.analysis.alascan.scan.process_residue_task", 
                 side_effect=mock_process_results)
    mock_cleanup = mocker.patch.object(scan_obj_parallel, "_cleanup_scheduler_resources")
    
    tasks = [
        ("A", 19, {"A-19": "THR"}, Path("native.pdb"), "ALA", -105.0, -30.5, -310.0, -10.5, 1100.0, "score_dir", False),
        ("A", 20, {"A-20": "ILE"}, Path("native.pdb"), "ALA", -101.0, -28.5, -306.0, -9.5, 1098.0, "score_dir", False)
    ]
    
    results = scan_obj_parallel._run_parallel_residues(tasks)
    
    assert len(results) == 2
    assert results[0][0] == "A"
    assert results[0][1] == 19
    mock_cleanup.assert_called_once()


def test_run_parallel_residues_exception(mocker, scan_obj_parallel, caplog):
    """Test parallel processing falls back to sequential on exception."""
    mock_scheduler = mocker.Mock()
    mock_scheduler.run.side_effect = Exception("Parallel execution failed")
    mocker.patch("haddock.modules.analysis.alascan.scan.Scheduler", return_value=mock_scheduler)
    mock_cleanup = mocker.patch.object(scan_obj_parallel, "_cleanup_scheduler_resources")
    # Mock process_residue_task to return valid data 
    expected = ["A", 19, "THR", "ALA", -100.0, -30.0, -300.0, -10.0, 1000.0, -5.0, -0.5, -10.0, -0.5, 100.0]
    mocker.patch("haddock.modules.analysis.alascan.scan.process_residue_task", return_value=expected)
    tasks = [
        ("A", 19, {"A-19": "THR"}, Path("native.pdb"), "ALA", 
         -105.0, -30.5, -310.0, -10.5, 1100.0, "score_dir", False)
    ]    
    result = scan_obj_parallel._run_parallel_residues(tasks)
    mock_cleanup.assert_called_once_with(mock_scheduler)
    assert "Parallel execution failed" in caplog.text
    assert "Falling back to sequential" in caplog.text
    assert result == [expected]


def test_run_parallel_residues_keyboard_interrupt(mocker, scan_obj_parallel, caplog, capsys):
    """Test parallel processing handles KeyboardInterrupt properly."""
    # Mock the Scheduler to raise KeyboardInterrupt when run() is called
    #  Mock at the module level where it's imported
    mock_scheduler = mocker.Mock()
    mock_scheduler.run.side_effect = KeyboardInterrupt("User interrupted")
    mocker.patch("haddock.modules.analysis.alascan.scan.Scheduler", return_value=mock_scheduler)
    mocker.patch("haddock.modules.analysis.alascan.scan.GenericTask")
    mocker.patch("time.sleep")
    # Mock the cleanup method to verify it's called
    mock_cleanup = mocker.patch.object(scan_obj_parallel, "_cleanup_scheduler_resources")
    tasks = [
        ("A", 19, {"A-19": "THR"}, Path("native.pdb"), "ALA",
         -105.0, -30.5, -310.0, -10.5, 1100.0, "score_dir", False)
    ]   
    # Call the method - this should trigger KeyboardInterrupt handling
    # (because mock_scheduler.run.side_effect)
    result = scan_obj_parallel._run_parallel_residues(tasks)
    # Verify the expected flow:
    # 1. KeyboardInterrupt was caught and logged
    assert "Parallel execution interrupted" in caplog.text
    # 2. _handle_interrupt_cleanup was called
    captured = capsys.readouterr()
    assert "Hold tight, interrupting run..." in captured.out
    # 3. _cleanup_scheduler_resources was called twice 
    # by _handle_interrupt_cleanup and the "finally" block in _run_parallel_residues
    assert mock_cleanup.call_count == 2
    mock_cleanup.assert_called_with(mock_scheduler)
    # 4. Method returns None (from _handle_interrupt_cleanup)
    assert result is None


def test_cleanup_scheduler_resources(mocker, scan_obj_parallel):
    """Test scheduler cleanup."""
    # Pretend workers, schedule and queue processing tacks
    mock_worker1 = mocker.Mock()
    mock_worker1.is_alive.return_value = True
    mock_worker1.join.return_value = None
    mock_worker1.terminate.return_value = None
    mock_worker2 = mocker.Mock()
    mock_worker2.is_alive.return_value = False
    mock_worker2.join.return_value = None
    mock_worker2.terminate.return_value = None
    mock_queue = mocker.Mock()
    mock_queue.empty.side_effect = [False, False, True] 
    mock_queue.get_nowait.return_value = "dummy_item"
    mock_queue.close.return_value = None
    mock_queue.join_thread.return_value = None
    mock_scheduler = mocker.Mock()
    mock_scheduler.worker_list = [mock_worker1, mock_worker2]
    mock_scheduler.queue = mock_queue    
    # Call the cleanup method
    try:
        scan_obj_parallel._cleanup_scheduler_resources(mock_scheduler)
    except Exception:
        pass
    mock_worker1.join.assert_called()
    mock_worker2.join.assert_not_called()
    mock_worker1.terminate.assert_called()
    mock_worker2.terminate.assert_not_called()


def test_handle_interrupt_cleanup(mocker, scan_obj_parallel, capsys):
    """Test cleanup handling."""
    mock_scheduler = mocker.Mock()
    mock_cleanup = mocker.patch.object(scan_obj_parallel, "_cleanup_scheduler_resources")
    mocker.patch("time.sleep")
    tasks = []
    result = scan_obj_parallel._handle_interrupt_cleanup(mock_scheduler, tasks)
    mock_cleanup.assert_called_once_with(mock_scheduler)
    captured = capsys.readouterr()
    assert "Hold tight, interrupting run..." in captured.out
    assert result is None


def test_handle_shutdown_signal():
    """Test that shutdown is requested when _handle_shutdown_signal called."""
    import haddock.modules.analysis.alascan.scan as scan_module
    scan_module._SHUTDOWN_REQUESTED = False
    _handle_shutdown_signal(signal.SIGINT, None)
    # Check that shutdown was requested
    assert scan_module._SHUTDOWN_REQUESTED is True
