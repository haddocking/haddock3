"""Test the Alascan module."""

import os
import logging
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
    add_delta_to_bfactor,
    add_zscores,
    alascan_cluster_analysis,
    calc_score,
    create_alascan_plots,
    generate_alascan_output,
    mutate,
    write_scan_out,
    MutationResult,
    InterfaceScanner,
    ModelPointMutation,
    )

from . import golden_data

# fixures
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


@pytest.fixture(name="interface_scanner")
def fixture_interface_scanner(protprot_model_list, params):
    """Interface scanner in default HADDOCK mode."""
    return InterfaceScanner(
        model=protprot_model_list[0],
        mutation_res="ALA",
        params=params, 
        library_mode=False
    )


@pytest.fixture(name="interface_scanner_library")
def fixture_interface_scanner_library(complex_pdb, params):
    """Interface scanner in library mode."""
    return InterfaceScanner(
        model=complex_pdb,
        mutation_res="ALA", 
        params=params,
        library_mode=True
    )


@pytest.fixture(name="mutation_job")
def fixture_mutation_job(complex_pdb):
    """Create a ModelPointMutation job.""" 
    native_scores = (-106.7, -29.6, -316.5, -13.8, 1494.7)
    return ModelPointMutation(
        model_path=complex_pdb,
        model_id="test_model",
        chain="A", 
        res_num=19,
        ori_resname="THR",
        target_resname="ALA",
        native_scores=native_scores,
        output_mutants=False
    )


@pytest.fixture(name="successful_mutation_result")
def fixture_successful_mutation_result():
    """Successful mutation result for A19 THR->ALA."""
    return MutationResult(
        model_id="protprot_complex_1",
        chain="A",
        res_num=19,
        ori_resname="THR",
        target_resname="ALA",
        mutant_scores=(-108.540, -43.234, -275.271, -10.252, 1589.260),
        delta_scores=(5.401, 0.119, 28.482, -0.414, 10.470),
        success=True
    )


@pytest.fixture(name="successful_mutation_result_2")
def fixture_successful_mutation_result_2():
    """Successful mutation result for A20 ILE->ALA."""
    return MutationResult(
        model_id="protprot_complex_1",
        chain="A",
        res_num=20,
        ori_resname="ILE",
        target_resname="ALA",
        mutant_scores=(-18.540, -41.234, -270.271, -11.252, 589.260),
        delta_scores=(4.401, 0.109, 27.482, -0.3, 0.470),
        success=True
    )


@pytest.fixture(name="failed_mutation_result")
def fixture_failed_mutation_result():
    """Failed mutation result."""
    return MutationResult(
        model_id="protprot_complex_1",
        chain="A",
        res_num=20,
        ori_resname="ILE",
        target_resname="ALA",
        mutant_scores=(0, 0, 0, 0, 0),
        delta_scores=(0, 0, 0, 0, 0),
        success=False,
        error_msg="Mutation failed"
    )

# custom mocks 
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


# tests starts here
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


"""Test MutationResult"""
def test_mutation_result_creation():
    """Test MutationResult dataclass creation and attributes."""
    # Basic creation with all required fields
    result = MutationResult(
        model_id="test.pdb",
        chain="A",
        res_num=19,
        ori_resname="THR",
        target_resname="ALA",
        mutant_scores=(-108.5, -43.2, -275.3, -10.3, 1589.3),
        delta_scores=(4.401, 0.109, 27.482, -0.3, 0.470),
        success=True
    )
    # Test all attributes
    assert result.model_id == "test.pdb"
    assert result.chain == "A"
    assert result.res_num == 19
    assert result.ori_resname == "THR"
    assert result.target_resname == "ALA"
    assert result.mutant_scores == (-108.5, -43.2, -275.3, -10.3, 1589.3)
    assert result.delta_scores == (4.401, 0.109, 27.482, -0.3, 0.470)
    assert result.success is True
    assert result.error_msg is None


def test_mutation_result_success_failure_states():
    """Test MutationResult success/failure states."""
    # test successful mutation result
    success_result = MutationResult(
        model_id="success_model",
        chain="A",
        res_num=10,
        ori_resname="VAL",
        target_resname="ALA",
        mutant_scores=(-95.2, -28.1, -290.5, -12.0, 1450.8),
        delta_scores=(8.3, 2.1, 15.7, -0.8, 45.2),
        success=True
    )
    assert success_result.success is True
    assert success_result.error_msg is None  
    assert all(isinstance(score, float) for score in success_result.mutant_scores)
    assert all(isinstance(delta, float) for delta in success_result.delta_scores)
    assert len(success_result.mutant_scores) == 5 
    assert len(success_result.delta_scores) == 5

    # test failed mutation result with error message
    failure_result = MutationResult(
        model_id="failure_model",
        chain="B",
        res_num=25,
        ori_resname="PHE",
        target_resname="ALA",
        mutant_scores=(0, 0, 0, 0, 0),
        delta_scores=(0, 0, 0, 0, 0),
        success=False,
        error_msg="Mutation failed"
    )
    assert failure_result.success is False
    assert failure_result.error_msg == "Mutation failed"


"""Test InterfaceScanner"""
def test_interface_scanner_init(protprot_model_list, params):
    """Test InterfaceScanner initialization in haddock mode."""
    model = protprot_model_list[0]
    scanner = InterfaceScanner(
        model=model,
        mutation_res="GLY",
        params=params,
        library_mode=False
    )
    assert scanner.model_path == model.rel_path
    assert scanner.model_id == model.file_name.removesuffix('.pdb')
    assert scanner.mutation_res == "GLY"
    assert scanner.library_mode is False 
    # extracts chain-specific residue lists from the parameters 
    # and arrange them into a dictionary for easy comparison with user input
    expected_filter = {key[-1]: value for key, value in params.items() if key.startswith("resdic")}
    assert scanner.filter_resdic == expected_filter


def test_interface_scanner_init_library_mode(complex_pdb, params):
    """Test InterfaceScanner initialization in library mode."""
    scanner = InterfaceScanner(
        model=complex_pdb,
        mutation_res="ALA",
        params=params,
        library_mode=True
    )    
    assert scanner.model_path == Path(complex_pdb)
    assert scanner.model_id == Path(complex_pdb).stem
    assert scanner.mutation_res == "ALA"
    assert scanner.library_mode is True
    assert scanner.params == params
    expected_filter = {key[-1]: value for key, value in params.items() if key.startswith("resdic")}
    assert scanner.filter_resdic == expected_filter


def test_interface_scanner_init_default_params():
    """Test InterfaceScanner initialization with default parameters."""
    scanner = InterfaceScanner(
        model=Path("test.pdb"),
        library_mode=False
    )
    # check default values 
    assert scanner.mutation_res == "ALA"
    assert scanner.params == {}
    assert scanner.filter_resdic == {}


def test_interface_scanner_run(mocker, interface_scanner):
    """Test InterfaceScanner.run() (in haddock mode, so returns mutation jobs)."""
    # mock stand-alone functions
    mock_calc_score = mocker.patch("haddock.modules.analysis.alascan.scan.calc_score",
        return_value=(-106.7, -29.6, -316.5, -13.8, 1494.7))
    mock_identify_interface = mocker.patch("haddock.modules.analysis.caprieval.capri.CAPRI.identify_interface",
        return_value={"A": [19, 20]})
    mock_get_atoms = mocker.patch("haddock.libs.libalign.get_atoms")
    mock_load_coords = mocker.patch("haddock.libs.libalign.load_coords",
        return_value=({
            ("A", 19, "CA", "THR"): [0, 0, 0], 
            ("A", 20, "CA", "ILE"): [1, 1, 1],})
    )
    mutation_jobs = interface_scanner.run()
    # should return 2 ModelPointMutation jobs
    assert isinstance(mutation_jobs, list)
    assert len(mutation_jobs) == 2
    assert all(isinstance(job, ModelPointMutation) for job in mutation_jobs)
    # verify jobs have correct attributes
    job = mutation_jobs[0]
    assert job.model_id == interface_scanner.model_id
    assert job.chain in ["A"]
    assert job.target_resname == "ALA"
    assert job.native_scores == (-106.7, -29.6, -316.5, -13.8, 1494.7)


def test_interface_scanner_run_empty_interface(mocker, interface_scanner):
    """Test InterfaceScanner.run() (in haddock mode) with no interface residues."""
    mocker.patch("haddock.modules.analysis.caprieval.capri.CAPRI.identify_interface",
        return_value={}) # no interface 
    mocker.patch("haddock.modules.analysis.alascan.scan.calc_score",
        return_value=(-106.7, -29.6, -316.5, -13.8, 1494.7)
        )
    mocker.patch("haddock.libs.libalign.get_atoms")
    mocker.patch("haddock.libs.libalign.load_coords", return_value=({}))
    mutation_jobs = interface_scanner.run()
    assert mutation_jobs == []
    

def test_interface_scanner_run_library_mode(mocker, interface_scanner_library, successful_mutation_result, successful_mutation_result_2, monkeypatch):
    """Test InterfaceScanner.run() in library mode, so executes mutations."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        # mock dependencies
        mock_calc_score = mocker.patch("haddock.modules.analysis.alascan.scan.calc_score",
            return_value=(-106.7, -29.6, -316.5, -13.8, 1494.7)
        )
        mock_identify_interface = mocker.patch("haddock.modules.analysis.caprieval.capri.CAPRI.identify_interface",
            return_value={"A": [19, 20]}
        )
        mock_get_atoms = mocker.patch("haddock.libs.libalign.get_atoms",
            return_value=[]
        )
        mock_load_coords = mocker.patch("haddock.libs.libalign.load_coords",
            return_value=({
                ("A", 19, "CA", "THR"): [0, 0, 0], 
                ("A", 20, "CA", "ILE"): [1, 1, 1]
            })
        )
        # mock stand-alone functions
        mock_mutation_run = mocker.patch("haddock.modules.analysis.alascan.scan.ModelPointMutation.run",
            side_effect=[successful_mutation_result, successful_mutation_result_2]
        )
        mock_write_scan_out = mocker.patch("haddock.modules.analysis.alascan.scan.write_scan_out")
        result = interface_scanner_library.run()
        assert result is None
        assert mock_mutation_run.call_count == 2
        mock_write_scan_out.assert_called_once()
        call_args = mock_write_scan_out.call_args
        results_passed = call_args[0][0] 
        model_id_passed = call_args[0][1]
        assert model_id_passed == interface_scanner_library.model_id
        assert len(results_passed) == 2
        assert results_passed[0] == successful_mutation_result
        assert results_passed[1] == successful_mutation_result_2
        assert all(result.success for result in results_passed)
        mock_identify_interface.assert_called_once_with(
            interface_scanner_library.model_path, 
            cutoff=interface_scanner_library.params.get("int_cutoff", 5.0))


def test_interface_scanner_run_library_mode_fails(mocker, interface_scanner_library, successful_mutation_result, failed_mutation_result, monkeypatch):
    """Test InterfaceScanner.run() in library mode with some failed mutations."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        mocker.patch("haddock.modules.analysis.alascan.scan.calc_score",
            return_value=(-106.7, -29.6, -316.5, -13.8, 1494.7)
        )
        mocker.patch("haddock.modules.analysis.caprieval.capri.CAPRI.identify_interface",
            return_value={"A": [19, 20]}
        )
        mocker.patch("haddock.libs.libalign.get_atoms", return_value=[])
        mocker.patch("haddock.libs.libalign.load_coords",
            return_value=({
                ("A", 19, "CA", "THR"): [0, 0, 0], 
                ("A", 20, "CA", "ILE"): [1, 1, 1]
            })
        )
        mock_write_scan_out = mocker.patch("haddock.modules.analysis.alascan.scan.write_scan_out")
        # mock ModelPointMutation.run() with mixed results
        mock_mutation_run = mocker.patch("haddock.modules.analysis.alascan.scan.ModelPointMutation.run",
            side_effect=[successful_mutation_result, failed_mutation_result]
        )
        result = interface_scanner_library.run()
        assert result is None
        assert mock_mutation_run.call_count == 2
        mock_write_scan_out.assert_called_once()
        results_passed = mock_write_scan_out.call_args[0][0]        
        assert len(results_passed) == 2
        assert results_passed[0].success is True
        assert results_passed[1].success is False
        assert results_passed[1].error_msg == "Mutation failed"


def test_interface_scanner_run_library_mode_no_mutations(mocker, interface_scanner_library, monkeypatch):
    """Test InterfaceScanner.run() in library mode when no mutations are created."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        mocker.patch("haddock.modules.analysis.alascan.scan.calc_score",
            return_value=(-113.941, -43.353, -303.753, -9.838, 1579.730)
        )
        mocker.patch("haddock.modules.analysis.caprieval.capri.CAPRI.identify_interface",
            return_value={}  # No interface residues
        )
        mocker.patch("haddock.libs.libalign.get_atoms", return_value=[])
        mocker.patch("haddock.libs.libalign.load_coords", return_value=({}, {}))
        mock_write_scan_out = mocker.patch("haddock.modules.analysis.alascan.scan.write_scan_out")
        result = interface_scanner_library.run()
        assert result is None
        mock_write_scan_out.assert_called_once_with([], interface_scanner_library.model_id)


def test_interface_scanner_chain_filtering(mocker, params, complex_pdb):
    """Test chain filtering in InterfaceScanner."""
    params_with_chains = {**params, "chains": ["A"]}
    mock_identify_interface = mocker.patch("haddock.modules.analysis.caprieval.capri.CAPRI.identify_interface",
        return_value={"A": [19, 20], "B": [30, 31]} )
    mock_calc_score = mocker.patch("haddock.modules.analysis.alascan.scan.calc_score",
        return_value=(-106.7, -29.6, -316.5, -13.8, 1494.7))
    mock_get_atoms = mocker.patch("haddock.libs.libalign.get_atoms")
    mock_load_coords = mocker.patch("haddock.libs.libalign.load_coords",
        return_value=({
            ("A", 19, "CA", "THR"): [0, 0, 0],
            ("A", 20, "CA", "ILE"): [1, 1, 1],
            ("B", 30, "CA", "VAL"): [2, 2, 2],
            ("B", 31, "CA", "LEU"): [3, 3, 3]
        })
    )
    scanner = InterfaceScanner(
        model=complex_pdb,
        params=params_with_chains,
        library_mode=False
    )
    mutation_jobs = scanner.run()
    # Should only include chain A mutations
    chains_in_jobs = {job.chain for job in mutation_jobs}
    assert chains_in_jobs == {"A"}
    assert len(mutation_jobs) == 2

def test_interface_scanner_residue_filtering(mocker, complex_pdb, params):
    """Test residue filtering with resdic parameters."""
    params_with_resdic = {**params, "resdic_A": [19], "resdic_B": [30]}
    mock_identify_interface = mocker.patch("haddock.modules.analysis.caprieval.capri.CAPRI.identify_interface",
        return_value={"A": [19, 20, 21], "B": [30, 31, 32]}
    )
    mock_calc_score = mocker.patch("haddock.modules.analysis.alascan.scan.calc_score",
        return_value=(-113.941, -43.353, -303.753, -9.838, 1579.730)
    )
    mock_get_atoms = mocker.patch("haddock.libs.libalign.get_atoms")
    mock_load_coords = mocker.patch("haddock.libs.libalign.load_coords",
        return_value=({
            ("A", 19, "CA", "THR"): [0, 0, 0],
            ("B", 30, "CA", "VAL"): [1, 1, 1]
        }, {})
    )
    scanner = InterfaceScanner(
        model=complex_pdb,
        params=params_with_resdic,
        library_mode=False
    )
    mutation_jobs = scanner.run()
    expected_residues = {(job.chain, job.res_num) for job in mutation_jobs}
    # only residues in resdic_*,  not entire inteface
    assert expected_residues == {("A", 19), ("B", 30)}


def test_interface_scanner_residue_filtering_not_in_interface(mocker, complex_pdb, params):
    """Test residue filtering when user residues are not in interface."""
    params_with_resdic = {**params, "resdic_A": [999], "resdic_B": [999]}  # Not in interface
    mock_identify_interface = mocker.patch("haddock.modules.analysis.caprieval.capri.CAPRI.identify_interface",
        return_value={"A": [19, 20], "B": [30, 31]}  
    )
    mock_calc_score = mocker.patch("haddock.modules.analysis.alascan.scan.calc_score",
        return_value=(-113.941, -43.353, -303.753, -9.838, 1579.730)
    )
    mock_get_atoms = mocker.patch("haddock.libs.libalign.get_atoms")
    mock_load_coords = mocker.patch(
        "haddock.libs.libalign.load_coords",
        return_value=({}, {})
    )
    scanner = InterfaceScanner(
        model=complex_pdb, 
        params=params_with_resdic,
        library_mode=False
    )
    mutation_jobs = scanner.run()
    # no mutatioins to return
    assert len(mutation_jobs) == 0


def test_interface_scanner_skip_same_residue_mutation(mocker, complex_pdb, params):
    """Test that mutations to same residue type are skipped."""
    mock_identify_interface = mocker.patch("haddock.modules.analysis.caprieval.capri.CAPRI.identify_interface",
        return_value={"A": [19, 20]} # here 19 is TRP, so should be not mutated 
    )
    mock_calc_score = mocker.patch("haddock.modules.analysis.alascan.scan.calc_score",
        return_value=(-113.941, -43.353, -303.753, -9.838, 1579.730)
    )
    scanner = InterfaceScanner(
        model=complex_pdb, 
        mutation_res="THR",
        params=params,
        library_mode=False
    )
    mutation_jobs = scanner.run() 
    assert len(mutation_jobs) == 1

"""Test ModelPointMutation"""
def test_model_point_mutation_init(complex_pdb):
    """Test ModelPointMutation initialization."""
    native_scores = (-106.7, -29.6, -316.5, -13.8, 1494.7)
    
    job = ModelPointMutation(
        model_path=complex_pdb,
        model_id="protprot_complex_1",
        chain="A",
        res_num=19,
        ori_resname="THR",
        target_resname="ALA",
        native_scores=native_scores,
        output_mutants=True
    )
    assert job.model_path == Path(complex_pdb)
    assert job.model_id == "protprot_complex_1"
    assert job.chain == "A"
    assert job.res_num == 19
    assert job.ori_resname == "THR"
    assert job.target_resname == "ALA"
    assert job.native_scores == native_scores
    assert job.output_mutants is True


def test_model_point_mutation_run(mocker, mutation_job, monkeypatch):
    """Test successful ModelPointMutation execution."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        mock_mutate = mocker.patch("haddock.modules.analysis.alascan.scan.mutate",
            return_value=Path("mutant.pdb")
        )
        mutant_scores = (-101.3, -28.5, -344.0, -14.2, 1484.2)
        mock_calc_score = mocker.patch("haddock.modules.analysis.alascan.scan.calc_score",
            return_value=mutant_scores
        )        
        # mock file operations
        mocker.patch("os.path.exists", return_value=True)
        mocker.patch("os.remove")
        mocker.patch("shutil.rmtree")
        # run ModelPointMutation.run()
        result = mutation_job.run()
        # test successful result
        assert isinstance(result, MutationResult)
        assert result.success is True
        assert result.model_id == "test_model"
        assert result.chain == "A"
        assert result.res_num == 19
        assert result.ori_resname == "THR"
        assert result.target_resname == "ALA"
        # check resulted mutant_scores aproximately match expected mutant_scores
        for result_score, expected_score in zip(result.mutant_scores, mutant_scores):
            assert pytest.approx(result_score, abs=0.2) == expected_score
        # test cleanup was called
        mock_mutate.assert_called_once_with(mutation_job.model_path, "A", 19, "ALA")
        mock_calc_score.assert_called_once()


def test_model_point_mutation_run_with_output_mutants_false(mocker, mutation_job, monkeypatch):
    """Test ModelPointMutation with output_mutants=False (default)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        mutant_pdb = Path("mutant.pdb")
        em_mutant_pdb = Path("mutant_hs.pdb")
        mutant_pdb.touch()
        em_mutant_pdb.touch()        
        mocker.patch("haddock.modules.analysis.alascan.scan.mutate",
            return_value=mutant_pdb
        )
        mocker.patch(
            "haddock.modules.analysis.alascan.scan.calc_score",
            return_value=(-101.3, -28.5, -344.0, -14.2, 1484.2)
        )
        # Mock file removal operations
        mock_remove = mocker.patch("os.remove")
        mocker.patch("shutil.rmtree")
        # run
        result = mutation_job.run()
        # Test that both files are since output_mutants=False
        assert result.success is True
        assert mock_remove.call_count >= 1


def test_model_point_mutation_run_with_output_mutants_true(mocker, monkeypatch):
    """Test ModelPointMutation with output_mutants=True."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        job = ModelPointMutation(
            model_path=Path("test.pdb"),
            model_id="test_model",
            chain="A",
            res_num=19,
            ori_resname="THR",
            target_resname="ALA",
            native_scores=(-106.7, -29.6, -316.5, -13.8, 1494.7),
            output_mutants=True # this one
        )
        mutant_pdb = Path("mutant.pdb")
        em_mutant_pdb = Path("mutant_hs.pdb")
        mutant_pdb.touch()
        em_mutant_pdb.touch()
        mocker.patch("haddock.modules.analysis.alascan.scan.mutate",
            return_value=mutant_pdb
        )
        mocker.patch("haddock.modules.analysis.alascan.scan.calc_score",
            return_value=(-101.3, -28.5, -344.0, -14.2, 1484.2)
        )
        # mock file operations
        mock_move = mocker.patch("shutil.move")
        mocker.patch("shutil.rmtree")
        # run 
        result = job.run()
        # shutil.move should be called to move em file over original
        if em_mutant_pdb.exists():
            mock_move.assert_called()
        assert result.success is True


def test_model_point_mutation_fail():
    """Test ModelPointMutation with non-existent model file."""
    job = ModelPointMutation(
        model_path=Path("non-pdb-here.pdb"),
        model_id="test",
        chain="A",
        res_num=1,
        ori_resname="ALA", 
        target_resname="TRP",
        native_scores=(0, 0, 0, 0, 0),
        output_mutants=False
    )
    result = job.run()
    # Should return failed result
    assert result.success is False


def test_model_point_mutation_cleanup(mocker, mutation_job, monkeypatch):
    """Test that temporary directories are cleaned up on success."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        mocker.patch("haddock.modules.analysis.alascan.scan.mutate",
            return_value=Path("mutant.pdb")
        )
        mocker.patch("haddock.modules.analysis.alascan.scan.calc_score",
            return_value=(-101.3, -28.5, -344.0, -14.2, 1484.2)
        )
        mocker.patch("os.path.exists", return_value=True)
        mocker.patch("os.remove")
        # mock cleanup calls
        mock_rmtree = mocker.patch("shutil.rmtree")
        # run
        result = mutation_job.run()
        # verify success + cleanup 
        assert result.success is True
        mock_rmtree.assert_called()


"""Test stand-alone functions"""
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


def test_alascan_cluster_analysis(protprot_input_list, scan_file, monkeypatch):
    """Test alascan_cluster_analysis."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        shutil.copy(scan_file, Path("scan_protprot_complex_1.tsv"))
        shutil.copy(scan_file, Path("scan_protprot_complex_2.tsv"))
        alascan_cluster_analysis(protprot_input_list)

        assert Path("scan_clt_unclustered.tsv").exists()

        protprot_input_list[1].clt_id = 1
        alascan_cluster_analysis(protprot_input_list)

        assert Path("scan_clt_1.tsv").exists()
        assert Path("scan_clt_unclustered.tsv").exists()


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


def test_write_scan_out_with_mutation_results(successful_mutation_result, successful_mutation_result_2, monkeypatch):
    """Test write_scan_out."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        results = [successful_mutation_result, successful_mutation_result_2]
        write_scan_out(results, "protprot_complex_1")        
        # test output file creation
        output_file = Path("scan_protprot_complex_1.tsv")
        assert output_file.exists()
        # test file content format
        content = output_file.read_text()
        assert "# `alascan` results for protprot_complex_1" in content
        assert "# native score =" in content
        assert "delta_score" in content
        assert "z_score" in content
        # test data content
        df = pd.read_csv(output_file, sep="\t", comment="#")
        assert len(df) == 2
        assert "z_score" in df.columns
        assert df["chain"].tolist() == ["A", "A"]
        assert df["res"].tolist() == [19, 20]
        assert df["ori_resname"].tolist() == ["THR", "ILE"]
        assert df["end_resname"].tolist() == ["ALA", "ALA"]