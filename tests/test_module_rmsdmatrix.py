"""Test the rmsdmatrix module."""
import os
from pathlib import Path

import numpy as np
import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.rmsdmatrix import DEFAULT_CONFIG as rmsd_pars
from haddock.modules.analysis.rmsdmatrix import HaddockModule
from haddock.modules.analysis.rmsdmatrix.rmsd import (
    RMSD,
    RMSDJob,
    get_pair,
    rmsd_dispatcher,
    )

from . import golden_data


@pytest.fixture
def input_protdna_models():
    """Prot-DNA models using for emscoring output."""
    return [
        PDBFile(
            Path(golden_data, "protdna_complex_1.pdb"),
            path=golden_data,
            score=42.0
            ),
        PDBFile(
            Path(golden_data, "protdna_complex_2.pdb"),
            path=golden_data,
            score=28.0
            )]


def test_get_pair():
    """Test the get_pair() function."""
    nmodels_vec = [10, 2, 200]
    idx_vec = [10, 0, 396]
    
    expexted_ivec = [1, 0, 1]
    expected_jvec = [3, 1, 199]
    for n in range(len(nmodels_vec)):
        observed_i, observed_j = get_pair(nmodels_vec[n], idx_vec[n])

        assert (observed_i, observed_j) == (expexted_ivec[n], expected_jvec[n])


def test_get_pair_error():
    """Test negative arguments to the get_pair() function."""
    neg_nmodels = -1
    neg_idx = -1
    with pytest.raises(Exception):
        get_pair(neg_nmodels, 1)
    with pytest.raises(Exception):
        get_pair(1, neg_idx)


def test_rmsd_dispatcher():
    """Test optimal dispatching."""
    nmodels_vec = [1, 10, 5]
    tot_npairs_vec = [1, 45, 10]
    ncores_vec = [1, 2, 10]
    
    expected_pairs_vec = [[1], [23, 22], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]
    expected_ref_vec = [[0], [0, 2], [0, 0, 0, 0, 1, 1, 1, 2, 2, 3]]
    expected_mod_vec = [[1], [1, 9], [1, 2, 3, 4, 2, 3, 4, 3, 4, 4]]

    for n in range(len(nmodels_vec)):
        dispatcher_output = rmsd_dispatcher(
            nmodels_vec[n],
            tot_npairs_vec[n],
            ncores_vec[n]
            )
        
        assert dispatcher_output[0] == expected_pairs_vec[n]

        assert dispatcher_output[1] == expected_ref_vec[n]

        assert dispatcher_output[2] == expected_mod_vec[n]


def test_overall_rmsd(input_protdna_models):
    """"""
    rmsd_module = HaddockModule(
        order=2,
        path=Path("2_rmsdmatrix"),
        initial_params=rmsd_pars
        )
    rmsd_module.previous_io.output = input_protdna_models
    rmsd_module._run()

    ls = os.listdir()

    assert "rmsd.matrix" in ls

    assert "rmsd_matrix.json" in ls

    # check correct rmsd matrix
    rmsd_matrix = open("rmsd.matrix").read()
    
    expected_rmsd_matrix = "1 2 2.257" + os.linesep

    assert rmsd_matrix == expected_rmsd_matrix

    os.unlink(Path("rmsd.matrix"))
    os.unlink(Path("rmsd_matrix.json"))
    os.unlink(Path("io.json"))


def test_RMSD_class(input_protdna_models):
    """Test focusing on the RMSD class."""
    params = {}
    rmsd_obj = RMSD(
        input_protdna_models,
        core=0,
        npairs=1,
        start_ref=0,
        start_mod=1,
        output_name="rmsd_0.matrix",
        path=Path("."),
        params=params
        )
    rmsd_obj.run()

    expected_data = np.array([[1, 2, 2.257]])

    np.testing.assert_allclose(rmsd_obj.data, expected_data, atol=0.001)


def test_RMSD_filter_resdic(input_protdna_models):
    """Test filter_resdic."""
    params = {"resdic_A": [1, 2, 3], "resdic_B": [4, 5, 6]}

    # this test only considers the initialisation of the object
    rmsd_obj = RMSD(
        input_protdna_models,
        core=0,
        npairs=1,
        start_ref=0,
        start_mod=1,
        output_name="rmsd_0.matrix",
        path=Path("."),
        params=params
        )

    expected_resdic = {"A": [1, 2, 3], "B": [4, 5, 6]}

    assert rmsd_obj.filter_resdic == expected_resdic


def test_RMSDJob(input_protdna_models):
    """Test the RMSD job."""
    rmsd_obj = RMSD(
        input_protdna_models,
        core=0,
        npairs=1,
        start_ref=0,
        start_mod=1,
        output_name="rmsd_0.matrix",
        path=Path("."),
        )

    job_f = "fake_rmsd.job"
    params = {}

    job = RMSDJob(
        job_f,
        params,
        rmsd_obj
        )
    
    assert job.rmsd_obj == rmsd_obj

    assert job.output == job_f
