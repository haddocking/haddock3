"""Test the rmsdmatrix module."""
import os
from pathlib import Path
import tempfile

import numpy as np
import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.rmsdmatrix import DEFAULT_CONFIG as rmsd_pars
from haddock.modules.analysis.rmsdmatrix import HaddockModule
from haddock.modules.analysis.rmsdmatrix.rmsd import (
    get_pair,
    rmsd_dispatcher,
    XYZWriter,
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
    with pytest.raises(ValueError):
        get_pair(neg_nmodels, 1)
    with pytest.raises(ValueError):
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
    """Test overall rmsdmatrix module."""
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


def test_xyzwriter(input_protdna_models):
    "test XYZWriter"
    with tempfile.TemporaryDirectory() as tmpdir:
        common_keys = [
            ('A', 10, 'N'),
            ('A', 10, 'CA'),
            ('A', 32, 'N'),
            ('B', 38, 'C6'),
        ]
        exp_output = Path(tmpdir, "test.xyz")
        xyzwriter_obj = XYZWriter(
            model_list=input_protdna_models,
            output_name=exp_output,
            core=1,
            n_atoms=4,
            common_keys=common_keys,
            filter_resdic=None,
            allatoms=False,
            )
        xyzwriter_obj.run()
        
        assert exp_output.exists()
        # check the content
        with open(exp_output) as f:
            content = f.read()
        exp_content_list = [
            "4",
            "",
            "A10N 12.163 -4.828 -1.492",
            "A10CA 11.392 -5.83 -0.759",
            "A32N 0.127 -0.923 -0.471",
            "B38C6 -6.564 -18.595 -5.571",
            "4",
            "",
            "A10N -11.179 7.766 -1.6",
            "A10CA -9.966 8.514 -1.292",
            "A32N -1.291 -0.108 -2.675",
            "B38C6 14.422 15.302 -5.743",
        ]

        exp_content = os.linesep.join(exp_content_list) + os.linesep
        assert content == exp_content
