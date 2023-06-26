"""Test haddock3-traceback client."""

import os
import shutil
from pathlib import Path

import pandas as pd
import pytest

from haddock.clis.cli_traceback import main

from . import golden_data


@pytest.fixture
def rigid_json():
    """Provide example rigidbody io.json file."""
    return Path(golden_data, "io_rigid.json")


@pytest.fixture
def flexref_json():
    """Provide example flexref io.json file."""
    return Path(golden_data, "io_flexref.json")


def test_main(rigid_json, flexref_json):
    """Test haddock3-traceback client."""
    # build fake run_dir
    run_dir = "example_dir"
    step_dirs = [os.path.join(run_dir, "1_rigidbody"),
                 os.path.join(run_dir, "4_flexref")]
    
    if os.path.isdir(run_dir):
        shutil.rmtree(run_dir)
    os.mkdir(run_dir)
    os.mkdir(step_dirs[0])
    os.mkdir(step_dirs[1])
    shutil.copy(rigid_json, os.path.join(step_dirs[0], "io.json"))
    shutil.copy(flexref_json, os.path.join(step_dirs[1], "io.json"))

    # run haddock3-traceback
    main(run_dir)

    # check traceback folder exists
    assert os.path.isdir(os.path.join(run_dir, "traceback"))

    # check traceback files exist
    tr_file = os.path.join(run_dir, "traceback", "traceback.tsv")
    assert os.path.isfile(tr_file)

    obs_tr = pd.read_csv(tr_file, sep="\t", dtype=str)
    exp_tr = [["00_topo1", "00_topo2", "1_rigidbody", "1_rigidbody_rank", "4_flexref", "4_flexref_rank"],  # noqa: E501
              ["4G6K_fv_haddock.psf", "4I1B-matched_haddock.psf", "rigidbody_2.pdb", "1", "flexref_1.pdb", "1"],  # noqa: E501
              ["4G6K_fv_haddock.psf", "4I1B-matched_haddock.psf", "rigidbody_4.pdb", "2", "flexref_2.pdb", "2"],  # noqa: E501
              ["4G6K_fv_haddock.psf", "4I1B-matched_haddock.psf", "rigidbody_1.pdb", "3", "-", "-"],  # noqa: E501
              ["4G6K_fv_haddock.psf", "4I1B-matched_haddock.psf", "rigidbody_3.pdb", "4", "-", "-"]]  # noqa: E501
    exp_tr_df = pd.DataFrame(exp_tr[1:], columns=exp_tr[0])

    assert obs_tr.columns.tolist() == exp_tr_df.columns.tolist()
    assert obs_tr.equals(exp_tr_df)

    # clean up
    shutil.rmtree(run_dir)
