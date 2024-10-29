"""Test the flexref module."""

import os
import tempfile
from pathlib import Path

import pytest

from haddock.modules.refinement.flexref import \
    DEFAULT_CONFIG as DEFAULT_FLEXREF_PARAMS
from haddock.modules.refinement.flexref import HaddockModule as Flexref

from . import golden_data


@pytest.fixture(name="flexref")
def fixture_fixture_flexref():
    """???"""
    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)
        yield Flexref(
            order=1,
            path=Path("."),
            initial_params=DEFAULT_FLEXREF_PARAMS,
        )


def test_prev_fnames(flexref, protdna_input_list):
    """Tests the correct retrieval of ambiguous restraints information."""

    protdna_input_list[0].restr_fname = Path(golden_data, "example_ambig_1.tbl")
    protdna_input_list[1].restr_fname = Path(golden_data, "example_ambig_2.tbl")

    prev_ambig_fnames = [model.restr_fname for model in protdna_input_list]
    obs_ambig_fnames = flexref.get_ambig_fnames(prev_ambig_fnames)
    assert obs_ambig_fnames is None
    # now setting previous to true
    flexref.params["previous_ambig"] = True
    obs_ambig_fnames = flexref.get_ambig_fnames(prev_ambig_fnames)
    assert obs_ambig_fnames == prev_ambig_fnames
    # removing restraints information from a model while previous is true
    #  should cause an exception
    protdna_input_list[0].restr_fname = None
    prev_ambig_fnames = [model.restr_fname for model in protdna_input_list]
    # FIXME: this should be a more specific exception
    with pytest.raises(Exception):  # noqa: B017
        obs_ambig_fnames = flexref.get_ambig_fnames(prev_ambig_fnames)
