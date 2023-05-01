"""Test the flexref module."""
from pathlib import Path

import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.base_cns_module import RestraintError
from haddock.modules.refinement.flexref import DEFAULT_CONFIG as flexref_pars
from haddock.modules.refinement.flexref import HaddockModule

from . import golden_data


@pytest.fixture
def input_models():
    """Prot-DNA models using for flexref input."""
    return [
        PDBFile(
            Path(golden_data, "protdna_complex_1.pdb"),
            path=golden_data,
            score=42.0,
            restr_fname=Path(golden_data, "example_ambig_1.tbl")
            ),
        PDBFile(
            Path(golden_data, "protdna_complex_2.pdb"),
            path=golden_data,
            score=28.0,
            restr_fname=Path(golden_data, "example_ambig_2.tbl")
            )]


def test_prev_fnames(input_models):
    """Tests the correct retrieval of ambiguous restraints information."""
    flexref_module = HaddockModule(
        order=1,
        path=Path("1_rigidbody"),
        initial_params=flexref_pars,
        )
    prev_ambig_fnames = [model.restr_fname for model in input_models]
    obs_ambig_fnames = flexref_module.get_ambig_fnames(prev_ambig_fnames)
    assert obs_ambig_fnames is None
    # now setting previous to true
    flexref_module.params["previous_ambig"] = True
    obs_ambig_fnames = flexref_module.get_ambig_fnames(prev_ambig_fnames)
    assert obs_ambig_fnames == prev_ambig_fnames
    # removing restraints information from a model while previous is true
    #  should cause an exception
    input_models[0].restr_fname = None
    prev_ambig_fnames = [model.restr_fname for model in input_models]
    with pytest.raises(RestraintError):
        obs_ambig_fnames = flexref_module.get_ambig_fnames(prev_ambig_fnames)
