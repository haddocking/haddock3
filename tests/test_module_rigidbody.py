"""Test the rigidbody module."""
from pathlib import Path

from haddock.modules.sampling.rigidbody import DEFAULT_CONFIG as rigidbody_pars
from haddock.modules.sampling.rigidbody import HaddockModule


def test_prev_fnames():
    """Tests the correct retrieval of ambiguous restraints information."""
    rbody_module = HaddockModule(
        order=1,
        path=Path("1_rigidbody"),
        initial_params=rigidbody_pars,
        )
    prev_ambig_fnames = [None for md in range(rbody_module.params["sampling"])]
    diff_ambig_fnames = rbody_module.get_ambig_fnames(prev_ambig_fnames)
    assert diff_ambig_fnames is None
