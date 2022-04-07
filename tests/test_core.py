"""Test the core module."""

from haddock.core.defaults import max_molecules_allowed


def test_max_molecules_allowed():
    assert max_molecules_allowed == 20
