"""Test the flexref module."""

import tempfile
import shutil
from pathlib import Path

import pytest

from haddock.libs.libio import extract_files_flat
from haddock.modules.refinement.flexref import \
    DEFAULT_CONFIG as DEFAULT_FLEXREF_PARAMS
from haddock.modules.refinement.flexref import HaddockModule as Flexref

from . import golden_data


@pytest.fixture(name="flexref")
def fixture_fixture_flexref(monkeypatch):
    """???"""
    with tempfile.TemporaryDirectory() as tempdir:
        monkeypatch.chdir(tempdir)
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


def test_multiple_ambigs(flexref):
    """Test usage of multiple ambigs."""
    # Copy ambig.tbl.tgz
    reference_archive = Path(golden_data, "ambig.tbl.tgz")
    shutil.copy(reference_archive, flexref.path)
    archive = Path(flexref.path, reference_archive.name)
    # First extract the file
    extract_files_flat(archive, flexref.path)
    # Set the parameter
    flexref.params["ambig_fname"] = archive
    ambig_fnames = flexref.get_ambig_fnames([])
    assert len(ambig_fnames) == 1
    # Remove file from path
    Path('ambig_1.tbl').unlink()
    with pytest.raises(Exception) as _e:
        ambig_fnames = flexref.get_ambig_fnames([])
        assert ambig_fnames is None
