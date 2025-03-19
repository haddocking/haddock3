"""Test libcns."""

import os
import random
import tempfile
from pathlib import Path

import pytest

from haddock import EmptyPath
from haddock.libs import libcns
from haddock.libs.libcns import (
    prepare_cns_input,
    prepare_expected_pdb,
    prepare_multiple_input,
)
from haddock.libs.libontology import Format, PDBFile, Persistent


@pytest.mark.parametrize(
    "value",
    [
        {"dic": 1},
        tuple([1, 2]),
        [1, 2],
        set([1, 2]),
    ],
)
def test_empty_vars_error(value):
    """Test empty vars of types that are not supported."""
    with pytest.raises(TypeError):
        libcns.filter_empty_vars(value)


@pytest.mark.parametrize(
    "value",
    [
        True,
        False,
        1,
        12,
        3453.543,
        "str",
        Path("path"),
        EmptyPath(),  # empty paths need to be written as ""
    ],
)
def test_empty_vars_True(value):
    """Test empty vars of types that are not supported."""
    result = libcns.filter_empty_vars(value)
    assert result is True
    assert type(result) is bool


@pytest.mark.parametrize(
    "value",
    [
        None,
        "",
        float("nan"),
    ],
)
def test_empty_vars_False(value):
    """Test empty vars of types that are not supported."""
    result = libcns.filter_empty_vars(value)
    assert result is False
    assert type(result) is bool


@pytest.fixture
def pdbfile():
    """Create a temporary file."""
    with tempfile.NamedTemporaryFile() as pdb_f, tempfile.NamedTemporaryFile() as top_f:

        pdb = PDBFile(
            file_name=pdb_f.name,
            path=Path(pdb_f.name).parent,
            topology=Persistent(
                file_name=top_f.name,
                path=Path(top_f.name).parent,
                file_type=Format.TOPOLOGY,
            ),
        )
        pdb.seed = 42  # type: ignore
        yield pdb


def test_load_workflow_params():
    """Test workflow params."""
    params = {
        "var1": 1,
        "var2": "some string",
        "var3": True,
        "var4": Path("some/path"),
        "var5": 5.5,
        "var6": "",
        "var7": None,
    }

    result = libcns.load_workflow_params(**params)

    expected = (
        f"{os.linesep}"
        f"! Parameters{os.linesep}"
        f"eval ($var1=1){os.linesep}"
        f'eval ($var2="some string"){os.linesep}'
        f"eval ($var3=true){os.linesep}"
        f'eval ($var4="some/path"){os.linesep}'
        f"eval ($var5=5.5){os.linesep}"
    )

    assert result == expected


def test_prepare_cns_input(pdbfile):

    # TODO: Improve this test to increase coverage AND refactor `prepare_cns_input`

    model_number = 1
    observed_cns_input = prepare_cns_input(
        model_number=model_number,
        input_element=pdbfile,
        step_path=Path("."),
        recipe_str="",
        defaults={},
        identifier="",
        ambig_fname="",
        native_segid=False,
        default_params_path=None,
        debug=False,
        seed=pdbfile.seed,
    )

    expected_cns_input = f"""
! Parameters
eval ($ambig_fname="")

! Input structure
structure
  @@{pdbfile.topology.rel_path}
end
coor @@{pdbfile.rel_path}
eval ($input_pdb_filename_1="{pdbfile.rel_path}")
eval ($ncomponents=0)
eval ($seed={pdbfile.seed})

! Output structure
eval ($output_pdb_filename="_{model_number}.pdb")
eval ($count=1)
"""

    assert observed_cns_input == expected_cns_input


def test_prepare_multiple_input(mocker):

    mocker.patch("haddock.libs.libpdb.identify_chainseg", return_value="A")

    observed_input = prepare_multiple_input(
        pdb_input_list=["pdb_1", "pdb_2"],
        psf_input_list=["psf_1"],
    )

    expected_input = """
! Input structure
structure
  @@psf_1
end
coor @@pdb_1
eval ($input_pdb_filename_1="pdb_1")
coor @@pdb_2
eval ($input_pdb_filename_2="pdb_2")
eval ($ncomponents=1)
"""

    assert observed_input == expected_input


def test_prepare_expected_pdb(pdbfile):

    model_number = random.randint(1, 100)
    output_pdb = prepare_expected_pdb(
        model_obj=pdbfile,
        model_nb=model_number,
        path=Path.cwd(),
        identifier="",
    )

    assert isinstance(output_pdb, PDBFile)
    assert output_pdb.file_name == f"_{model_number}.pdb"
    assert output_pdb.topology is not None
    assert output_pdb.topology.file_name == pdbfile.topology.file_name
