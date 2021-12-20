"""Test libcns."""
import os
from pathlib import Path

import pytest

from haddock.libs import libcns


@pytest.mark.parametrize(
    'value',
    [
        {'dic': 1},
        tuple([1, 2]),
        [1, 2],
        set([1, 2]),
        ]
    )
def test_empty_vars_error(value):
    """Test empty vars of types that are not supported."""
    with pytest.raises(TypeError):
        libcns.filter_empty_vars(value)


@pytest.mark.parametrize(
    'value',
    [
        True,
        False,
        1,
        12,
        3453.543,
        'str',
        Path('path'),
        '',
        ]
    )
def test_empty_vars_True(value):
    """Test empty vars of types that are not supported."""
    result = libcns.filter_empty_vars(value)
    assert result is True
    assert type(result) is bool


@pytest.mark.parametrize(
    'value',
    [None],
    )
def test_empty_vars_False(value):
    """Test empty vars of types that are not supported."""
    result = libcns.filter_empty_vars(value)
    assert result is False
    assert type(result) is bool


def test_load_workflow_params():
    """Test workflow params."""
    params = {
        'var1': 1,
        'var2': 'some string',
        'var3': True,
        'var4': Path('some/path'),
        'var5': 5.5,
        'var6': '',
        'var7': None,
        }

    result = libcns.load_workflow_params(params)

    expected = (
        f'{os.linesep}'
        f'! Parameters{os.linesep}'
        f'eval ($var1=1){os.linesep}'
        f'eval ($var2="some string"){os.linesep}'
        f'eval ($var3=true){os.linesep}'
        f'eval ($var4="some/path"){os.linesep}'
        f'eval ($var5=5.5){os.linesep}'
        f'eval ($var6=""){os.linesep}'
        )

    assert result == expected


def test_load_link():
    """Test the loading of the .link file."""
    mol_link = "/some/path"
    result = libcns.load_link(mol_link)

    expected = (
        f'{os.linesep}'
        f'! Link file{os.linesep}'
        f'eval ($link_file = "/some/path" ){os.linesep}'
        )

    assert result == expected


def test_load_trans_vectors():
    pass


def test_load_axis():
    pass


def test_load_waterbox():
    pass


def test_load_ambig():
    pass


def test_load_unambig():
    pass


def test_load_hbond():
    pass


def test_load_tensor_tbl():
    pass


def test_prepare_output():
    pass


def test_load_protonation_state():
    pass


def test_prepare_multiple_input():
    pass


def test_prepare_single_input():
    pass


def test_prepare_cns_input():
    pass


def test_prepare_expected_pdb():
    pass
