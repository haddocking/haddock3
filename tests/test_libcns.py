"""Test libcns."""
import os
from pathlib import Path

import pytest

from haddock import EmptyPath
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
        EmptyPath(),  # empty paths need to be written as ""
        ]
    )
def test_empty_vars_True(value):
    """Test empty vars of types that are not supported."""
    result = libcns.filter_empty_vars(value)
    assert result is True
    assert type(result) is bool


@pytest.mark.parametrize(
    'value',
    [
        None,
        '',
        float('nan'),
        ],
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

    result = libcns.load_workflow_params(**params)

    expected = (
        f'{os.linesep}'
        f'! Parameters{os.linesep}'
        f'eval ($var1=1){os.linesep}'
        f'eval ($var2="some string"){os.linesep}'
        f'eval ($var3=true){os.linesep}'
        f'eval ($var4="some/path"){os.linesep}'
        f'eval ($var5=5.5){os.linesep}'
        )

    assert result == expected
