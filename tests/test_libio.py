"""Test libio."""
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libio import (
    clean_suffix,
    dot_suffix,
    file_exists,
    folder_exists,
    read_from_yaml,
    write_dic_to_file,
    write_nested_dic_to_file,
    )

from . import emptycfg, haddock3_yaml_cfg_examples


@pytest.mark.parametrize(
    "cfg",
    [
        emptycfg,
        haddock3_yaml_cfg_examples,
        ],
    )
def test_read_from_yaml(cfg):
    """Test read from yaml file."""
    result = read_from_yaml(cfg)
    assert isinstance(result, dict)


def test_write_nested_dic_to_file():
    """Test write nested dictionary to file."""
    f = tempfile.NamedTemporaryFile(delete=False)
    write_nested_dic_to_file(
        data_dict={1: {"something": "something"}},
        output_fname=f.name)

    assert Path(f.name).exists()
    assert Path(f.name).stat().st_size != 0

    Path(f.name).unlink()


def test_write_dic_to_file():
    """Test write dictionary to file."""
    f = tempfile.NamedTemporaryFile(delete=False)
    write_dic_to_file(
        data_dict={"something": "something"},
        output_fname=f.name)

    assert Path(f.name).exists()
    assert Path(f.name).stat().st_size != 0

    Path(f.name).unlink()


@pytest.mark.parametrize(
    "in_,expected",
    [
        (".ext", ".ext"),
        ("ext", ".ext"),
        (".out.gz", ".out.gz"),
        ("out.gz", ".out.gz"),
        ]
    )
def test_dot_suffix(in_, expected):
    result = dot_suffix(in_)
    assert result == expected


@pytest.mark.parametrize(
    "in_,expected",
    [
        (".ext", "ext"),
        ("ext", "ext"),
        (".out.gz", "out.gz"),
        ("out.gz", "out.gz"),
        ]
    )
def test_clean_suffix(in_, expected):
    result = clean_suffix(in_)
    assert result == expected


@pytest.mark.parametrize(
    'i,expected',
    [
        (Path(__file__), Path(__file__)),
        (str(Path(__file__)), Path(__file__)),
        ],
    )
def test_file_exists(i, expected):
    """."""
    r = file_exists(i)
    assert r == expected


@pytest.mark.parametrize(
    'i',
    [
        'some_bad_path',
        Path(__file__).parent,  # this is a folder
        ],
    )
def test_file_exists_wrong(i):
    """."""
    with pytest.raises(ValueError):
        file_exists(i)


def test_folder_exists():
    """."""
    r = folder_exists(Path(__file__).parent)
    assert r == Path(__file__).parent


@pytest.mark.parametrize(
    'i',
    [
        'some_bad_path',
        Path(__file__),  # this is a file
        str(Path(__file__)),  # this is a file
        ],
    )
def test_folder_exists_wrong(i):
    """."""
    with pytest.raises(ValueError):
        folder_exists(i)


def test_folder_exists_wrong_othererror():
    with pytest.raises(TypeError):
        folder_exists("some_bad_path", exception=TypeError)
