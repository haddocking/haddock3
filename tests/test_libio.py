"""Test libio."""
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libio import (
    parse_suffix,
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
        (".ext", "ext"),
        ("ext", "ext"),
        ]
    )
def test_parse_suffix(in_, expected):
    result = parse_suffix(in_)
    assert result == expected
