"""Test libutil."""
from os import cpu_count
from pathlib import Path

import pytest

from haddock import EmptyPath
from haddock.libs.libutil import (
    extract_keys_recursive,
    get_number_from_path_stem,
    non_negative_int,
    parse_ncores,
    recursive_convert_paths_to_strings,
    recursive_dict_update,
    sort_numbered_paths,
    transform_to_list,
    )


@pytest.mark.parametrize(
    'i,expected',
    [
        (1, 1),
        ('1', 1),
        (0, 0),
        ('0', 0),
        ],
    )
def test_non_negative_int(i, expected):
    """."""
    r = non_negative_int(i)
    assert r == expected


@pytest.mark.parametrize('i', [-1, -12, -14])
def test_non_negative_int_error(i):
    """."""
    with pytest.raises(ValueError):
        non_negative_int(i)


@pytest.mark.parametrize(
    'in1,expected',
    [
        ('pdb_1.pdb', 1),
        ('pdb2.pdb', 2),
        ('pdb2.pdb', 2),
        ('pdb_3.pdb', 3),
        ('pdb_1231.pdb', 1231),
        ('pdb_0011.pdb', 11),
        ('pdb_1_.pdb', 1),
        ('pdb_1', 1),
        ('5', 5),
        ('pdb_20200101_1.pdb', 1),
        ],
    )
def test_get_number(in1, expected):
    """Test get number from path."""
    result = get_number_from_path_stem(in1)
    assert result == expected


@pytest.mark.parametrize(
    'in1,expected',
    [
        (
            ['f_1.pdb', 'f_11.pdb', 'f_2.pdb'],
            ['f_1.pdb', 'f_2.pdb', 'f_11.pdb'],
            ),
        (
            ['b.pdb', 'c.pdb', 'a.pdb'],
            ['a.pdb', 'b.pdb', 'c.pdb']),
        ],
    )
def test_sort_numbered_input_1(in1, expected):
    """Test sort numbered inputs."""
    result = sort_numbered_paths(*in1)
    assert result == expected


@pytest.mark.parametrize(
    'in1,error',
    [
        (['f_1.pdb', 'f_11.pdb', 'f_2.pdb'], TypeError),
        ]
    )
def test_sort_numbered_inputs_error(in1, error):
    """Test sort numbered inputs raised Errors."""
    with pytest.raises(error):
        sort_numbered_paths(in1)


def test_recursive_dict_update():
    """Test recursive dict update."""
    a = {"a": 1, "b": {"c": 2, "d": {"e": 3}}}
    _list = list(range(10))
    b = {"a": 2, "b": {"d": {"e": 4}}, "z": {"z1": _list}}
    c = recursive_dict_update(a, b)
    assert a is not c
    assert a["b"] is not c["b"]
    assert a["b"]["d"] is not c["b"]["d"]
    assert b["z"]["z1"] is not c["z"]["z1"]
    assert c == {"a": 2, "b": {"c": 2, "d": {"e": 4}}, "z": {"z1": _list}}


def test_recursive_dict_update_empty():
    """Test recursive dict update."""
    a = {"a": 1, "b": {"c": 2, "d": {"e": 3}}}
    c = recursive_dict_update(a, {})
    assert a is not c
    assert a == c


a = {
    "param1": 1,
    "param2": {"param3": 4, "param4": {"param5": {"param6": 7, "param7": 8}}},
    }
b = {"param1", "param3", "param6", "param7"}


@pytest.mark.parametrize(
    "inp,expected",
    [
        (a, b),
        ]
    )
def test_extract_keys_recursive(inp, expected):
    """Test extract keys recursive."""
    result = set(extract_keys_recursive(inp))
    assert result == expected


@pytest.mark.parametrize(
    "value,expected",
    [
        [1, [1]],
        [Path("a"), [Path("a")]],
        [list(range(10)), list(range(10))],
        [tuple(range(10)), tuple(range(10))],
        [1.1, [1.1]],
        ["a", ["a"]],
        [set([1, 2, 3]), [1, 2, 3]],
        [{"a": 1}, ["a"]],
        [None, [None]],
        ]
    )
def test_transform_to_list(value, expected):
    result = transform_to_list(value)
    assert result == expected


def test_convert_paths_to_strings_recursive():
    """Test converts paths to strings."""
    i = {
        "a": 1,
        "p1": Path("file"),
        "v": [Path("file1"), Path("file2")],
        "v2": {
            "v3": EmptyPath(),
            "v4": Path("file2"),
            }
        }

    e = {
        "a": 1,
        "p1": "file",
        "v": ["file1", "file2"],
        "v2": {
            "v3": "",
            "v4": "file2",
            }
        }

    r = recursive_convert_paths_to_strings(i)
    assert r == e


@pytest.mark.parametrize(
    'n,njobs,maxcpus,expected',
    [
        (10, 10, 10, 10),
        (10, 10, 5, 5),
        (1000, 1000, None, cpu_count() - 1),
        (1000, 1000, False, cpu_count() - 1),
        (1000, 1000, True, cpu_count()),
        (1000, 1, True, 1),
        (1000, 1, False, 1),
        (5, 10, False, min(5, cpu_count() - 1)),
        (None, None, None, cpu_count() - 1),
        ]
    )
def test_parse_ncores(n, njobs, maxcpus, expected):
    """Test parse_ncores function."""
    r = parse_ncores(n, njobs, maxcpus)
    assert r == expected


@pytest.mark.parametrize(
    'maxcpus',
    [
        1.1,
        {},
        (1, 2),
        set([1, 3]),
        '1',
        ]
    )
def test_parse_ncores_error(maxcpus):
    """Test parse_ncores function."""
    with pytest.raises(TypeError):
        parse_ncores(max_cpus=maxcpus)
