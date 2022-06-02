"""Test pretty table module."""
# noqa: W291
from pathlib import Path

import pytest

from haddock.gear.pretty_table import (
    TableFormatError,
    convert_row_to_column_table,
    create_human_readable_table,
    format_value,
    value_to_str,
    )
from haddock.libs.libontology import PDBFile


table1 = """col1      col2       col3 
 0      haddock      1.001
 1         is      123.543
 2      amazing    654.100
 3         -        99.000

"""

table2 = """# some header
col1      col2       col3 
 0      haddock      1.001
 1         is      123.543
 2      amazing    654.100
 3         -        99.000

"""


def test_create_table_1():
    d = {}
    d["col1"] = list(range(4))
    d["col2"] = ["haddock", "is", "amazing", None]
    d["col3"] = [1.001, 123.543, 654.1, 98.999999]
    result = create_human_readable_table(d)
    assert result == table1


def test_create_table_2():
    d = {}
    d["col1"] = list(range(4))
    d["col2"] = ["haddock", "is", "amazing", None]
    d["col3"] = [1.001, 123.543, 654.1, 98.999999]
    result = create_human_readable_table(d, header="# some header")
    assert result == table2


def test_create_table_error():
    d = {}
    d["col1"] = list(range(5))
    d["col2"] = ["haddock", "is", "amazing", None]
    d["col3"] = [1.001, 123.543, 654.1, 98.999999]
    with pytest.raises(TableFormatError):
        create_human_readable_table(d)


def test_create_table_error_2():
    d = {}
    d["col1"] = list(range(4))
    d["col2"] = ["haddock", "is", "amazing"]
    d["col3"] = [1.001, 123.543, 654.1, 98.999999]
    with pytest.raises(TableFormatError):
        create_human_readable_table(d)


@pytest.mark.parametrize(
    "in_, expected",
    [
        (PDBFile("some.pdb"), str(Path("..", "haddock3", "some.pdb"))),
        (Path("some", "path", "file.txt"), str(Path("some", "path", "file.txt"))),  # noqa: E501
        ("some string", "some string"),
        (123.1234, "123.123"),
        (10.99999, "11.000"),
        ]
    )
def test_value_PDFFile_to_str(in_, expected):
    result = value_to_str(in_)
    assert result == expected


@pytest.mark.parametrize(
    "in_, expected, spacing",
    [
        (PDBFile("some.pdb"), str(Path("..", "haddock3", "some.pdb")) + "  ", 22),  # noqa: E501
        (Path("some", "path", "file.txt"), str(Path("some", "path", "file.txt")) + "    ", 22),  # noqa: E501
        ("some string", "  some string  ", 15),
        (123.1234, "  123.123", 9),
        (10.99999, "   11.000", 9),
        ]
    )
def test_format_value_str(in_, expected, spacing):
    result = format_value(in_, spacing, " ")
    assert result == expected


def test_convert_row_table():
    d = {}
    d["col1"] = list(range(4))
    d["col2"] = ["haddock", "is", "amazing", None]

    r = {}
    r[0] = {
        "col1": 0,
        "col2": "haddock",
        }
    r[1] = {
        "col1": 1,
        "col2": "is",
        }
    r[2] = {
        "col1": 2,
        "col2": "amazing",
        }
    r[3] = {
        "col1": 3,
        "col2": None,
        }

    result = convert_row_to_column_table(r)
    assert result == d
