"""Test tables gear."""
from math import isclose, isnan
from pathlib import Path

import pytest

from haddock.gear.tables import (
    TableFormatError,
    convert_row_to_column_table,
    convert_ssc_csv,
    convert_ssc_tsv,
    create_adjusted_col_width_table,
    create_table,
    format_value,
    parse_table_to_data_dict,
    read_table_to_data_dict,
    value_to_str,
    )
from haddock.libs.libontology import PDBFile

from . import adjusted_col_width_table


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
    result = create_adjusted_col_width_table(d)
    assert result == table1


def test_create_table_2():
    d = {}
    d["col1"] = list(range(4))
    d["col2"] = ["haddock", "is", "amazing", None]
    d["col3"] = [1.001, 123.543, 654.1, 98.999999]
    result = create_adjusted_col_width_table(d, header="# some header")
    assert result == table2


def test_create_table_error():
    d = {}
    d["col1"] = list(range(5))
    d["col2"] = ["haddock", "is", "amazing", None]
    d["col3"] = [1.001, 123.543, 654.1, 98.999999]
    with pytest.raises(TableFormatError):
        create_adjusted_col_width_table(d)


def test_create_table_error_2():
    d = {}
    d["col1"] = list(range(4))
    d["col2"] = ["haddock", "is", "amazing"]
    d["col3"] = [1.001, 123.543, 654.1, 98.999999]
    with pytest.raises(TableFormatError):
        create_adjusted_col_width_table(d)


@pytest.mark.parametrize('table', [table1, table2])
def test_parse_table(table):
    d = {}
    d["col1"] = list(range(4))
    d["col2"] = ["haddock", "is", "amazing", None]
    d["col3"] = [1.001, 123.543, 654.1, 99.000]
    result = parse_table_to_data_dict(table)
    assert result.keys() == d.keys()
    assert result["col1"] == d["col1"]
    assert result["col2"] == d["col2"]
    assert all(isinstance(v, float) for v in result["col3"])
    assert all(str(i) == str(j) for i, j in zip(result["col1"], d["col1"]))


def test_read_table():
    d = {}
    d["col_name_1"] = ["value1", float('nan'), "value7"]
    d["col_name_2"] = [None, "value4", None]
    d["col_name_3"] = ["value5", 6, 1.2]
    result = read_table_to_data_dict(adjusted_col_width_table)
    assert result.keys() == d.keys()
    assert result["col_name_1"][0] == "value1"
    assert isnan(result["col_name_1"][1])
    assert result["col_name_1"][2] == "value7"
    assert result["col_name_2"][0] is None
    assert result["col_name_2"][1] == "value4"
    assert result["col_name_2"][2] is None
    assert result["col_name_3"][0] == "value5"
    assert result["col_name_3"][1] == 6
    assert isclose(result["col_name_3"][2], 1.2)


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


tsv = """# some comment
col_name_1\tcol_name_2\tcol_name_3
value1\tnone\tvalue5
nan\tvalue4\t6
value7\t-\t1.2"""

csv = """# some comment
col_name_1,col_name_2,col_name_3
value1,none,value5
nan,value4,6
value7,-,1.2"""


def test_convert_ssc_tsv():
    result = convert_ssc_tsv(adjusted_col_width_table)
    assert result == tsv


def test_convert_ssc_csv():
    result = convert_ssc_csv(adjusted_col_width_table)
    assert result == csv


def test_create_table_csv():
    d = {}
    d["col_name_1"] = ["value1", "none", "value5"]
    d["col_name_2"] = [float('nan'), "value4", 6]
    d["col_name_3"] = ["value7", None, 1.2]
    result = create_table(
        d,
        sep=",",
        header="# some comment",
        float_fmt="{:.1f}",
        )
    assert result == csv
