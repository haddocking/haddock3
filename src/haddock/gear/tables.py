"""
Tools to create visually appealing tables.

This ``gear`` provide several tools to create pure text tables, from various
python data structures.
"""
import os
from copy import copy
from pathlib import Path

from haddock.core.exceptions import HaddockError
from haddock.libs.libontology import PDBFile
from haddock.libs.libutil import transform_to_list
from haddock.libs.libfunc import true, is_str_int, is_str_float, nan, none


class TableFormatError(HaddockError):
    """Error related to table format."""

    pass


def create_adjusted_col_width_table(
        data,
        header=None,
        column_spacing="    ",
        **fmt_args,
        ):
    """
    Create a space-separated text table with column width ajusted to values.

    Given a table in dictionary format:

    .. code:: python

            {
                "col_name_1": [value1, value3],
                "col_name_2": [value2, value4],
                "col_name_3": [value5, value6],
                }

    Creates the following table in text::

        col_name_1    col_name_2    col_name_3
          value1        value2        value5
          value3        value4        valut6

    The returned table can be saved to a file, for example:

    .. code:: python

        from pathlib import Path

        table = create_adjusted_col_width_table(table_data)
        Path('table.txt').write_text(table)

    Developers should decide on the suffix of the table file, ``txt`` in the
    above example.  We suggest ``txt``.

    Examples, on how to parse theses tables:

    These tables are space-separated values, where the number of spaces
    varies to make proper alignment of the text. If you are using Python,
    you can parse these tables with, for example:

    .. code:: python

        import os
        with open('table', 'r') as fin:
            lines = fin.read().strip(os.linesep).split(os.linesep)
            table = [l.split() for l in lines]

        headers = table[0]
        values = table[1:]

    ``values`` is a list of lists where each index is a row.

    Alternatively, you can use ``numpy.loadtxt`` and work from there:

    .. code:: python

        import numpy as np
        np.loadtxt('table.txt', dtype=str)

    Parameters
    ----------
    data : dict
        A dictionary where keys are the column names and values are the
        column values for each row. ``data`` dictionary values can be in
        the form of a list if there are multiple rows, or a single value
        if there is only one row.

    header : None or str, optional
        The header text to write before the table. Usually this is a commented
        text consisting of a small report or instructions. Defaults to
        ``None``, no header is written.

    column_spacing : str, optional
        The spacing between columns. Defaults to four empty spaces.

    fmt_args : arguments
        Named arguments to pass to
        :py:func:`haddock.gear.tables.format_value`.

    Returns
    -------
    str
        Table formatted string.

    See Also
    --------
    :py:func:`haddock.gear.tables.read_table_to_data_dict`
    :py:func:`haddock.gear.tables.convert_adjusted_width_to_tsv`
    :py:func:`haddock.gear.tables.convert_adjusted_width_to_csv`
    """
    header = "" if header is None else header
    col_names = list(data.keys())
    formatted_col_names = [None] * len(col_names)

    for i, col_name in enumerate(copy(col_names)):

        longest_ = max(
            len(value_to_str(value))
            for value in transform_to_list(data[col_name])
            )
        longest = max(longest_, len(col_name))

        formatted_col_names[i] = col_name.center(longest, " ")

        data[col_name] = [
            format_value(value, longest, " ", **fmt_args)
            for value in transform_to_list(data[col_name])
            ]

    col_name_row = column_spacing.join(formatted_col_names)

    # all columns should have the same lenght
    num_rows = len(data[col_names[0]])
    rows = [None] * num_rows
    for i in range(num_rows):
        rows[i] = []
        for col_name in col_names:
            try:
                rows[i].append(data[col_name][i])
            except IndexError as err:
                emsg = (
                    f"Column {col_name!r} as less items than the "
                    "expected number of rows."
                    )
                raise TableFormatError(emsg) from err
        rows[i] = column_spacing.join(rows[i])
    table_data = os.linesep.join(rows)

    # to situations to save tables with and header to avoid having extra
    # empty lines.
    if header:
        table = os.linesep.join([
            header,
            col_name_row,
            table_data,
            os.linesep,
            ])
    else:
        table = os.linesep.join([
            col_name_row,
            table_data,
            os.linesep,
            ])

    return table


def read_table_to_data_dict(fname, **kwargs):
    """
    Read a table to a data dictionary.

    Parameters
    ----------
    fname : str or :external:py:class:`pathlib.Path`
        Path to the table file.

    kwargs : named parameters
        As described in
        :py:func:`haddock.gear.tables.parse_table_to_data_dict`.

    Returns
    -------
    dict
        The data dictionary in the form of:

    .. code:: python

            {
                "col_name_1": [value1, value3],
                "col_name_2": [value2, value4],
                "col_name_3": [value5, value6],
                }

    See Also
    --------
    :py:func:`haddock.gear.tables.parse_table_to_data_dict`
    """
    txt = Path(fname).read_text()
    return parse_table_to_data_dict(txt)


def parse_table_to_data_dict(table_txt, comment="#", sep=None):
    """
    Read a adjusted with table to a data dictionary format.

    Does the inverse process of :py:func:`create_adjusted_col_width_table`.

    Given a table text::

        col_name_1    col_name_2    col_name_3
          value1        value2        value5
          value3        value4        valut6

    Returns a dictionary in the form of:

    .. code:: python

            {
                "col_name_1": [value1, value3],
                "col_name_2": [value2, value4],
                "col_name_3": [value5, value6],
                }

    Parameters
    ----------
    table_txt : str
        The string containing the table.

    comment : str
        Ignores line starting with comment character. Deafults to ``#``.

    sep : str
        Value separator. Defaults to ``None``, empty spaces.

    See Also
    --------
    :py:func:`haddock.gear.tables.create_adjusted_col_width_table`
    :py:func:`haddock.gear.tables.read_table_to_data_dict`
    """
    lines = table_txt.strip(os.linesep).split(os.linesep)
    table = [
        line.split(sep)
        for line in lines
        if line and not line.startswith(comment)
        ]
    table_dict = {}
    headers = table[0]
    print(headers)
    for col_num, header in enumerate(headers):
        table_dict[header] = [str_to_value(line[col_num]) for line in table[1:]]
    return table_dict


def str_to_value(value, missing_chars=('-',)):
    """
    Convert a string to a value.

    Implemented values:
        * int
        * float
        * None
        * nan
    """
    conversions = [
        (lambda x: x in ("None", "none"), none),
        (lambda x: x in missing_chars, none),
        (lambda x: x=="nan", nan),
        (is_str_int, int),
        (is_str_float, float),
        (true, str),
        ]

    for validate, func in conversions:
        if validate(value):
            return func(value)


def format_value(value, spacing, char, float_fmt="{:.3f}"):
    """
    Format the value for the table depending on its type.

    Usually this means applying decimal number reduction, left or right
    justifying depending on the type, and similar operatations.

    Parameters
    ----------
    value : anything
        The value to process.

    spacing : int
        The spacing of the cell.

    char : str
        The char to space the value.

    float_fmt : str
        The float format string. Defaults to ``{:.3f}``.

    Returns
    -------
    str
        The formatted value to fill the table.
    """
    if isinstance(value, float):
        return float_fmt.format(value).rjust(spacing, char)
    elif isinstance(value, (PDBFile, Path)) \
            or isinstance(value, str) and value.endswith(".pdb"):
        return value_to_str(value).ljust(spacing, " ")
    else:
        return value_to_str(value).center(spacing, char)


def value_to_str(value, float_fmt="{:.3f}"):
    """
    Convert value to string according to its type.

    Applies additional style and format convertion specific for
    HADDOCK3.

    Parameters
    ----------
    value : anything
        The value to process.

    Returns
    -------
    str
        The formatted value to fill the table.
    """
    if isinstance(value, float):
        return float_fmt.format(value)
    elif value is None:
        return "-"
    elif isinstance(value, PDBFile):
        return str(value.rel_path)
    elif isinstance(value, Path):
        return str(value)
    else:
        return str(value)


def convert_row_to_column_table(data):
    """
    Convert row based dictionaries to column based dictionaries.

    In other words, convert dictionaries in the form of:

    .. code:: python

        {
            1: {
                "col_name_1": value1,
                "col_name_2": value2,
                },

            2: {
                "col_name_1": value3,
                "col_name_2": value4,
                },
            }

    To dictionaries compatible with
    :py:func:`adjusted_col_width_table`.

    .. code:: python

        {
            "col_name_1": [value1, value3],
            "col_name_2": [value2, value4],
            }

    Returns
    -------
    dict
        The parsed data dictionary.

    See Also
    --------
    :py:func:`haddock.gear.tables.create_adjusted_col_width_table`
    """
    data2 = {}
    for row in data:
        for element in data[row]:
            value = data[row][element]
            values = data2.setdefault(element, [])
            values.append(value)
    return data2


def convert_sep_to_sep(fname, fout=None, in_sep=None, out_sep="\t", comment="#"):
    """
    Convert a text based table from one sep to another.

    For example, convert a space separated table to a tab separated
    table.

    Parameters
    ----------
    fname : str or :external:py:class:`pathlib.Path`
        The path to the table file.

    fout : str or :external:py:class:`pathlib.Path`
        If provided, saves the table to the file path specified in
        ``fout``. If ``None`` returns the table string.

    comment : str
        Reproduces lines starting with comment char unmodified.
        Defaults to ``#``.

    in_sep : str
        The value separator for the input table. Defaults to ``None``,
        empty spaces.

    out_sep : str
        The value separator for the output table. Defaults to empty spaces.

    Returns
    -------
    str
        The new table text if ``fout`` is set to ``None``.

    None
        ``None``. Writes the table to a file if ``fout`` is given.

    See Also
    --------
    :py:func:`haddock.gear.tables.convert_ssc_tsv`
    :py:func:`haddock.gear.tables.convert_ssc_csv`
    """
    table_text = Path(fname).read_text()
    lines = table_text.strip(os.linesep).split(os.linesep)
    new_text = []
    for line in lines:
        if line.startswith('#'):
            new_text.append(line)
        else:
            new_text.append(out_sep.join(line.split(in_sep)))

    return os.linesep.join(new_text)


def convert_ssc_tsv(*args, **kwargs):
    r"""
    Convert a space-separated table to a tab-separated table.

    Uses :py:func:`haddock.gear.convert_sep_to_sep` with ``in_sep=" "``
    and ``out_sep="\t"``. Other parameters are as explained there.
    """
    return convert_sep_to_sep(*args, in_sep=None, out_sep="\t", **kwargs)


def convert_ssc_csv(*args, **kwargs):
    """
    Convert a space-separated table to a comma-separated table.

    Uses :py:func:`haddock.gear.convert_sep_to_sep` with ``in_sep=" "``
    and ``out_sep=","``. Other parameters are as explained there.
    """
    return convert_sep_to_sep(*args, in_sep=None, out_sep=",", **kwargs)
