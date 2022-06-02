"""
Tools to help the different modules creating pretty tables.

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

See:

* :py:func:`create_human_readable_table`
* :py:func:`convert_row_to_column_table`
"""
import os
from copy import copy
from pathlib import Path

from haddock.core.exceptions import HaddockError
from haddock.libs.libontology import PDBFile
from haddock.libs.libutil import transform_to_list


class TableFormatError(HaddockError):
    """Error related to table format."""

    pass


def create_human_readable_table(
        data,
        header=None,
        column_spacing="    ",
        ):
    """
    Create a human-readable table.

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

    Parameters
    ----------
    data : dict
        A dictionary where keys are the column names and values are the
        column values for each row. ``data`` dictionary values can be in
        the form of a list if there are multiple rows, or a single value
        if there is only one row.

    header : None or str, optional
        The header to write before the column. Usually this is a
        commented text with a small report or any instructions. Defaults
        to ``None``, no header is written.

    column_spacing : str, optional
        The spacing between columns. Defaults to four empty spaces.

    Returns
    -------
    str
        Table formatted string.
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
            format_value(value, longest, " ")
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


def format_value(value, spacing, char):
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

    Returns
    -------
    str
        The formatted value to fill the table.
    """
    if isinstance(value, float):
        return f"{value:.3f}".rjust(spacing, char)
    elif isinstance(value, (PDBFile, Path)) \
            or isinstance(value, str) and value.endswith(".pdb"):
        return value_to_str(value).ljust(spacing, " ")
    else:
        return value_to_str(value).center(spacing, char)


def value_to_str(value):
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
        return f"{value:.3f}"
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
    :py:func:`create_human_readable_table`.

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
    :py:func:`create_human_readable_table`
    """
    data2 = {}
    for row in data:
        for element in data[row]:
            value = data[row][element]
            values = data2.setdefault(element, [])
            values.append(value)
    return data2
