"""
Adjusted column width tables.

These functions create a tables of variable column width to accommodate the
width of the different column values. In this way, columns are easily
redable by humans and easily parseable by computers.
"""
import os


def create_adjusted_col_width_table(
        data,
        header=None,
        column_spacing="    ",
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

    Returns
    -------
    str
        Table formatted string.

    See Also
    --------
    :py:func:`haddock.gear.tables.read_table_to_data_dict`
    :py:func:`haddock.gear.tables.convert_adjusted_with_to_tsv`
    :py:func:`haddock.gear.tables.convert_adjusted_with_to_csv`
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


def read_table_to_data_dict(fname, comment="#", sep=" "):
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
    fname : str or :external:py:class:`pathlib.Path`
        The path to the table file.

    comment : str
        Ignores line starting with comment character. Deafults to ``#``.

    sep : str
        Value separator. Defaults to empty spaces.
    """
    table_text = Path(fname).read_text()
    lines = table_text.read().strip(os.linesep).split(os.linesep)
    table = [l.split(sep) for l in lines if not l.startwith(comment)]
    table_dict = {}
    headers = table[0]
    for col_num, header in enumerate(headers):
        table_dict[header] = [line[col_num] for line in table[1:]]
    return table_dict
