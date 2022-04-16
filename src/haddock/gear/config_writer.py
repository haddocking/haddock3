"""Writes to HADDOCK3 config files."""
import collections.abc
import os
from pathlib import Path

from haddock import EmptyPath
from haddock.gear.config_reader import get_module_name
from haddock.modules import modules_category


def write_config(d):
    """Write a config dictionary to a HADDOCK3 config."""
    for key, value in d.items():
        if get_module_name(key) in modules_category:
            yield f"[{get_module_name}]"
            continue

        elif isinstance(d[key], collections.abc.Mapping):
            yield f"[{key}]"
            yield from write_config(d[key])
            continue

        # TODO: this has to be recursive in multiline
        # and consider paths, emptypaths, strings, and numbers
        elif isinstance(value, (list, tuple)):
            yield f"{key} = ["
            yield from _list_by_value(value)
            yield "    ]"
            # multiline lists in haddock3 configuration files
            # need to be followed by an empty line.
            yield os.linesep
            continue

        else:
            yield _convert_value_to_config_string(value, key + " = {}")
            continue

        emsg = (
            f"Can't convert {key!r} and value {value} "
            "to HADDOCK3 config parameter. "
            "Code shouldn't reach this state."
            )
        raise AssertionError(emsg)


def _list_by_value(values):
    for value in values[:-1]:
        yield _convert_value_to_config_string(value, fmt="    {},")

    yield _convert_value_to_config_string(values[-1], fmt="    {}")


def _convert_value_to_config_string(value, fmt):
    """Convert a value to its string representation for HADDOCK3 config file."""
    if isinstance(value, bool):
        value = str(value).lower()
        return fmt.format(value)

    elif isinstance(value, str):
        value = '"' + value + '"'
        return fmt.format(value)

    elif isinstance(value, Path):
        value = '"' + str(value.resolve()) + '"'
        return fmt.format(value)

    elif isinstance(value, EmptyPath):
        return fmt.format('""')

    elif isinstance(value, (int, float)):
        return fmt.format(value)

    else:
        raise AssertionError(
            "There should be no values of this type: "
            f"{type(value)}"
            )


def save_config(params, path):
    """Save a dictionary to a HADDOCK3 config file."""
    with open(path, "w") as fout:
        _ = (line for line in write_config(params))
        ostr = os.linesep.join(_)
        fout.write(ostr)
