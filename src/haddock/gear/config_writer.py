"""Writes to HADDOCK3 config files."""
import collections.abc
import os
from copy import deepcopy
from pathlib import Path

from haddock import EmptyPath
from haddock.gear.config_reader import get_module_name
from haddock.modules.definitions import (
    modules_category,
    non_mandatory_general_parameters_defaults,
    )


def convert_config(d, ignore_params=None):
    """
    Convert parameter dictionary to a HADDOCK3 user configuration file.

    This function is a generator.

    Examples
    --------
    >>> gen = convert_config(d)
    >>> text = os.linesep.join(gen)
    >>> with open("params.cfg", "w") as fout:
    >>>     fout.write(text)

    Parameters
    ----------
    d : dictionary
        The dictionary containing the parameters.

    ignore_params : set, list, tuple
        List of parameters to be ignored. These parameters won't be
        yielded.

    Yields
    ------
    str
        Line by line for the HADDOCK3 user configuration file.
    """
    ignore_params = set() or ignore_params

    valid_keys = (
        (k, v)
        for k, v in d.items()
        if k not in non_mandatory_general_parameters_defaults
        )

    # place parameters that are dictionaries at the end of the config file
    # this is important because toml-like config are to be read line by line
    # this was made first for compatibility with the topoaa `molecules` and
    # `mol1` like parameters, which are dictionaries. Note `molecules` and
    # `mol*` start with `mol`
    sorted_keys = sorted(valid_keys, key=_is_dict)

    for key, value in sorted_keys:

        if get_module_name(key) in modules_category:
            yield f"[{get_module_name(key)}]"
            yield from convert_config(value, ignore_params=ignore_params)
            continue

        elif isinstance(d[key], collections.abc.Mapping):
            yield f"[{key}]"
            yield from convert_config(d[key])
            continue

        elif isinstance(value, (list, tuple)):
            yield f"{key} = ["
            yield from _list_by_value(value)
            yield "    ]"
            # multiline lists in haddock3 configuration files
            # need to be followed by an empty line.
            yield os.linesep
            continue

        yield _convert_value_to_config_string(value, key + " = {}")


def _is_dict(t):
    """
    Return 1 if value is type Mapping (dict) or 0 otherwise.

    Parameters
    ----------
    t : tuple
        key and value pair.
    """
    return 1 if isinstance(t[1], collections.abc.Mapping) else 0


def _list_by_value(values):
    for value in values[:-1]:
        yield _convert_value_to_config_string(value, fmt="    {},")

    yield _convert_value_to_config_string(values[-1], fmt="    {}")


def _convert_value_to_config_string(value, fmt):
    """Convert a value to its string representation in an HADDOCK3 config file."""
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
            f"{value!r} of type: {type(value)!r}"
            )


def save_config(params, path, module_name=None):
    """
    Save a dictionary to a HADDOCK3 config file.

    Parameters
    ----------
    params : dict
        The dictionary containing the parameters.

    path : str or pathlib.Path
        File name where to save the configuration file.

    module_name : str or None
        Provide a name to the module to add the module's name to the
        output CFG file. This applies for cases where the parameters
        dictionary is a flat dictionary defining only the parameters of
        a specific module without referring to the module's name.
    """
    if module_name:
        params = {module_name: deepcopy(params)}

    ostr = os.linesep.join(convert_config(params))
    with open(path, "w") as fout:
        fout.write(ostr)
