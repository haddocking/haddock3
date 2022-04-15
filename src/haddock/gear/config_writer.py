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

        elif isinstance(d[key], collections.abc.Mapping):
            yield f"[{key}]"
            yield from write_config(d[key])

        # TODO: this has to be recursive in multiline
        # and consider paths, emptypaths, strings, and numbers
        elif isinstance(value, (list, tuple)):
            yield f"{key} = {list(map(str, value))!r}"

        elif isinstance(value, (EmptyPath, Path)):
            yield f"{key} = {str(value)!r}"

        else:
            yield f"{key} = {value!r}"


def save_config(params, path):
    """Save a dictionary to a HADDOCK3 config file."""
    with open(path, "w") as fout:
        _ = (line for line in write_config(params))
        ostr = os.linesep.join(_)
        fout.write(ostr)
