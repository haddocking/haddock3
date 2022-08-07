"""
Convert HADDOCK3 config dictionaries to user configuration files.

The functions implemented here have a general character. For them to be
functional within the HADDOCK3 scope you need to provide additional
input arguments. All the details are explained in the function's
docstrings (help).

If you look for equal functions but already pre-prepared for HADDOCK3
modules, see the functions:

* :py:func:`haddock.modules.convert_config`
* :py:func:`haddock.modules.save_config`
"""
import collections.abc
import os
from copy import deepcopy
from pathlib import Path

import toml

from haddock import EmptyPath
from haddock.gear.config_reader import get_module_name, _main_quoted_header_re, _sub_quoted_header_re


def save_config(params, path):
    """
    Save a dictionary to a HADDOCK3 config file.

    Parameters
    ----------
    params : dict
        The dictionary containing the parameters.

    path : str or pathlib.Path
        File name where to save the configuration file.
    """
    params = recursive_convert_paths_to_strings(params)

    #ostr = os.linesep.join(convert_config(params, **kwargs))
    cfg_str = toml.dumps(params).split(os.linesep)

    new_lines = []
    for line in cfg_str:
        if group := _main_quoted_header_re.match(line):
            name = group[1]
            new_line = f"[{name}]"
        elif group := _sub_quoted_header_re.match(line):
            name = group[1]
            new_line = f"[{name}{group[2]}]"
        else:
            new_line = line
        new_lines.append(new_line)

    cfg_str = os.linesep.join(new_lines)

    with open(path, "w") as fout:
        fout.write(cfg_str)


def recursive_convert_paths_to_strings(params):
    for param, value in params.items():
        if isinstance(value, (Path, EmptyPath)):
            params[param] = str(value)
        elif isinstance(value, collections.abc.Mapping):
            params[param] = recursive_convert_paths_to_strings(value)
        elif isinstance(value, (tuple, list)):
            for i, v in enumerate(value):
                if isinstance(v, (Path, EmptyPath)):
                    value[i] = str(v)
            params[param] = value
            print(value)

    return params
