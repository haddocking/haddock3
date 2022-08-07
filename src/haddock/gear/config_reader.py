"""
Implementation of a TOML-like configuration file for HADDOCK3.

HADDOCK3 user configuration files follow TOML syntax plus additional
features required for HADDOCK3. Therefore, we have implemented our own
TOML-like parser to accommodate the extra features needed for HADDOCK3.

The most relevant features of HADDOCK3 user configuration files are:

**Accepts repeated keys:** blocks can have repeated names, for example:

.. code:: toml

        [topoaa]
        # parameters here...

        [rigidbody]
        # parameters here...

        [caprieval]
        # parameters here...

        [flexref]
        # parameters here...

        [caprieval]
        # parameters here...

**Allows in-line comments:**

.. code:: toml

    [module]
    parameter = value # some comment

**The following types are implemented**:

* strings
* numbers
* null/none
* ``nan``
* lists defined in one lines
* lists defined in multiple lines (require empty line after the list
  definition)
* `datetime.fromisoformat`

**To efficiently use this module, see functions:**

* :py:func:`haddock.gear.config_reader.read_config`
* :py:func:`haddock.gear.config_reader.get_module_name`
"""
import collections.abc
import os
import re
from contextlib import suppress
from pathlib import Path

import toml

from haddock import EmptyPath
from haddock.core.defaults import RUNDIR
from haddock.core.exceptions import ConfigurationError


# the re.ASCII parameter makes sure non-ascii chars are rejected in the \w key

# Captures the main headers.
# https://regex101.com/r/9urqti/1
_main_header_re = re.compile(r'^ *\[(\w+)\]', re.ASCII)

# regex by @sverhoeven
# Matches ['<name>.<digit>']
_main_quoted_header_re = re.compile(r'^ *\[\'(\w+)\.\d+\'\]', re.ASCII)

# regex by @sverhoeven
_sub_quoted_header_re = re.compile(
    r'^ *\[\'(\w+)\.\d+\'((?:\.\w+)+)\]',
    re.ASCII,
    )

# Captures sub-headers
# https://regex101.com/r/6OpJJ8/1
# thanks https://stackoverflow.com/questions/39158902
_sub_header_re = re.compile(r'^ *\[(\w+)((?:\.\w+)+)\]', re.ASCII)

# Some parameters are paths. The configuration reader reads those as
# pathlib.Paths so that they are injected properly in the rest of the
# workflow. In general, any parameter ending with `_fname` is a path,
# but there are also other parameters that are paths. Those need to be
# added to this list bellow:
_keys_that_accept_files = [
    "cns_exec",
    "executable",
    RUNDIR,
    ]


class ConfigFormatError(Exception):
    """Exception if there is a format error."""

    pass


class DuplicatedParameterError(Exception):
    """Exception if duplicated parameters are found."""

    pass


# main public API
def read_config(fpath):
    """
    Read HADDOCK3 configure file to a dictionary.

    Parameters
    ----------
    fpath : str or :external:py:class:`pathlib.Path`
        Path to user configuration file.

    Returns
    -------
    dictionary
        Representing the user configuration file where first level of
        keys are the module names. Step keys will have a numeric
        suffix, for example: ``module.1``.

    .. see-also::
        * :py:func:`read_config_str`
    """
    try:
        return read_config_str(Path(fpath).read_text())
    except Exception as err:
        raise ConfigurationError('Something is wrong with the config file.') from err  # noqa: E501


def read_config_str(cfg_str):
    """
    Read a string representing a config file to a parameter dictionary.

    Config strings are converted to toml-compatible format and finally
    read by the toml library.

    Parameters
    ----------
    cfg_str : str
        The string reprensenting the config file. Accepted formats are
        the HADDOCK3 config file or pure `toml` syntax.

    Returns
    -------
    dict
        The config in the form of a dictionary.
        The order of the keys matters as it defines the order of the
        steps in the workflow.
    """
    # first, attempts to read toml directly
    with suppress(toml.TomlDecodeError):
        return toml.loads(cfg_str)

    new_lines = []
    cfg_lines = cfg_str.split(os.linesep)

    counter = {}

    for line in cfg_lines:

        if group := _main_header_re.match(line):
            name = group[1]
            counter.setdefault(name, 0)
            counter[name] += 1
            count = counter[name]
            new_line = f"['{name}.{count}']"

        elif group := _main_quoted_header_re.match(line):
            name = group[1]
            counter.setdefault(name, 0)
            counter[name] += 1
            count = counter[name]
            new_line = f"['{name}.{count}']"

        elif group := _sub_header_re.match(line):
            name = group[1]
            count = counter[name]  # name should be already defined here
            new_line = f"['{name}.{count}'{group[2]}]"

        elif group := _sub_quoted_header_re.match(line):
            name = group[1]
            count = counter[name]  # name should be already defined here
            new_line = f"['{name}.{count}'{group[2]}]"

        else:
            new_line = line

        new_lines.append(new_line)

    cfg = os.linesep.join(new_lines)

    cfg = toml.loads(cfg)

    cfg = convert_variables_to_paths(cfg)

    return cfg


def convert_variables_to_paths(cfg):
    """
    Convert config variables to paths.

    Applies only to those variables HADDOCK3 needs as paths.
    """
    for param, value in cfg.items():
        if param == "molecules":
            cfg[param] = [convert_to_path(v) for v in value]

        elif match_path_criteria(param):
            cfg[param] = convert_to_path(value)

        elif isinstance(value, collections.abc.Mapping):
            cfg[param] = convert_variables_to_paths(value)

    return cfg


def convert_to_path(value):
    """Convert a variable to Path if variable evaluates to true."""
    if value:
        return Path(value)
    else:
        return EmptyPath()


def match_path_criteria(param):
    """Confirm a parameter should be converted to path."""
    return param.endswith('_fname') or param in _keys_that_accept_files


def get_module_name(name):
    """Get the name according to the config parser."""
    return name.split('.')[0]


def write_config(params, path, toml=False):
    """
    Write a dictionary to a HADDOCK3 config file.

    Parameters
    ----------
    params : dict
        The dictionary containing the parameters.

    path : str or pathlib.Path
        File name where to save the configuration file.

    toml : bool
        Save config in pure toml format.
    """
    params = recursive_convert_paths_to_strings(params)

    if toml:
        with open(path, 'w') as fout:
            toml.dump(params, fout)
        return

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
    """Convert paths to strings recursively over a dictionary."""
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

    return params
