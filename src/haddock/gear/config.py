"""
Module to read and write HADDOCK3 configuration files.

HADDOCK3 user configuration files follow TOML syntax plus additional
features required for HADDOCK3. Therefore, we have implemented our own
TOML-like parser to accommodate the extra features needed for HADDOCK3.
However, HADDOCK3 can also read pure TOML files.

The most relevant features of HADDOCK3 user configuration files are the
hability to prepare TOML-like files but with repeated header names, for
example:

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

From the user perspective, everything else just behaves as in TOML.
Actually, HADDOCK3 uses TOML to finally read the configuration files.
Users can also provide pure TOML compatible files, for example:

.. code:: toml

        [topoaa]
        # parameters here...

        [rigidbody]
        # parameters here...

        [caprieval]
        # parameters here...

        [flexref]
        # parameters here...

        ['caprieval.1']
        # parameters here...

Because HADDOCK3 uses TOML in the end to read the configuration files,
all TOML rules apply regarding defining ``parameter = value`` pairs.

Read more about TOML for Python here: https://pypi.org/project/toml/
"""

import collections.abc
import os
import re
from pathlib import Path

import toml

from haddock import EmptyPath
from haddock.core.defaults import RUNDIR
from haddock.core.exceptions import ConfigurationError
from haddock.core.typing import FilePath, ParamDict, ParamMap, ParamMapT, Union
from haddock.libs.libutil import (
    recursive_convert_paths_to_strings,
    transform_to_list,
)


# the re.ASCII parameter makes sure non-ascii chars are rejected in the \w key

# Captures the main headers.
# https://regex101.com/r/9urqti/1
_main_header_re = re.compile(r"^ *\[(\w+)\]", re.ASCII)

# regex by @sverhoeven
# Matches ['<name>.<digit>']
_main_quoted_header_re = re.compile(r"^ *\[\'(\w+)\.\d+\'\]", re.ASCII)

# Matches a quoted bare header ['<name>'] (pure TOML, no numeric suffix).
# Normalized like a main header so all instances of a module carry a `.N`.
_main_quoted_bare_header_re = re.compile(r"^ *\[\'(\w+)\'\]", re.ASCII)

# Captures sub-headers
# https://regex101.com/r/6OpJJ8/1
# thanks https://stackoverflow.com/questions/39158902
_sub_header_re = re.compile(r"^ *\[(\w+)((?:\.\w+)+)\]", re.ASCII)

# regex by @sverhoeven
_sub_quoted_header_re = re.compile(
    r"^ *\[\'(\w+)\.\d+\'((?:\.\w+)+)\]",
    re.ASCII,
)

# Captures parameter uppercase boolean
_uppercase_bool_re = re.compile(r"(_?\w+((_?\w+?)+)?\s*=\s*)(True|False)", re.ASCII)

# Some parameters are paths. The configuration reader reads those as
# pathlib.Path so that they are injected properly in the rest of the
# workflow. In general, any parameter ending with `_fname` is a path,
# but there are also other parameters that are paths. Those need to be
# added to this list bellow:
_keys_that_accept_files = [
    "cns_exec",
    "executable",
    RUNDIR,
]


class ConfigFormatError(ConfigurationError):
    """Exception raised when incompatible config header formats are mixed."""

    pass


# main public API
def load(fpath: FilePath, strict: bool = True) -> ParamDict:
    """
    Load an HADDOCK3 configuration file to a dictionary.

    Accepts HADDOCK3 ``cfg`` files or pure ``toml`` files.

    Parameters
    ----------
    fpath : str or :external:py:class:`pathlib.Path`
        Path to user configuration file.

    strict : bool
        When ``True`` (default), reject numeric dotted headers such as
        ``[module.1]`` which silently merge into sub-tables and drop
        workflow steps silently. Set to ``False`` only when reloading internally
        generated ``params.cfg`` files, which legitimately contain such
        numeric sub-tables. See :py:func:`loads`.

    Returns
    -------
    dictionary
        Representing the user configuration file where first level of
        keys are the module names. Step keys will have a numeric
        suffix, for example: ``module.1``.

    .. see-also::
        * :py:func:`loads`
    """
    try:
        return loads(Path(fpath).read_text(), strict=strict)
    except ConfigFormatError:
        # raise instead of masking
        raise
    except Exception as err:
        raise ConfigurationError("Something is wrong with the config file.") from err  # noqa: E501


def loads(cfg_str: str, strict: bool = True) -> ParamDict:
    """
    Read a string representing a config file to a dictionary.

    Config strings are converted to toml-compatible format and finally
    read by the toml library.

    All headers (dictionary keys) will be suffixed by an integer
    starting at ``1``. For example: ``topoaa.1``. If the key is
    repeated, ``2`` is appended, and so forth. Even if specific
    integers are provided by the user, the suffix integers will be
    normalized.

    Parameters
    ----------
    cfg_str : str
        The string representing the config file. Accepted formats are
        the HADDOCK3 config file or pure `toml` syntax.

    strict : bool
        When ``True`` (default), a numeric dotted header such as
        ``[module.1]`` raises :py:class:`ConfigFormatError`. Repeated
        HADDOCK3 modules must be declared with repeated bare headers
        (``[module]`` written several times); the dotted form is a TOML
        sub-table and would silently drop the intended step. Set to
        ``False`` only when reloading internally generated ``params.cfg``
        files, which legitimately store numeric sub-tables.

    Returns
    -------
    all_configs : dict
        A dictionary holding all the configuration file steps:

        - 'raw_input': Original input file as provided by user.
        - 'cleaned_input': Regex cleaned input file.
        - 'loaded_cleaned_input': Dict of toml loaded cleaned input.
        - 'final_cfg': The config in the form of a dictionary. In which
          the order of the keys matters as it defines the order of the
          steps in the workflow.
    """
    new_lines: list[str] = []
    cfg_lines = cfg_str.split(os.linesep)
    counter: dict[str, int] = {}
    styles: dict[str, str] = {}

    def register_style(name: str, style: str) -> None:
        """Reject mixing numbered and unnumbered instances of a module."""
        previous = styles.setdefault(name, style)
        if strict and previous != style:
            raise ConfigFormatError(
                f"Inconsistent headers for module '{name}': a numbered "
                f"instance (['{name}.N']) is mixed with an unnumbered one "
                f"([{name}] or ['{name}']). Use one style consistently for "
                "every instance of a module."
            )

    # this for-loop normalizes all headers regardless of their input format.
    for line in cfg_lines:
        if group := _main_header_re.match(line):
            name = group[1]
            register_style(name, "bare")
            counter.setdefault(name, 0)
            counter[name] += 1
            count = counter[name]
            new_line = f"['{name}.{count}']"

        elif group := _main_quoted_header_re.match(line):
            name = group[1]
            register_style(name, "numbered")
            counter.setdefault(name, 0)
            counter[name] += 1
            count = counter[name]
            new_line = f"['{name}.{count}']"

        elif group := _main_quoted_bare_header_re.match(line):
            name = group[1]
            register_style(name, "bare")
            counter.setdefault(name, 0)
            counter[name] += 1
            count = counter[name]
            new_line = f"['{name}.{count}']"

        elif group := _sub_header_re.match(line):
            name = group[1]
            # A numeric first sub-segment (e.g. `[caprieval.1]`) is the
            # ambiguous mix of the `.cfg` repeated-header format and TOML
            # NOTE: For some reason internally haddock will generate this
            # mix of toml-cfg - here we reject it unless we are reloading
            # an internally generated `params.cfg` (strict=False)
            if strict and group[2].split(".")[1].isdigit():
                raise ConfigFormatError(
                    f"Ambiguous config header '[{name}{group[2]}]': numeric "
                    "dotted headers are not valid HADDOCK3 module syntax. "
                    "To declare repeated modules use bare repeated headers: "
                    f"'[{name}]' written several times or '[\"{name}.N\"]'"
                )
            count = counter[name]  # name should be already defined here
            new_line = f"['{name}.{count}'{group[2]}]"

        elif group := _sub_quoted_header_re.match(line):
            name = group[1]
            count = counter[name]  # name should be already defined here
            new_line = f"['{name}.{count}'{group[2]}]"

        elif group := _uppercase_bool_re.match(line):
            param = group[1]  # Catches 'param = '
            uppercase_bool = group[4]
            new_line = f"{param}{uppercase_bool.lower()}"  # Lowercase bool

        else:
            new_line = line

        new_lines.append(new_line)

    # Re-build workflow configuration file
    cfg = os.linesep.join(new_lines)

    try:
        cfg_dict = toml.loads(cfg)  # Try to load it with the toml library
    except Exception as err:
        raise ConfigurationError(
            f"Some thing is wrong with the config file: {str(err)}"
        ) from err

    final_cfg = convert_variables_to_paths(cfg_dict)

    all_configs = {
        "raw_input": cfg_str,
        "cleaned_input": cfg,
        "loaded_cleaned_input": cfg_dict,
        "final_cfg": final_cfg,
    }

    # Return all configuration steps
    return all_configs


def convert_variables_to_paths(cfg: ParamMapT) -> ParamMapT:
    """
    Convert config variables to paths.

    Applies only to those variables HADDOCK3 needs as paths.

    Parameters
    ----------
    cfg : dictionary
        A dictionary containing the HADDOCK3 parameters.

    .. see-also::
        :py:func:`match_path_criteria`
    """
    for param, value in cfg.items():
        if param == "molecules":
            # `transform_to_list` is needed to allow strings OR lists in the
            # `molecule` parameter
            cfg[param] = [convert_to_path(v) for v in transform_to_list(value)]

        elif match_path_criteria(param):
            cfg[param] = convert_to_path(value)

        elif isinstance(value, collections.abc.Mapping):
            cfg[param] = convert_variables_to_paths(value)  # type: ignore

    return cfg


def convert_to_path(value: FilePath) -> Union[Path, EmptyPath]:
    """Convert a variable to Path if variable evaluates to true."""
    if value:
        return Path(value)
    else:
        return EmptyPath()


def match_path_criteria(param: str) -> bool:
    """
    Confirm a parameter should be converted to path.

    This function contains the rules defining which parameter need to be
    converted to Path before be sent to HADDOCK3.
    """
    return param.endswith("_fname") or param in _keys_that_accept_files


def get_module_name(name: str) -> str:
    """
    Get the module name according to the config parser.

    Get the name without the integer suffix.
    """
    return name.split(".")[0]


def save(params: ParamMap, path: FilePath, pure_toml: bool = False) -> None:
    """
    Write a dictionary to a HADDOCK3 config file.

    Write the HADDOCK3 parameter dictionary to a `.cfg` file. There is
    also the option to write in pure TOML format. Both are compatible with
    HADDOCK3.

    Parameters
    ----------
    params : dict
        The dictionary containing the parameters.

    path : str or pathlib.Path
        File name where to save the configuration file.

    pure_toml : bool
        Save config in pure toml format.
    """
    params = recursive_convert_paths_to_strings(params)

    if pure_toml:
        with open(path, "w") as fout:
            toml.dump(params, fout)
        return

    cfg_str = toml.dumps(params).split(os.linesep)

    new_lines: list[str] = []
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

    cfg_str = os.linesep.join(new_lines)  # type: ignore

    with open(path, "w") as fout:
        fout.write(cfg_str)  # type: ignore
