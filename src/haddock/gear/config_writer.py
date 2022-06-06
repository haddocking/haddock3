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

from haddock import EmptyPath
from haddock.gear.config_reader import get_module_name


def convert_config(
        params,
        ignore_params=None,
        module_names=None,
        _module_key=None,
        ):
    """
    Convert parameter dictionary to a HADDOCK3 user configuration file.

    This function is a generator.

    Examples
    --------
    >>> gen = convert_config(parameters_dict)
    >>> text = os.linesep.join(gen)
    >>> with open("params.cfg", "w") as fout:
    >>>     fout.write(text)

    Note you might want to give the list of the available module names,
    and remove the non mandatory general parameters.

    >>> from haddock.modules import modules_category
    >>> from haddock.modules import non_mandatory_general_parameters_defaults
    >>> gen = convert_config(
        parameters_dict,
        module_names=modules_category,
        ignore_param=non_mandatory_general_parameters_defaults,
        )
    >>> text = os.linesep.join(gen)
    >>> with open("params.cfg", "w") as fout:
    >>>     fout.write(text)

    Parameters
    ----------
    params : dictionary
        The dictionary containing the parameters.

    ignore_params : set, list, tuple
        List of parameters to be ignored. These parameters won't be
        yielded.

    module_names : set, list, tuple
        A list with the names of possible headings. Headings usually are
        module names.

    _module_key : internal parameter
        This function uses itself recursively. `_module_key` is for
        internal usage only.

    Yields
    ------
    str
        Line by line for the HADDOCK3 user configuration file.

    See Also
    --------
    :py:func:`haddock.modules.convert_config`.
    """
    # can't do
    # ignore_params = set() or ignore_params
    # because set() also evaluates to false
    if ignore_params is None:
        ignore_params = set()

    if module_names is None:
        module_names = set()

    valid_keys = (
        (k, v)
        for k, v in params.items()
        if k not in ignore_params
        )

    # place parameters that are dictionaries at the end of the config file
    # this is important because toml-like config are to be read line by line
    # this was made first for compatibility with the topoaa `molecules` and
    # `mol1` like parameters, which are dictionaries. Note `molecules` and
    # `mol*` start with `mol`
    sorted_keys = sorted(valid_keys, key=_is_dict)

    for i, (key, value) in enumerate(sorted_keys):

        module_key = get_module_name(key)
        if module_key in module_names:

            # used to get an extra line when writing the second module
            # header
            if i > 0:
                yield ""

            yield f"[{module_key}]"
            yield from convert_config(
                value,
                ignore_params=ignore_params,
                module_names=module_names,
                _module_key=module_key,
                )
            continue

        elif isinstance(value, collections.abc.Mapping):
            if _module_key is None:
                yield ""
                yield f"[{key}]"
            else:
                # module's subkeys should follow without an extra
                # newline
                yield f"[{_module_key}.{key}]"
            yield from convert_config(value)
            continue

        elif isinstance(value, (list, tuple)):
            if not value:
                yield f"{key} = []"
            else:
                yield f"{key} = ["
                yield from _list_by_value(value)
                yield "    ]"
                # multiline lists in haddock3 configuration files
                # need to be followed by an empty line.
                yield ""
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
    """Convert values in a list."""
    if not values:
        # empty lists are address at convert_config
        raise AssertionError("Empty lists should enter this function.")

    for value in values[:-1]:
        yield _convert_value_to_config_string(value, fmt="    {},")

    # note the comma in the `fmt` parameter!
    yield _convert_value_to_config_string(values[-1], fmt="    {}")


def _convert_value_to_config_string(value, fmt):
    """Convert a value to its string representation in an config file."""
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

    # nan values also enter here because they are floats
    # and get converted to nan strings and saved to the config without quotes
    elif isinstance(value, (int, float)):
        return fmt.format(value)

    elif value is None:
        return fmt.format("none")

    else:
        raise AssertionError(
            "There should be no values of this type: "
            f"{value!r} of type: {type(value)!r}"
            )


def save_config(params, path, module_name=None, **kwargs):
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

    kwargs
        Other parameters are as defined in `convert_config`. Take special
        attention to `module_names` parameter.

    See Also
    --------
    :py:func:`haddock.modules.save_config`.
    """
    if module_name:
        if not isinstance(module_name, str):
            raise TypeError("`module_name` must be str: {type(module_name)} given") # noqa E501
        params = {module_name: deepcopy(params)}
        # add module name to module names parameter
        mn = kwargs.setdefault("module_names", [])
        kwargs["module_names"] = list(mn) + [module_name]

    ostr = os.linesep.join(convert_config(params, **kwargs))
    with open(path, "w") as fout:
        fout.write(ostr)
