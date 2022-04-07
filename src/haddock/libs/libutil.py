"""General utilities."""
import collections.abc
import contextlib
import re
import shutil
import subprocess
import sys
from copy import deepcopy
from functools import partial
from math import ceil, log10
from os import cpu_count
from pathlib import Path

from haddock import log
from haddock.core.exceptions import SetupError
from haddock.gear.greetings import get_goodbye_help


check_subprocess = partial(
    subprocess.run,
    shell=True,
    check=True,
    stdout=subprocess.DEVNULL,
    )


def get_result_or_same_in_list(function, value):
    """
    Return the result if True or the value within a list.

    Applies `function` to `value` and returns its result if it evaluates
    to True. Otherwise, return the value within a list.

    `function` should receive a single argument, the `value`.
    """
    result = function(value)
    return result if result else [value]


def make_list_if_string(item):
    """Put `item` into a list."""
    if isinstance(item, str):
        return [item]
    return item


def transform_to_list(item):
    """
    Put `item` into a list if not a list already.

    If it is set, transforms the set into a list.

    If it is a dict, returns a list of the keys.

    If it is tuple, returns the tuple.

    If a list, returns the same.

    Everything else returns `item` inside a one element list.
    """
    if isinstance(item, set):
        return list(item)

    if isinstance(item, dict):
        return list(item.keys())

    if isinstance(item, (list, tuple)):
        return item

    return [item]


def copy_files_to_dir(paths, directory):
    """
    Copy files to directory.

    Parameters
    ----------
    paths : iterable of paths
        Source files.

    directory : path
        Where to copy files to.
    """
    for path in paths:
        shutil.copy(path, directory)


def zero_fill(number, digits=2):
    """Make a number string zero filled to the left."""
    return str(number).zfill(digits)


def remove_folder(folder):
    """Remove a folder if it exists."""
    if folder.exists():
        log.warning(f'{folder} exists and it will be REMOVED!')
        shutil.rmtree(folder)


def remove_dict_keys(d, keys):
    """
    Remove `keys` from dictionary (`d`).

    Return
    ------
    dict
        A copy of `d` dictionary without the `keys`.
    """
    return {k: deepcopy(v) for k, v in d.items() if k not in keys}


def parse_ncores(n=None, njobs=None, max_cpus=None):
    """
    Check the number of cores according to HADDOCK3 architecture.

    Parameters
    ----------
    n : int or str
        The desired number of cores. If `None` is given, returns the
        maximum number of cores allowed, see `max_cpus`.

    njobs : int
        The number of jobs to execute. Optional. The number of cores
        will be compared to `njobs`.

    max_cpus : int
        The maximum number of CPUs allowed. If not specified, defaults
        to the available CPUs minus one.

    Raises
    ------
    SetupError
        If `n` is not positive or not convertable to `int`.

    Returns
    -------
    int
        A correct number of cores according to specifications.
    """
    max_cpus = max_cpus or max(cpu_count() - 1, 1)

    if n is None:
        return max_cpus

    try:
        n = int(n)
    except (TypeError, ValueError) as err:
        _msg = f"`n` must be `int` or `int`-convertable `str`: {n!r} given."
        raise SetupError(_msg) from err

    if n < 1:
        _msg = f"`n` is not positive, this is not possible: {n!r}"
        raise SetupError(_msg)

    if njobs:
        ncores = min(n, njobs, max_cpus)
        log.info(
            f"Selected {ncores} cores to process {njobs} jobs, with {max_cpus} "
            "maximum available cores."
            )
        return ncores

    log.debug(f"`njobs` not specified, evaluating initial value {n}...")
    ncores = min(n, max_cpus)
    log.debug(f"Selected {ncores} for a maximum of {max_cpus} CPUs")
    return ncores


def non_negative_int(
        n,
        exception=ValueError,
        emsg="`n` do not satisfies",
        ):
    """
    Transform `n` in int and returns if `compare` evaluates to True.

    Parameters
    ----------
    n : int-convertable
        Something that can be converted to int.

    exception : Exception
        The Exception to raise in case `n` is not a positive integer.

    emsg : str
        The error message to give to `exception`. May accept formatting
        to pass `n`.

    Raises
    ------
    ValueError, TypeError
        If `n` cannot be converted to `int`
    """
    n1 = int(n)
    if n1 >= 0:
        return n1

    # don't change to f-strings, .format has a purpose
    raise exception(emsg.format(n))


def file_exists(
        path,
        exception=ValueError,
        emsg="`path` is not a file or does not exist",
        ):
    """
    Assert if file exist.

    Parameters
    ----------
    path : str or pathlib.Path
        The file path.

    exception : Exception
        The Exception to raise in case `path` is not file or does not
        exist.

    emsg : str
        The error message to give to `exception`. May accept formatting
        to pass `path`.

    Raises
    ------
    Exception
        Any exception that pathlib.Path can raise.
    """
    p = Path(path)

    valid = [p.exists, p.is_file]

    if all(f() for f in valid):
        return p

    # don't change to f-strings, .format has a purpose
    raise exception(emsg.format(str(path)))


def recursive_dict_update(d, u):
    """
    Update dictionary `d` according to `u` recursively.

    https://stackoverflow.com/questions/3232943

    Returns
    -------
    dict
        A new dict object with updated key: values. The original dictionaries
        are not modified.
    """
    def _recurse(d_, u_):
        for k, v in u_.items():
            if isinstance(v, collections.abc.Mapping):
                d_[k] = _recurse(d_.get(k, {}), v)
            else:
                d_[k] = deepcopy(v)  # in case these are also lists

        return d_

    new = deepcopy(d)
    _recurse(new, u)

    return new


def get_number_from_path_stem(path):
    """
    Extract tail number from path.

    Examples
    --------

        >>> get_number_from_path_stem('src/file_1.pdb')
        >>> 1

        >>> get_number_from_path_stem('src/file_3.pdb')
        >>> 3

        >>> get_number_from_path_stem('file_1231.pdb')
        >>> 1231

        >>> get_number_from_path_stem('src/file11')
        >>> 11

        >>> get_number_from_path_stem('src/file_1234_1.pdb')
        >>> 1

    Parameters
    ----------
    path : str or Path obj
        The path to evaluate.

    Returns
    -------
    int
        The tail integer of the path.
    """
    stem = Path(path).stem
    number = re.findall(r'\d+', stem)[-1]
    return int(number)


def get_number_of_digits(num):
    """
    Get the number of digits of a number.

    10 has two digits.
    100 has three digits.
    """
    # also: return len(str(num)) :-)
    return max(ceil(log10(num + 1)), 1)


def get_zerofill_for_modules(modules):
    """
    Get a the prefix zerofill for modules.

    If there are 5 modules, zerofill digits is 1.
    If there are 10 modules, zerofill digits is 1 because counting
    starts at 0 (for topoaa).
    If there are 100 modules, zerofill digits is 2 because counting
    starts at 0 (for topoaa).

    This function is used in combination with `zero_fill`.
    """
    return get_number_of_digits(len(modules) - 1)


def sort_numbered_paths(*paths):
    """
    Sort input paths to tail number.

    If possible, sort criteria is provided by :py:func:`get_number`.
    If paths do not have a numbered tag, sort paths alphabetically.

    Parameters
    ----------
    *inputs : str or pathlib.Path
        Paths to files.

    Returns
    -------
    list
        The sorted pathlist. The original types are not modified. If
        strings are given, strings are returns, if Paths are given
        paths are returned.
    """
    try:
        return sorted(paths, key=get_number_from_path_stem)
    except TypeError as err:
        log.exception(err)
        emsg = (
            "Mind the packing *argument, input should be strings or Paths, "
            "not a list."
            )
        raise TypeError(emsg)
    except IndexError:
        return sorted(paths, key=lambda x: Path(x).stem)


@contextlib.contextmanager
def log_error_and_exit():
    """Exit with exception."""
    try:
        yield
    except Exception as err:
        log.exception(err)
        log.error(err)
        log.error(
            'An error has occurred, see log file. '
            'And contact the developers if needed.'
            )
        log.info(get_goodbye_help())
        sys.exit(1)


def convert_seconds_to_min_sec(seconds):
    """
    Convert seconds to min&sec.

    Examples
    --------
    >>> convert_seconds_to_min_sec(60)
    1 minute
    >>> convert_seconds_to_min_sec(120)
    2 minutes and 0 seconds
    >>> convert_seconds_to_min_sec(179)
    2 minutes and 59 seconds

    Parameters
    ----------
    seconds : int
        The elapsed time in seconds. Seconds are round to integers.

    Returns
    -------
    str
    """
    seconds = int(round(seconds, 0))
    hours = seconds // 3600
    minutes = (seconds) // 60 % 60
    seconds = (seconds) % 60

    if hours:
        return f"{hours}h{minutes}m{seconds}s"

    if minutes:
        s = "" if minutes == 1 else "s"
        return f"{minutes} minute{s} and {seconds} seconds"

    s = "" if seconds == 1 else "s"
    return f"{seconds} seconds"


def extract_keys_recursive(config):
    """Extract keys recursively for the needed modules."""
    for param_name, value in config.items():
        if isinstance(value, collections.abc.Mapping):
            yield from extract_keys_recursive(value)
        else:
            yield param_name
