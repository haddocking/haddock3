"""General utilities."""
import collections.abc
import re
import shutil
import subprocess
from copy import deepcopy
from functools import partial
from os import cpu_count
from pathlib import Path

from haddock import log
from haddock.core.exceptions import SetupError


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
    Update dictionary recursively.

    https://stackoverflow.com/questions/3232943
    """
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = recursive_dict_update(d.get(k, {}), v)
        else:
            d[k] = v
    return d


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
