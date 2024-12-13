"""General utilities."""
import collections.abc
import contextlib
import os
import re
import shutil
import subprocess
import sys
from copy import deepcopy
from functools import partial
from pathlib import Path

from haddock import EmptyPath, log
from haddock.core.exceptions import SetupError
from haddock.core.typing import (
    PT,
    Any,
    AnyT,
    Callable,
    Container,
    FilePath,
    FilePathT,
    Generator,
    Iterable,
    Optional,
    ParamDict,
    ParamMap,
    ParamMapT,
    Union,
)
from haddock.gear.greetings import get_goodbye_help


check_subprocess = partial(
    subprocess.run,
    shell=True,
    check=True,
    stdout=subprocess.DEVNULL,
)


def get_result_or_same_in_list(
    function: Callable[[PT], AnyT], value: PT
) -> Union[AnyT, list[PT]]:
    """
    Return the result if True or the value within a list.

    Applies `function` to `value` and returns its result if it evaluates
    to True. Otherwise, return the value within a list.

    `function` should receive a single argument, the `value`.
    """
    result = function(value)
    return result if result else [value]


def make_list_if_string(item: Union[str, list[str]]) -> list[str]:
    """Put `item` into a list."""
    if isinstance(item, str):
        return [item]
    return item


def transform_to_list(
    item: Union[Iterable[AnyT], AnyT]
) -> Union[list[AnyT], tuple[AnyT, ...]]:
    """
    Put `item` into a list if not a list already.

    If it is set, transforms the set into a list.

    If it is a dict, returns a list of the keys.

    If it is tuple, returns the tuple.

    If a list, returns the same.

    Everything else returns `item` inside a one element list.
    """
    if isinstance(item, (set, dict)):
        return list(item)

    if isinstance(item, (list, tuple)):
        return item

    return [item]


def copy_files_to_dir(paths: Iterable[FilePath], directory: FilePath) -> None:
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


def remove_folder(folder: FilePath) -> None:
    """
    Remove a folder if it exists.

    Parameters
    ----------
    folder : str or Path
        Path to folder to remove.
    """
    if Path(folder).exists():
        log.warning(f"{folder} exists and it will be REMOVED!")
        shutil.rmtree(folder)


def remove_dict_keys(d: ParamMap, keys: Container[str]) -> ParamDict:
    """
    Remove `keys` from dictionary (`d`).

    Return
    ------
    dict
        A copy of `d` dictionary without the `keys`.
    """
    return {k: deepcopy(v) for k, v in d.items() if k not in keys}


def cpu_count() -> int:
    """Count number of available CPU for the process.
    
    User suggestion, by https://github.com/EricDeveaud

    Note: pid 0 == current process

    FIXME: from python3.13, better to use os.process_cpu_count()

    Returns
    -------
    process_ncores : int
        Number of cores allocated to the process pid.
    """
    def _cpu_count() -> int:
        """Detect number of cores available.
        
        Returns
        -------
        ncores : int
            Maximum number of cores that can be used.
        """
        try:
            # Try to obtain the number of cpu allocated to the process
            _process_ncores = os.sched_getaffinity(0)
            if isinstance(_process_ncores, set):
                ncores = len(_process_ncores)
            else:
                ncores = int(_process_ncores)
        except AttributeError:
            # If unsucessful, return the number cores detected in the machine
            # Note: this can happen on MacOS, where the os.sched_getaffinity 
            # may not be defined/exist.
            ncores = int(os.cpu_count())
        return ncores

    try:
        process_ncores = int(os.process_cpu_count())
    except AttributeError:
        process_ncores = _cpu_count()
    return process_ncores


def parse_ncores(
    n: Optional[Union[int, str]] = None,
    njobs: Optional[int] = None,
    max_cpus: Optional[bool] = None,
) -> int:
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
    if max_cpus is None or max_cpus is False:
        max_cpus = max(cpu_count() - 1, 1)  # type: ignore
    if max_cpus is True:
        max_cpus = cpu_count()  # type: ignore
    elif not isinstance(max_cpus, int):
        raise TypeError(f"`max_cpus` not of valid type: {type(max_cpus)}")

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

    if njobs is not None:
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
    n: Any,
    exception: type[Exception] = ValueError,
    emsg: str = "`n` do not satisfies",
) -> int:
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
    try:
        n1 = int(n)
        if n1 >= 0:
            return n1
    except Exception as e:
        raise e
    else:
        # don't change to f-strings, .format has a purpose
        raise exception(emsg.format(n))


def recursive_dict_update(d: ParamMapT, u: ParamMap) -> ParamMapT:
    """
    Update dictionary `d` according to `u` recursively.

    https://stackoverflow.com/questions/3232943

    Returns
    -------
    dict
        A new dict object with updated key: values. The original dictionaries
        are not modified.
    """

    def _recurse(d_: ParamMapT, u_: ParamMap) -> ParamMapT:
        for k, v in u_.items():
            if isinstance(v, collections.abc.Mapping):
                d_[k] = _recurse(d_.get(k, {}), v)
            else:
                d_[k] = deepcopy(v)  # in case these are also lists

        return d_

    new = deepcopy(d)
    _recurse(new, u)

    return new


def get_number_from_path_stem(path: FilePath) -> int:
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
    number = re.findall(r"\d+", stem)[-1]
    return int(number)


def sort_numbered_paths(*paths: FilePathT) -> list[FilePathT]:
    """
    Sort input paths to tail number.

    If possible, sort criteria is provided by
    :py:func:`get_number_from_path_stem`.
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
def log_error_and_exit() -> Generator[None, None, None]:
    """Exit with exception."""
    try:
        yield
    except Exception as err:
        log.exception(err)
        log.error(err)
        log.error(
            "An error has occurred, see log file. "
            "And contact the developers if needed."
        )
        log.info(get_goodbye_help())
        sys.exit(1)


def extract_keys_recursive(config: ParamMap) -> Generator[str, None, None]:
    """Extract keys recursively for the needed modules."""
    for param_name, value in config.items():
        if isinstance(value, collections.abc.Mapping):
            yield from extract_keys_recursive(value)
        else:
            yield param_name


def recursive_convert_paths_to_strings(params: ParamMapT) -> ParamMapT:
    """
    Convert paths to strings recursively over a dictionary.

    Parameters
    ----------
    params : dictionary

    Returns
    -------
    dictionary
        A copy of the original dictionary with paths converted to strings.
    """
    params = deepcopy(params)
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
