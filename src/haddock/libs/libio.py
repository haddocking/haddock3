"""Lib I/O."""
import contextlib
import os
from pathlib import Path

import yaml


def read_from_yaml(yaml_file):
    """
    Read a YAML file to a dictionary.

    Used internally to read HADDOCK3's default configuration files.

    Parameters
    ----------
    yaml_file : str or Path
        Path to the YAML file.

    Returns
    -------
    dict
        Always returns a dictionary.
        Returns empty dictionary if yaml_file is empty.
    """
    with open(yaml_file, 'r') as fin:
        ycfg = yaml.safe_load(fin)

    # ycfg is None if yaml_file is empty
    # returns an empty dictionary to comply with HADDOCK workflow
    if ycfg is None:
        return {}

    assert isinstance(ycfg, dict), type(ycfg)
    return ycfg


def open_files_to_lines(*files):
    """
    Open files to lines.

    New-lines are stripped.

    Returns
    -------
    list of lists of strings
        The lines of the files.
        Input order is maintained.
    """
    f_paths = map(Path, files)
    return [f.read_text().split(os.linesep) for f in f_paths]


def save_lines_to_files(files, lines):
    """
    Save a list of list of lines to files.

    The first list of strings in `lines` will be saved in the first file
    of `files`, and so on.

    Lines are saved using `pathlib.Path.write_text` function.

    Parameters
    ----------
    files : list
        The list of file names to save.

    lines : list of lists of str
        A list containing lists of lines that are the file contents.
        Must be synched with `files`.
    """
    for file_, content in zip(files, lines):
        Path(file_).write_text(os.linesep.join(content) + os.linesep)

    return


def add_suffix_to_files(files, suffix):
    """
    Add a suffix to file paths.

    Yields
    ------
    pathlib.Path objects
        Exhausts when files exhaust.
    """
    for file_ in files:
        p = Path(file_)
        folder = p.parent
        psuffix = p.suffix
        name = p.stem + suffix + psuffix
        path = Path(folder, name)
        yield path


# thanks to @brianjimenez
@contextlib.contextmanager
def working_directory(path):
    """Change working directory and returns to previous on exit."""
    prev_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)
