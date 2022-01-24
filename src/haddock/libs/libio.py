"""Lib I/O."""
import contextlib
import os
from pathlib import Path


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
        Path(file_).write_text(content)

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
