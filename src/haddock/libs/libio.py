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
