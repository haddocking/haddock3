"""Lib I/O."""
import contextlib
import os
from pathlib import Path


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


