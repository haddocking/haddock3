"""Lib I/O."""
import contextlib
import os
from pathlib import Path

import yaml


def read_from_yaml(yaml_file):
    """Read a configuration from a yaml file."""
    with open(cfg, 'r') as fin:
        ycfg = yaml.safe_load(fin)
    return ycfg


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
