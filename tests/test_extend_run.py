"""Test extend run gear."""
import argparse
from pathlib import Path

import pytest

from haddock.gear.extend_run import add_extend_run


@pytest.mark.parametrize(
    'path,expected',
    [
        ('path', Path('path')),
        ('.', Path('.')),
        ]
    )
def test_add_argument_copy_args(path, expected):
    ap = argparse.ArgumentParser()
    add_extend_run(ap)
    cmd = ap.parse_args(f'--extend-run {path}'.split())
    assert cmd.extend_run == expected


def test_add_argument_copy_none():
    ap = argparse.ArgumentParser()
    add_extend_run(ap)
    cmd = ap.parse_args([])
    assert cmd.extend_run is None
