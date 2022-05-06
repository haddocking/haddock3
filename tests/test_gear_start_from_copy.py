"""Test start from copy gear."""
import argparse
from pathlib import Path

import pytest

from haddock.gear.start_from_copy import (
    add_start_from_copy,
    )


@pytest.mark.parametrize(
    'path,expected',
    [
        ('path', Path('path')),
        ('.', Path('.')),
        ]
    )
def test_add_argument_copy_args(path, expected):
    ap = argparse.ArgumentParser()
    add_start_from_copy(ap)
    cmd = ap.parse_args(f'--start-from-copy {path}'.split())
    assert cmd.start_from_copy == expected


def test_add_argument_copy_none():
    ap = argparse.ArgumentParser()
    add_start_from_copy(ap)
    cmd = ap.parse_args([])
    assert cmd.start_from_copy is None
