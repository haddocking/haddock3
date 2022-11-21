"""Test lib client."""
import argparse
from pathlib import Path

import pytest

from haddock.libs import libcli


@pytest.mark.parametrize(
    'cmd,expected',
    [
        ('-odir my_folder', Path("my_folder")),
        ('--output-directory my_folder/other', Path("my_folder", "other")),
        ]
    )
def test_output_dir(cmd, expected):
    """Test output directory."""
    parser = argparse.ArgumentParser()
    libcli.add_output_dir_arg(parser)
    v = vars(parser.parse_args(cmd.split()))
    assert v['output_directory'] == expected


def test_output_dir_default():
    """Test output directory."""
    parser = argparse.ArgumentParser()
    libcli.add_output_dir_arg(parser)
    v = vars(parser.parse_args([]))
    assert v['output_directory'] == Path.cwd()


@pytest.mark.parametrize(
    'cmd,expected',
    [
        ('-n 4', 4),
        ('--ncores 2', 2),
        ('', 1),
        ('-n', None),
        ('--ncores', None),
        ],
    )
def test_ncores(cmd, expected):
    """Test output directory."""
    parser = argparse.ArgumentParser()
    libcli.add_ncores_arg(parser)
    v = vars(parser.parse_args(cmd.split()))
    assert v['ncores'] == expected
