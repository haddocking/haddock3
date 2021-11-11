"""Test gear.restart_run."""
import argparse

import pytest

from haddock.gear import restart_run


def test_has_help():
    """Assert module has _help_cli variable."""
    assert restart_run._help_cli


@pytest.mark.parametrize(
    'i,expected',
    [
        ('0', 0),
        ('1', 1),
        ('57', 57),
        (100, 100),
        ]
    )
def test_non_neg_arg(i, expected):
    """Test non negative arg type."""
    r = restart_run._arg_non_neg_int(i)
    assert r == expected


@pytest.mark.parametrize(
    'i,expected',
    [
        ('0', 0),
        ('1', 1),
        ('57', 57),
        (100, 100),
        ]
    )
def test_restart_cli(i, expected):
    """Test non negative arg type."""
    ap = argparse.ArgumentParser()
    restart_run.add_restart_arg(ap)
    cmd = ap.parse_args(f'--restart {i}'.split())
    assert cmd.restart == expected


@pytest.mark.parametrize(
    'n',
    (-1, -10, '-1230', -50000),
    )
def test_arg_non_neg_error(n):
    """Test non-negative int."""
    with pytest.raises(argparse.ArgumentTypeError):
        restart_run._arg_non_neg_int(n)


@pytest.mark.parametrize(
    'n',
    (-1, -10, '-1230', -50000),
    )
def test_restart_cli_error(n):
    """Test --restart number with error."""
    ap = argparse.ArgumentParser()
    restart_run.add_restart_arg(ap)
    with pytest.raises(SystemExit) as exit:
        ap.parse_args(f'--restart {n}'.split())
    assert exit.type == SystemExit
    assert exit.value.code == 2
