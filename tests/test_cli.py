"""Test main haddock3 client."""
from pathlib import Path

import pytest

from haddock.clis import cli

from . import configs_data


@pytest.fixture(name="workflow")
def fixture_workflow():
    yield Path(configs_data, "recipe.cfg")


def test_cli_has_maincli():
    """
    Test maincli func in CLI.

    maincli is used in setup.py.
    """
    assert cli.maincli


def test_ap_recipe_does_not_exist():
    """Test raise error if workflow does not exist."""
    with pytest.raises(SystemExit) as exit:
        cli.ap.parse_args('does_not_exit.cfg'.split())
    assert exit.type == SystemExit
    assert exit.value.code == 2


def test_ap_workflow_exists(workflow):
    """Test reading workflows."""
    cmd = cli.ap.parse_args(str(workflow).split())
    with open(cmd.workflow) as fin:
        fin.readlines()


def test_ap_setup_true(workflow):
    """Test --setup flag."""
    cmd = cli.ap.parse_args(f'{workflow} --setup'.split())
    assert cmd.setup_only is True


def test_ap_setup_false(workflow):
    """Test setup only default."""
    cmd = cli.ap.parse_args(str(workflow).split())
    assert cmd.setup_only is False


def test_ap_version():
    """Test -v version flag."""
    with pytest.raises(SystemExit) as exit:
        cli.ap.parse_args('-v'.split())
    assert exit.type == SystemExit
    assert exit.value.code == 0


@pytest.mark.parametrize(
    'level',
    ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
    )
def test_ap_log_level(workflow, level):
    """Test --log-level correct."""
    cmd = cli.ap.parse_args(f'{workflow} --log-level {level}'.split())
    assert cmd.log_level == level


def test_ap_log_level_error(workflow):
    """Test --log-level error with bad input."""
    with pytest.raises(SystemExit) as exit:
        cli.ap.parse_args(f'{workflow} --log-level BAD'.split())
    assert exit.type == SystemExit
    assert exit.value.code == 2
