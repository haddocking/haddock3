from pathlib import Path

import pytest

from haddock.clis.cli import ap


recipe = Path(Path(__file__).resolve().parent, 'recipe.toml')


def test_ap_recipe_does_not_exist():
    with pytest.raises(SystemExit) as exit:
        ap.parse_args('does_not_exit.toml'.split())
    assert exit.type == SystemExit
    assert exit.value.code == 2


def test_ap_recipe_exists():
    cmd = ap.parse_args(str(recipe).split())
    with open(cmd.recipe) as fin:
        fin.readlines()


def test_ap_setup_true():
    cmd = ap.parse_args(f'{recipe} --setup'.split())
    assert cmd.setup_only == True


def test_ap_setup_false():
    cmd = ap.parse_args(str(recipe).split())
    assert cmd.setup_only == False


def test_ap_version():
    with pytest.raises(SystemExit) as exit:
        ap.parse_args('-v'.split())
    assert exit.type == SystemExit
    assert exit.value.code == 0


@pytest.mark.parametrize(
    'level',
    ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
    )
def test_ap_log_level(level):
    cmd = ap.parse_args(f'{recipe} --log-level {level}'.split())
    assert cmd.log_level == level


def test_ap_log_level_error():
    with pytest.raises(SystemExit) as exit:
        ap.parse_args(f'{recipe} --log-level BAD'.split())
    assert exit.type == SystemExit
    assert exit.value.code == 2
