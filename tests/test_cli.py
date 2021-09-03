import pytest

from haddock.clis.cli import ap


def test_ap_recipe_does_not_exist():
    with pytest.raises(SystemExit) as exit:
        ap.parse_args('does_not_exit.toml'.split())
    assert exit.type == SystemExit
    assert exit.value.code == 2

def test_ap_recipe_exists():
    cmd = ap.parse_args('tests/recipe.toml'.split())
    cmd.recipe.close()


def test_ap_3():
    cmd = ap.parse_args('tests/recipe.toml --setup'.split())
    assert cmd.setup_only == True

def test_ap_4():
    with pytest.raises(SystemExit) as exit:
        ap.parse_args('-v'.split())
    assert exit.type == SystemExit
    assert exit.value.code == 0
