"""Test gear zero fill."""
import pytest

from haddock.gear.zerofill import (
    _ZeroFill,
    get_number_of_digits,
    get_zerofill_for_modules,
    make_zero_fill,
    zero_fill,
    )


@pytest.mark.parametrize(
    "num,expected",
    [
        (0, 1),
        (1, 1),
        (9, 1),
        (10, 2),
        (22, 2),
        (99, 2),
        (100, 3),
        (335, 3),
        ]
    )
def test_get_number_of_digits(num, expected):
    assert get_number_of_digits(num) == expected


@pytest.mark.parametrize(
    "mods,expected",
    [
        [[f"mod{i}" for i in range(1)], 1],
        [[f"mod{i}" for i in range(5)], 1],
        [[f"mod{i}" for i in range(10)], 1],
        [[f"mod{i}" for i in range(11)], 2],
        [[f"mod{i}" for i in range(99)], 2],
        [[f"mod{i}" for i in range(100)], 2],
        [[f"mod{i}" for i in range(101)], 3],
        [[f"mod{i}" for i in range(1000)], 3],
        [[f"mod{i}" for i in range(1001)], 4],
        ]
    )
def test_get_zerofill_for_modules(mods, expected):
    """
    Test zerofill for modules.

    Here "mod#" strings represent modules. The number prefix means nothing.
    """
    assert get_zerofill_for_modules(mods) == expected


@pytest.mark.parametrize(
    "num,dig,exp",
    [
        (1, 1, "1"),
        (10, 1, "10"),
        (0, 2, "00"),
        (3, 2, "03"),
        (50, 3, "050"),
        ],
    )
def test_make_zero_fill(num, dig, exp):
    result = make_zero_fill(num, dig)
    assert result == exp


@pytest.mark.parametrize(
    "modules, exp_num, exp_folder",
    [
        ([f"m{i}" for i in range(10)], 1, "5_module"),
        ([f"m{i}" for i in range(20)], 2, "05_module"),
        ]
    )
def test_zerofill_class(modules, exp_num, exp_folder):
    zf = _ZeroFill()
    zf.read(modules)
    assert zf.zfnum == exp_num
    fresult = zf.fill("module", 5)
    assert fresult == exp_folder


def test_zerofill_singleton():
    assert isinstance(zero_fill, _ZeroFill)
