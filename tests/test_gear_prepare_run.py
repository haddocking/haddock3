"""Test prepare run module."""
import pytest

from haddock.gear.prepare_run import (
    get_blocks_single_index,
    read_single_index_blocks,
    )


@pytest.mark.parametrize(
    "inp, expected",
    [
        (
            {
                "par1": 1,
                "par_1": 1,
                "this_is_not_a_group_1": None,
                "group_something_1": None,
                "group_else_1": None,
                "par2": None,
                "group2_else_other_1": None,
                "group2_else_many_1": None,
                },
            {
                ("group", "1"): {"something", "else"},
                ("group2", "1"): {"else_other", "else_many"},
                },
            )
        ]
    )
def test_get_blocks(inp, expected):
    """Test get blocks."""
    result = get_blocks_single_index(inp)
    assert result == expected


@pytest.mark.parametrize(
    "default, user, expected",
    [
        (
            {
                ("group", "1"): {"something", "else"},
                ("group2", "1"): {"else_other", "else_many"},
                },
            {
                "par1": 1,
                "par_1": 1,
                "group_something_1": None,
                "group_else_1": None,
                "group_something_2": None,
                "group_else_2": None,
                "group2_else_other_1": None,
                "group2_else_many_1": None,
                "group2_else_other_2": None,
                "group2_else_many_2": None,
                "group2_else_other_3": None,
                "group2_else_many_3": None,
                },
            {
                "group_something_1",
                "group_else_1",
                "group_something_2",
                "group_else_2",
                "group2_else_other_1",
                "group2_else_many_1",
                "group2_else_other_2",
                "group2_else_many_2",
                "group2_else_other_3",
                "group2_else_many_3",
                }
            )
        ]
    )
def test_read_single_index_blocks(default, user, expected):
    """Test read_blocks."""
    result = read_single_index_blocks(default, user)
    assert result == expected
