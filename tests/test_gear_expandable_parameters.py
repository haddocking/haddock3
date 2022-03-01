"""Test expandable parameter modules."""
import pytest

from haddock.gear.expandable_parameters import (

    belongs_to_multiple_index,
    belongs_to_single_index,
    extract_multiple_index_params,
    extract_single_index_params,
    rejoin_parts_multiple_index,
    rejoin_parts_single_index,
    make_param_name_single_index,
    make_param_name_multiple_index,

    )


@pytest.mark.parametrize(
    "param,expected",
    (
        ("param_tag_1", True),
        ("param_tag_else_1", True),
        ("param_tag_else2_1", True),
        ("param", False),
        ("param_1", False),
        ("param_1_2", False),
        ("param_other", False),
        ("param_other_1_2", False),
        )
    )
def test_belongs_to_single_index(param, expected):
    """Test if param belongs to single index."""
    result = belongs_to_single_index(param.split("_"))
    assert result == expected


@pytest.mark.parametrize(
    "param,expected",
    (
        ("param_tag_1", False),
        ("param_tag_else_1", False),
        ("param_tag_else2_1", False),
        ("param", False),
        ("param_1", False),
        ("param_1_2", False),
        ("param_other", False),
        ("param_other_1_2", True),
        ("param_other_else_2_1", True),
        ("param_other2_else_2_1", True),
        )
    )
def test_belongs_to_multiple_index(param, expected):
    """Test if param belongs to single index."""
    result = belongs_to_multiple_index(param.split("_"))
    assert result == expected


@pytest.mark.parametrize(
    "parts,expected",
    (
        (("param", "other", "1"), "other"),
        (("param", "other", "else", "1"), "other_else"),
        (("param", "other", "else2", "2"), "other_else2"),
        (("param", "other", "else2", "3", "2"), "other_else2_3"),
        )
    )
def test_rejoin_parts_single_index(parts, expected):
    """Test join part of parameter names in single index group."""
    result = rejoin_parts_single_index(parts)
    assert result == expected


@pytest.mark.parametrize(
    "parts,expected",
    (
        (("param", "other", "1", "1"), "other"),
        (("param", "other", "else", "1", "1"), "other_else"),
        (("param", "other", "else2", "2", "2"), "other_else2"),
        (("param", "other", "else2", "3", "2"), "other_else2"),
        )
    )
def test_rejoin_parts_multiple_index(parts, expected):
    """Test join part of parameter names in single index group."""
    result = rejoin_parts_multiple_index(parts)
    assert result == expected


@pytest.mark.parametrize(
    "user_config,param_name,group_idx,expected",
    (
        (
            {"param_some_1", "param_other_1", "param2", "ppp_ooo"},
            "param", "1",
            {"param_some_1", "param_other_1"},
            ),
        (
            {"param_some_1", "param_other_1", "param_some_2", "param_other_2"},
            "param", "2",
            {"param_some_2", "param_other_2"},
            ),
        (
            {"param_some_1", "param_other_1", "par_some_2_2", "par_other_2_2"},
            "param", "1",
            {"param_some_1", "param_other_1"},
            ),
        (
            {"param_some_1", "param_other_1", "par_some_2_2", "par_other_2_2"},
            "par", "2",
            {"par_some_2_2", "par_other_2_2"},
            ),
        (
            {"param_some_1", "param_other_1", "par_some_2_2", "par_other_2_2"},
            "por", "2",
            set(),
            ),
        )
    )
def test_extract_single_index_params(
        user_config,
        param_name,
        group_idx,
        expected,
        ):
    """Test extract single index params from config."""
    result = extract_single_index_params(user_config, param_name, group_idx)
    assert result == expected


@pytest.mark.parametrize(
    "user_config,param_name,group_idx,expected",
    (
        (
            {"param_some_1", "param_other_1", "param2", "ppp_ooo"},
            "param", "1",
            set(),
            ),
        (
            {"param_some_1", "param_other_1", "param_some_2", "param_other_2"},
            "param", "2",
            set(),
            ),
        (
            {"param_some_1", "param_other_1", "par_some_2_2", "par_other_2_2"},
            "par2", "2",
            {"par_some_2_2", "par_other_2_2"},
            ),
        (
            {"param_some_1", "param_other_1", "par_some_2_2", "par_other_2_2"},
            "par", "1",
            set(),
            ),
        (
            {"param_some_1", "param_other_1", "par_some_2_2", "par_other_2_2"},
            "por", "2",
            set(),
            ),
        )
    )
def test_extract_multiple_index_params(
        user_config,
        param_name,
        group_idx,
        expected,
        ):
    """Test extract single index params from config."""
    result = extract_multiple_index_params(user_config, param_name, group_idx)
    assert result == expected


@pytest.mark.parametrize(
    "parts,expected",
    (
        (("param", "other", "1"), ("param", "1")),
        (("param", "other", "else", "1"), ("param", "1")),
        (("param", "other", "else", "2"), ("param", "2")),
        )
    )
def test_make_param_name_single(parts, expected):
    """Test make param namel single."""
    result = make_param_name_single_index(parts)
    assert result == expected


@pytest.mark.parametrize(
    "parts,expected",
    (
        # this is senseless for the logic but it works
        (("param", "other", "1"), ("paramother", "1")),
        # logic examples
        (("param", "other", "else", "3", "1"), ("param3", "1")),
        (("param", "other", "else", "2", "2"), ("param2", "2")),
        )
    )
def test_make_param_name_multiple(parts, expected):
    """Test make param namel single."""
    result = make_param_name_multiple_index(parts)
    assert result == expected
