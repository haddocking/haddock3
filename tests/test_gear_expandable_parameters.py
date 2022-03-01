"""Test expandable parameter modules."""
import importlib

import pytest

from haddock.core.exceptions import ConfigurationError
from haddock.gear.expandable_parameters import (
    belongs_to_multiple_index,
    belongs_to_single_index,
    extract_multiple_index_params,
    extract_single_index_params,
    get_multiple_index_groups,
    get_single_index_groups,
    make_param_name_multiple_index,
    make_param_name_single_index,
    read_multiple_idx_groups_user_config,
    read_single_idx_groups_user_config,
    rejoin_parts_multiple_index,
    rejoin_parts_single_index,
    type_simplest_ep,
    )
from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.libs.libutil import extract_keys_recursive
from haddock.modules import modules_category


def test_type_simplest_ep():
    """Test all keys are modules."""
    for key in type_simplest_ep:
        assert key in modules_category


@pytest.mark.parametrize(
    "mod_cat",
    [(k, v) for k, v in modules_category.items() if k in type_simplest_ep],
    )
def test_type_simplest_keys(mod_cat):
    """Test all parameters are defined in defaults."""
    module_name, category = mod_cat
    mod = ".".join(['haddock', 'modules', category, module_name])
    module = importlib.import_module(mod)

    config = read_from_yaml_config(module.DEFAULT_CONFIG)
    params = set(extract_keys_recursive(config))
    module_expandable = type_simplest_ep[module_name]

    # using for-loop instead of subset to give better error messages
    for _param in module_expandable:
        assert f"{_param}_1" in params, f"{_param} not in {module_name}"


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
            set(),
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


@pytest.mark.parametrize(
    "config,expected",
    # here only the keys are represented
    [
        [
            {
                "param_tag_1",
                "param_teg_1",
                "param_tig_1",
                "param_tog_1",
                "param_tug_1",
                "param_tag_2",
                "param_teg_2",
                "param_tig_2",
                "param_tog_2",
                "param_tug_2",
                "param_1",
                "other_sta_1",
                "other_end_1",
                "param_newtag_1_1",
                "param_newteg_1_1",
                "param_newtig_1_1",
                },
            {
                ("param", "1"): {"tag", "teg", "tig", "tog", "tug"},
                ("other", "1"): {"sta", "end"},
                },
            ],
        ],
    )
def test_get_single_index_group(config, expected):
    """Test get single index group."""
    result = get_single_index_groups(config)
    assert result == expected


@pytest.mark.parametrize(
    "config,expected",
    # here only the keys are represented
    [
        [
            {
                "param_tig_1",
                "param_tog_1",
                "param_tug_1",
                "param_tag_2",
                "param_teg_2",
                "param_tig_2",
                "param_1",
                "other_sta_1",
                "other_end_1",
                "param_newtag_1_1",
                "param_newteg_1_1",
                "param_newtig_1_1",
                "other_new_2_1",
                "other_mew_2_1",
                },
            {
                ("param1", "1"): {"newtag", "newteg", "newtig"},
                ("other2", "1"): {"new", "mew"},
                },
            ],
        ],
    )
def test_get_multiple_index_group(config, expected):
    """Test get single index group."""
    result = get_multiple_index_groups(config)
    assert result == expected


@pytest.mark.parametrize(
    "user_config,default_groups,expected",
    [
        (
            {"param_other_1", "param_else_1", "param_1",
             "param_other_2", "param_else_2", "param_other_1_1"},
            {("param", "1"): {"other", "else"}},
            {"param_other_1", "param_else_1", "param_other_2", "param_else_2"},
            ),
        ]
    )
def test_read_single_idx_groups_user_config(
        user_config,
        default_groups,
        expected,
        ):
    """Test read single index groups in user config."""
    result = read_single_idx_groups_user_config(user_config, default_groups)
    assert result == expected


@pytest.mark.parametrize(
    "user_config,default_groups,expected",
    [
        (
            {"param_other_1", "param_else_1",
             "param_other_2", "param_else_2",
             "param_1",
             "param_sam_1_1", "param_frodo_1_1",
             "param_sam_1_2", "param_frodo_1_2",
             "param_sam_2_1", "param_frodo_2_1",
             },
            {("param1", "1"): {"sam", "frodo"},
             ("param2", "1"): {"sam", "frodo"}},
            {"param_frodo_1_1", "param_sam_1_1",
             "param_frodo_1_2", "param_sam_1_2",
             "param_sam_2_1", "param_frodo_2_1",
             },
            ),
        ]
    )
def test_read_multiple_idx_groups_user_config(
        user_config,
        default_groups,
        expected,
        ):
    """Test read single index groups in user config."""
    result = read_multiple_idx_groups_user_config(user_config, default_groups)
    assert result == expected


@pytest.mark.parametrize(
    "user_config,default_groups",
    [
        (
            {"param_other_1", "param_else_1", "param_1",
             "param_other_2", "param_else_2", "param_other_1_1",
             "pirim_some_1", "pirim_other_1"},
            {("param", "1"): {"other", "else"}},
            ),
        ]
    )
def test_read_single_idx_groups_user_config_no_group_error(
        user_config,
        default_groups,
        ):
    """Test read single index groups in user config."""
    with pytest.raises(ConfigurationError):
        read_single_idx_groups_user_config(user_config, default_groups)


@pytest.mark.parametrize(
    "user_config,default_groups",
    [
        (
            {"param_other_1", "param_else_1", "param_bad_1", "param_1",
             "param_other_2", "param_else_2", "param_other_1_1"},
            {("param", "1"): {"other", "else"}},
            ),
        ]
    )
def test_read_single_idx_groups_user_config_unexpeceted_error(
        user_config,
        default_groups,
        ):
    """Test read single index groups in user config."""
    with pytest.raises(ConfigurationError):
        read_single_idx_groups_user_config(user_config, default_groups)


@pytest.mark.parametrize(
    "user_config,default_groups",
    [
        (
            {"param_other_1", "param_1",
             "param_other_2", "param_else_2", "param_other_1_1"},
            {("param", "1"): {"other", "else"}},
            ),
        ]
    )
def test_read_single_idx_groups_user_config_less_error(
        user_config,
        default_groups,
        ):
    """Test read single index groups in user config."""
    with pytest.raises(ConfigurationError):
        read_single_idx_groups_user_config(user_config, default_groups)


@pytest.mark.parametrize(
    "user_config,default_groups",
    [
        (
            {"param_other_1", "param_else_1",
             "param_other_2", "param_else_2",
             "param_1",
             "param_sam_1_1", "param_frodo_1_1",
             "param_sam_1_2", "param_frodo_1_2",
             "param_sam_2_1", "param_frodo_2_1",
             "pirim_sta_1_1", "pirim_end_1_1",
             },
            {("param1", "1"): {"sam", "frodo"},
             ("param2", "1"): {"sam", "frodo"}},
            ),
        ]
    )
def test_read_multiple_idx_groups_user_config_no_group(
        user_config,
        default_groups,
        ):
    """Test read single index groups in user config."""
    with pytest.raises(ConfigurationError):
        read_multiple_idx_groups_user_config(user_config, default_groups)


@pytest.mark.parametrize(
    "user_config,default_groups",
    [
        (
            {"param_other_1", "param_else_1",
             "param_other_2", "param_else_2",
             "param_1",
             "param_sam_1_1", "param_frodo_1_1",
             "param_baggins_1_1",
             "param_sam_1_2", "param_frodo_1_2",
             "param_sam_2_1", "param_frodo_2_1",
             },
            {("param1", "1"): {"sam", "frodo"},
             ("param2", "1"): {"sam", "frodo"}},
            ),
        ]
    )
def test_read_multiple_idx_groups_user_config_unexpected(
        user_config,
        default_groups,
        ):
    """Test read single index groups in user config."""
    with pytest.raises(ConfigurationError):
        read_multiple_idx_groups_user_config(user_config, default_groups)


@pytest.mark.parametrize(
    "user_config,default_groups",
    [
        (
            {"param_other_1", "param_else_1",
             "param_other_2", "param_else_2",
             "param_1",
             "param_sam_1_1",
             "param_sam_1_2", "param_frodo_1_2",
             "param_sam_2_1", "param_frodo_2_1",
             },
            {("param1", "1"): {"sam", "frodo"},
             ("param2", "1"): {"sam", "frodo"}},
            ),
        ]
    )
def test_read_multiple_idx_groups_user_config_less(
        user_config,
        default_groups,
        ):
    """Test read single index groups in user config."""
    with pytest.raises(ConfigurationError):
        read_multiple_idx_groups_user_config(user_config, default_groups)
