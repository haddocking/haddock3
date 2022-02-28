"""
Expandable parameter module.

This module contains all logic reading and preparing expandable parameters.

Expandable parameters are parameters that follow a specific name structure.
They can be repeated multiple times (in groups) in the user configuration file
while being defined only once in the `defaults.yaml` files. For example, in the
`refinement.flexref` module there are the parameters:

- `fle_sta_1_1`
- `fle_end_1_1`

For these cases, you can also define in the user configuration file the
following:

- `fle_sta_1_2`
- `fle_end_1_2`

The above parameter will be accepted and used in the CNS modules.

Currently, there are three types of expandable parameters:

1 - Parameters following the name structure `name_something_1`,
`name_something_else_1`. These two parameters belong to the parameter
name `name` and group `1`. The two parameters have the name `something`
and `something_else`. Therefore, you could define also:

- `name_something_2`
- `name_something_else_2`

Numbers are allowed as long as in combination with word characters, for
example: `name_something1_2`, but not `name_something_1_1`. The latter
is of the second type.

2 - Parameters following the name structure `name_something_X_Y`. In
these cases, the parameter name is `nameX`, where `X` is an integer. `Y`
is the number of the group, and `something` is the name of the
parameters. Taking above example, the parameter name is `fle1`; `sta`
and `end` are the parameters of the group, and the group is either `1`
or `2`.

3 - The simplest parameters in the form of `param_1`, where you could
define `param_2`, `param_3`, etc. However, in these cases, the
expandable rules apply only to the parameters manually defined in this
file.
"""
from functools import partial

from haddock.core.exceptions import ConfigurationError


# errors messages for single indexed groups
_emsg_no_group_single = \
    "The parameter block '{}_*_{}' is not a valid expandable parameter."
_emsg_unexpected_params_single = \
    "These parameters do not belong to the block '{}_*_{}': {!r}."

# errors messages for multiple indexed groups
_emsg_no_group_multiple = \
    "The parameter block '{}_*_{}_*' is not a valid expandable parameter."
_emsg_unexpected_params_multiple = \
    "These parameters do not belong to the block '{}_*_{}_*': {!r}."

# common error messages
_emsg_num_differ = (
    "The parameter block {!r} expects "
    "{} parameters, but {} are present "
    "in the user configuration file."
    )


def _get_groups(
        config,
        belongs_to_group,
        make_param_name,
        rejoin_parts,
        minimum=2,
        ):
    """
    Identify expandable parameter groups - main engine function.

    Parameters
    ----------
    config : dict
        The configuration file, where keys are the parameter names.

    belongs_to_group : callable
        A callable receiving the parameter name after str.split("_")
        that returns a boolean identifying if the parameter belongs a
        given expadanble group.

    make_param_name : callable
        Recreates the parameter preffix from the parameter parts.

    rejoin_parts : callable
        Recreates the parameter name from the parts. For example:
        `something_else` in `param_something_else_1`.

    minimum : int
        Consider only the groups with at least `minimum` parameters.
    """
    splitted = (parameter.split("_") for parameter in config)
    parts = (_parts for _parts in splitted if belongs_to_group(_parts))

    # identify groups in config
    groups = {}
    for p in parts:
        group_identity = make_param_name(p)  # ("param", "1")
        new = groups.setdefault(group_identity, {})

        # counts the number of parameters within the block
        new.setdefault("counts", 0)
        new["counts"] += 1

        # registers the parameters belonging to the group
        new.setdefault("mid", set())
        new["mid"].add(rejoin_parts(p))

    # creates the final dictionary
    final_groups = {
        k: v["mid"]
        for k, v in groups.items()
        if v["counts"] >= minimum
        }

    return final_groups


def _read_groups_in_user_config(
        user_config,
        default_groups,
        get_user_groups,
        extract_params,
        _emsg_no_group="Parameter block is not a valid expandable group",
        _emsg_unexpected_params="Unexpected params for group",
        ):
    """
    Read groups in user config.

    This is an engine-like function that serves the different expandable
    parameter rules. This function should be parametrized with the correct
    functions that provide functionality. See `get_` and `read_` functions
    defined in this module.

    Parameters
    ----------
    user_config : dict
        The user configuration dictionary.

    default_groups : dict
        The groups present in the default configuration file for the
        specific module. This dictionary is that created by,
        `get_single_index_groups`, or `get_multiple_index_groups`.

    get_user_groups : callable
        A function to extract groups from a parameter dictionary.
        See `get_single_index_groups` and `get_multiple_index_groups`.

    extract_params : callable
        A function to extract parameters from a parameter dictionary
        according to the selected group. See `extract_single_index_params`,
        `extract_multiple_index_params`.

    Returns
    -------
    set
        A set with the new parameters in the user configuration file
        that are acceptable according to the expandable rules.
    """
    user_groups = get_user_groups(user_config, limit=1)
    default_group_names = [group[0] for group in default_groups]

    new = set()
    for (param_name, group_idx), params in user_groups.items():

        if param_name in default_group_names:
            emsg = _emsg_no_group.format(param_name, group_idx)
            raise ConfigurationError(emsg)

        default_key = (param_name, "1")
        expected_params = default_groups[default_key]
        diff = params.difference(expected_params)
        if diff:
            emsg = _emsg_unexpected_params.format(
                param_name,
                group_idx,
                ", ".join(diff),
                )
            raise ConfigurationError(emsg)

        num_found = len(params)
        num_expected = len(expected_params)

        if num_found != num_expected:
            emsg = _emsg_num_differ.format(
                param_name,
                group_idx,
                num_expected,
                num_found,
                )
            raise ConfigurationError(emsg)

        related_params = extract_params(user_config, param_name, group_idx)
        new.update(related_params)

    return new


def make_param_name_single_index(param_parts):
    """
    Make the key name from param parts.

    For example, ("param", "tag", "1") -> ("param", "1").
    """
    return (param_parts[0], param_parts[-1])


def make_param_name_multiple_index(param_parts):
    """
    Make the key name from param parts.

    For example, ("param", "tag", "2", "1") -> ("param2", "1").
    """
    return (param_parts[0] + param_parts[-2], param_parts[-1])


def belongs_to_single_index(param_parts):
    """
    Assert param belong to single index group.

    Parameter
    ---------
    param_parts : list of str
        The parameter name parts after `str.split("_")`.

    Returns
    -------
    boolean
        True if the parameter name follows the rules of single index
        expandable parameters.

        Rules:
        - param names have more than 2 parts
        - the last part is a digit
        - the part before last is not a digit
    """
    return \
        len(param_parts) > 2 \
        and param_parts[-1].isdigit() \
        and not param_parts[-2].isdigit()


def belongs_to_multiple_index(param_parts):
    """
    Assert param belong to multiple index group.

    Parameter
    ---------
    param_parts : list of str
        The parameter name parts after `str.split("_")`.

    Returns
    -------
    boolean
        True if the parameter name follows the rules of multiple index
        expandable parameters.

        Rules:
        - param names have more than 4 parts
        - the last part is a digit
        - the part before last is a digit
    """
    return \
        len(param_parts) > 4 \
        and param_parts[-1].isdigit() \
        and param_parts[-2].isdigit()


def rejoin_parts_single_index(param_parts):
    """Join parameter name parts."""
    return "_".join(param_parts[1: -1])


def rejoin_parts_multiple_index(param_parts):
    """Join parameter name parts."""
    return "_".join(param_parts[1: -2])


def extract_single_index_params(user_config, param_name, group_idx):
    """Extract the parameters belonging to a group."""
    new = set()
    for param in user_config:
        try:
            pname, *_, pidx = param.split("_")
        except ValueError:
            continue
        if pname == param_name and pidx == group_idx:
            new.add(param)
    return new


def extract_multiple_index_params(user_config, param_name, group_idx):
    """Extract the parameters belonging to a group."""
    new = set()
    for param in user_config:
        try:
            pname, *_, molidx, pidx = param.split("_")
        except ValueError:
            continue
        if pname + molidx == param_name and pidx == group_idx:
            new.add(param)
    return new


get_single_index_groups = partial(
    _get_groups,
    belongs_to_group=belongs_to_single_index,
    make_param_name=make_param_name_single_index,
    rejoin_parts=rejoin_parts_single_index,
    _emsg_no_group=_emsg_no_group_single,
    _emsg_unexpected_params=_emsg_unexpected_params_single,
    )
"""
Get single indexed blocks from a configuration.

Block parameters follow the rule <preffix>_<something>_<group>.

These belong to the same parameter (`param`) of the group `1`:

- param_something1_1
- param_something2_1
- param_something_else_1
- param_something4_1

You could expand these with:

- param_something1_2
- param_something2_2
- param_something_else_2
- param_something4_2

When used to read the modules' default configuration we expect the
<group> to be only "_1". But having <group> allows to identify
blocks from the user configuration.

When we execute this function, we want to know:

- param, the parameter name preffix
- how many parameters exist of the block (`something1`, `something2`, ...)
- _1, defines the number of the block

Parameters
----------
config : dictionary
    Where keys are the parameter names. The values are not important.

minimum : int
    Consider only the groups with at least `minimum` parameters.

Returns
-------
dictionary
    In the form of:
    {
        ("param", "1"): {
            "something1", "something2",
            "something_else", "something4"
            },
        ("otherparam", "1"): {"somename", "othername"},
        }
"""

get_multiple_index_groups = partial(
    _get_groups,
    belong_to_group=belongs_to_multiple_index,
    make_param_name=make_param_name_multiple_index,
    rejoin_parts=rejoin_parts_multiple_index,
    _emsg_no_group=_emsg_no_group_multiple,
    _emsg_unexpected_params=_emsg_unexpected_params_multiple,
    )
"""
Get parameter blocks with multiple indexes.

These blocks of params apply, for example, to different molecules.
That's why they have two looping indexes: one for the molecule and
other for the group.

This block of parameters follow the rule:

<param>_<something>_<N>_<G>

Where the actual param name is `paramN`, where `N` is an integer.

<something> be can any combination of alphanumeric chars and
underscores.

<G> is the number of the group, the will be incremented as we
expand.

These belong to the same parameter (`param1`) of the group `1`:
- param_something_1_1
- param_something2_1_1
- param_something_else_1_1
- param_something4_1

You could expand these with:
- param_something_1_2
- param_something2_1_2
- param_something_else_1_2
- param_something4_2

When used to read the modules' default configuration we expect the
<digit> to be only "_1". But having <digit> allows to identify
blocks from the user configuration.

We want to know:

Returns
-------
dictionary
    In the form of:
    {"param1": {
        "something1", "something2", "something_else",
        "something4"}
"""

read_single_groups_user_config = partial(
    _read_groups_in_user_config,
    get_user_groups=get_single_index_groups,
    )
"""
Read single indexed groups in user config.

Parameters
----------
user_config : dict
    The user configuration dictionary.

default_groups : dict
    The groups present in the default configuration file for the
    specific module. This dictionary is that created by,
    `get_single_index_groups`, or `get_multiple_index_groups`.

Returns
-------
set
    A set with the new parameters in the user configuration file
    that are acceptable according to the expandable rules.
"""

read_multiple_groups_user_config = partial(
    _read_groups_in_user_config,
    get_user_groups=get_multiple_index_groups,
    )
"""
Read multiple indexed groups in user config.

Parameters
----------
user_config : dict
    The user configuration dictionary.

default_groups : dict
    The groups present in the default configuration file for the
    specific module. This dictionary is that created by,
    `get_single_index_groups`, or `get_multiple_index_groups`.

Returns
-------
set
    A set with the new parameters in the user configuration file
    that are acceptable according to the expandable rules.
"""
