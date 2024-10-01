"""
Expandable parameter module.

This module contains all logic for reading and preparing expandable parameters.

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

The above parameters will be accepted and used in the CNS modules. We
say they are the second groupd of the `fle1` group. We explain in detail
below.

Currently, there are three types of expandable parameters:

1 - Parameters following the name structure:

- `name_something_1`,
- `name_something_else_1`

These two parameters belong to the parameter name `name` and group `1`.
The two parameters have the name `something` and `something_else`.
Therefore, you could define also:

- `name_something_2`
- `name_something_else_2`

Numbers are allowed as long as in combination with word characters, for
example: `name_something1_2`, but not `name_something_1_1` because the
latter is of the second type.

2 - Parameters following the name structure:

- `name_something_X_Y`
- `name_else_X_Y`

In these cases, the parameter name is `nameX`, where `X` is an integer.
`Y` is the number of the group, and `something` and `else` are the name
of the parameters. Taking example:

- `fle_sta_1_1`
- `fle_end_1_1`

The parameter name is `fle1` and the parameters belonging to the groups
are `sta` and `end`, and the group is either `1`.

3 - The simplest parameters in the form of:

- `param_1`

In these cases, you can define `param_2`, `param_3`, etc, in the user
configuration file despite those are not defined in the `defaults.yaml`.
However, because `<name>_<integer>` is too much of a simple rule, we
need to define in this module which parameters are actualy expandable.
If you are developing here look for the `type_simplest_ep` variable.

4. Parameters that are expandable to the max number of molecules:

Those are the parameters starting with `mol`. For example:
`mol_fix_origin_1`, which refers to the `fix_origin` parameter for
molecule 1. These parameters are allowed to expand only to the maximum
of input molecules, and at most to the max number of molecules allowed.
"""
import itertools as it
from copy import deepcopy
from functools import partial

from haddock.core.defaults import max_molecules_allowed
from haddock.core.exceptions import ConfigurationError
from haddock.core.typing import (
    Any,
    AnyT,
    Callable,
    Iterable,
    Optional,
    ParamMap,
    Sequence,
    SupportsAddT,
)


# error messages used in _read_groups_in_user_config function
# errors messages for single indexed groups
_emsg_no_group_single = (
    "The parameter block '{}_*_{}' is not a valid expandable parameter."
)
_emsg_unexpected_params_single = (
    "These parameters do not belong to the block '{}_*_{}': {!r}."
)

# errors messages for multiple indexed groups
_emsg_no_group_multiple = (
    "The parameter block '{}_*_{}_*' is not a valid expandable parameter."
)
_emsg_unexpected_params_multiple = (
    "These parameters do not belong to the block '{}_*_{}_*': {!r}."
)

# common error messages
_emsg_num_differ = (
    "The parameter block {!r} expects "
    "{} parameters ({}), but {} are present "
    "in the user configuration file: {}. "
    "Parameter(s) {} are missing."
)


# this dictionary defines which parameters of the form "param_1" are
# expandable, and for which modules. We cannot define a general
# expandable rule for such a simple parameter name because there can be
# many other parameters with the same structure but that are not
# supposed to be expandable.
type_simplest_ep = {
    "topoaa": {
        "hisd",  # hisd_1 should be defined in the `defaults.yaml`
        "hise",  # hise_1 should be defined in the `defaults.yaml`
    },
}


GroupDict = dict[tuple[str, str], Any]


# this is an engine function that must be populated with the respective
# functions the specify its behaviour. This engine is used in
# `get_single_index_groups` and `get_multiple_index_groups`, for
# example. The latter are defined later on with `partial`.
def _get_groups(
    config: ParamMap,
    belongs_to_group: Callable[[Sequence[str]], bool],
    make_param_name: Callable[[Sequence[str]], tuple[str, str]],
    rejoin_parts: Callable[[Sequence[str]], str],
    minimum: int = 2,
    reference: bool = True,
) -> GroupDict:
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

    reference : bool
        Whether reading a default config file or not. If true, only
        groups ending with "1" will be considered. This avoids capturing
        parameters for molecules, such as "mol_fix_origin_1",
        "mol_fix_origin_2", etc.
    """
    splitted = (parameter.split("_") for parameter in config)
    parts = (_parts for _parts in splitted if belongs_to_group(_parts))

    # identify groups in config
    groups: GroupDict = {}
    for p in parts:
        group_identity = make_param_name(p)  # ("param", "1")

        new = groups.setdefault(group_identity, {})

        # counts the number of parameters within the block
        new.setdefault("counts", 0)
        new["counts"] += 1

        # registers the parameters belonging to the group
        new.setdefault("mid", set())
        new["mid"].add(rejoin_parts(p))

    if reference:
        groups = remove_ghost_groups(groups)

    # creates the final dictionary
    final_groups = {k: v["mid"] for k, v in groups.items() if v["counts"] >= minimum}

    return final_groups


def remove_ghost_groups(groups: GroupDict) -> GroupDict:
    """
    Remove ghost groups from dictionary.

    Ghost groups are parameters that are defined in the defaults
    in a manner that they look they are already expanded. For example:

    - mol_shape_1
    - mol_shape_2
    - mol_shape_3
    - mol_fix_origin_1
    - mol_fix_origin_2
    - mol_fix_origin_3

    These integer suffix refer to the input molecules and not to
    expandable groups.

    Parameters
    ----------
    groups : dict
        A dictionary containing the groups. For example as those created
        by `get_single_index_groups`.
    """
    counter: dict[str, int] = {}
    for key in groups:
        counter.setdefault(key[0], 0)
        counter[key[0]] += 1

    new_groups = deepcopy(groups)

    for key, value in counter.items():  # type: ignore
        if value > 1:
            for gk in groups.keys():
                if gk[0] == key:
                    new_groups.pop(gk)

    return new_groups


# this is an engine function that must be populated with the respective
# functions the specify its behaviour. This engine is used in
# `read_single_idx_groups_user_config` and
# `read_multiple_idx_groups_user_config`, for example. The latter are
# defined later on with `partial`.
def _read_groups_in_user_config(
    user_config: Iterable[str],
    default_groups: GroupDict,
    get_user_groups: Callable[..., GroupDict],
    extract_params: Callable[..., set[str]],
    _emsg_no_group: str = "Parameter block is not a valid expandable group",
    _emsg_unexpected_params: str = "Unexpected params for group",
) -> tuple[set[str], dict[str, int]]:
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
    new : set[str]
        A set with the new parameters in the user configuration file
        that are acceptable according to the expandable rules.
    param_name_counts : dict[str, int]
        Count of expendable parameter parameter names
    """
    # minimum=1 is used to capture groups with missing parameters
    user_groups = get_user_groups(user_config, minimum=1, reference=False)
    default_group_names = set(group[0] for group in default_groups)

    new: set[str] = set()
    param_name_counts: dict[str, int] = {}
    for (param_name, group_idx), params in user_groups.items():
        if param_name not in default_group_names:
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

        # it should only trigger when num_found < num_expected
        # when num_found > num_expected the error about in `diff` is
        # triggered instead
        if num_found != num_expected:
            expected_missings = set(expected_params).difference(params)
            emsg = _emsg_num_differ.format(
                param_name,
                num_expected,
                ", ".join(expected_params),
                num_found,
                ", ".join(params),
                ", ".join([f"`{p}`" for p in expected_missings]),
            )
            raise ConfigurationError(emsg)

        related_params = extract_params(user_config, param_name, group_idx)
        # Increment counts for this parameter
        if related_params:
            param_name_counts.setdefault(param_name, 0)
            param_name_counts[param_name] += 1
        new.update(related_params)
    return new, param_name_counts


def read_simplest_expandable(
    config: Iterable[str], expparams: Iterable[str], 
) -> set[str]:
    """
    Read expandable parameters from config file of the type `param_1`.

    Parameters
    ----------
    config : dict, dict.keys, set, or alike
        The user configuration file.
    expparams : dict, dict.keys, set, or alike
        The parameter names that should be considered as expandable.
        Usually, this is a module subdictionary of `type_simplest_ep`.

    Returns
    -------
    set of str
        The parameters in `config` that comply with `expparams`.
    """
    new: set[str] = set()
    for param in config:
        try:
            name, idx = param.split("_")
        except ValueError:
            continue
        if idx.isdigit() and name in expparams:
            new.add(param)
    return new


def get_mol_parameters(config: Iterable[str]) -> set[str]:
    """Identify expandable `mol` parameters."""
    return set(param for param in config if is_mol_parameter(param))


def is_mol_parameter(param: str) -> bool:
    """Identify if a parameter is a `mol` parameter."""
    parts = param.split("_")
    return param.startswith("mol_") and parts[-1].isdigit() and len(parts) > 2


def read_mol_parameters(
    user_config: ParamMap,
    default_groups: Iterable[str],
    max_mols: int = max_molecules_allowed,
) -> set[str]:
    """
    Read the mol parameters in the user_config following expectations.

    Parameters
    ----------
    user_config : dict
        The user configuration dictionary.

    default_groups : dict or set.
        The mol parameters present in the default configuration file for
        the specific module. These are defined by `get_mol_parameters`.

    max_mols : int
        HADDOCK3 has a limit in the number of different molecules it
        accepts for a calculation. Expandable parameters affecting molecules
        should not be allowed to go beyond that number. Defaults to
        `core.default.max_molecules_allowed`.

    Returns
    -------
    set
        The allowed parameters according to the default config and the
        max allowed molecules.
    """
    # removes the integer suffix from the default mol parameters
    default_names = [remove_trail_idx(p) for p in default_groups]

    new: set[str] = set()
    for param in get_mol_parameters(user_config):
        param_name = remove_trail_idx(param)
        param_idx = param.split("_")[-1]
        if param_name in default_names and int(param_idx) <= max_mols:
            new.add(param)
    return new


def remove_trail_idx(param: str) -> str:
    """
    Remove the trailing integer from a parameter.

    If trail is defined by an underscore "_".

    Examples
    --------
    remove_trail_idx('some_parameter_1')
    >>> 'some_parameter'

    remove_trail_idx('some_parameter')
    >>> 'some_parameter'

    remove_trail_idx('parameter')
    >>> 'parameter'

    Parameters
    ----------
    param : str
        A parameter name.
    """
    under_parts = param.rpartition("_")
    if under_parts[2].isdigit():
        return under_parts[0]
    return param


def get_trail_index(param: str) -> Optional[str]:
    """
    Get the trail index if underscored.

    Examples
    --------
    has_trail_index('some_parameter_1')
    >>> '1'

    has_trail_index('some_parameter1')
    >>> None

    has_trail_index('some_parameter_1_1')
    >>> '1'

    has_trail_index('some_parameter_1-1')
    >>> None

    has_trail_index('some_parameter')
    >>> None

    Parameters
    ----------
    param : str
        A parameter name.

    Returns
    -------
    str or None
        The trail index if exist. ``None`` if not found.
    """
    under_parts = param.rpartition("_")
    if under_parts[2].isdigit():
        return under_parts[2]
    return None


def make_param_name_single_index(param_parts: Sequence[AnyT]) -> tuple[AnyT, AnyT]:
    """
    Make the key name from param parts.

    For example, ("param", "tag", "1") -> ("param", "1").
    """
    return (param_parts[0], param_parts[-1])


def make_param_name_multiple_index(
    param_parts: Sequence[SupportsAddT],
) -> tuple[SupportsAddT, SupportsAddT]:
    """
    Make the key name from param parts.

    For example, ("param", "tag", "2", "1") -> ("param2", "1").
    """
    return (param_parts[0] + param_parts[-2], param_parts[-1])


def belongs_to_single_index(param_parts: Sequence[str]) -> bool:
    """
    Assert param belong to single index group.

    Parameters
    ----------
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
    return (
        len(param_parts) > 2
        and param_parts[-1].isdigit()
        and not param_parts[-2].isdigit()
        and not param_parts[0].startswith("mol")
    )


def belongs_to_multiple_index(param_parts: Sequence[str]) -> bool:
    """
    Assert param belong to multiple index group.

    Parameters
    ----------
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
    return (
        len(param_parts) >= 4
        and param_parts[-1].isdigit()
        and param_parts[-2].isdigit()
        and not param_parts[0].startswith("mol")
    )


def rejoin_parts_single_index(param_parts: Sequence[str]) -> str:
    """Join parameter name parts."""
    return "_".join(param_parts[1:-1])


def rejoin_parts_multiple_index(param_parts: Sequence[str]) -> str:
    """Join parameter name parts."""
    return "_".join(param_parts[1:-2])


def extract_single_index_params(
    user_config: Iterable[str], param_name: str, group_idx: str
) -> set[str]:
    """
    Extract the parameters belonging to a group.

    Descriminates between parameters of other expandable groups.
    See: `belongs_to_single_index`.

    Examples
    --------
    >>> ES = extract_single_index_params

    >>> ES({"param_some_1", "param_other_1", "param2", "ppp_ooo"}, "param", "1")
        {"param_some_1", "param_other_1"}

    >>> ES(
        {"param_some_1", "param_other_1", "param_some_2", "param_other_2"},
        "param", "2")
        {"param_some_2", "param_other_2"}

    >>> ES(
        {"param_some_1", "param_other_1", "par_some_2_2", "par_other_2_2"},
        "param", "1")
        {"param_some_1", "param_other_1"}

    >>> ES(
        {"param_some_1", "param_other_1", "par_some_2_2", "par_other_2_2"},
        "par", "2")
        set()

    >>> ES(
        {"param_some_1", "param_other_1", "par_some_2_2", "par_other_2_2"},
        "por", "2")
        set()

    Parameters
    ----------
    user_config : dict, or dict.keys(), or similar
        The user configuration file - parameter keys should be parameter
        names.

    param_name : str
        The parameter name of the group to select. For example: `c3sym`.

    group_idx : str
        The number of the group: "1", "2", "3", etc. Only that group
        number will be selected.

    Returns
    -------
    set of strings
        A set with the selected parameters.
    """
    new: set[str] = set()
    for param in user_config:
        parts = param.split("_")
        try:
            pname, *_, pidx = parts
        except ValueError:
            continue
        if pname == param_name and pidx == group_idx and belongs_to_single_index(parts):
            new.add(param)
    return new


def extract_multiple_index_params(
    user_config: Iterable[str], param_name: str, group_idx: str
) -> set[str]:
    """
    Extract the parameters belonging to a group.

    See also: `belongs_to_multiple_index`.

    Examples
    --------
    >>> EM = extract_multiple_index_params

    >>> EM({"param_some_1", "param_other_1", "param2", "ppp_ooo"}, "param", "1")
        set()

    >>> EM(
        {"param_some_1", "param_other_1", "param_some_2", "param_other_2"},
        "param", "2")
        set()

    >>> EM(
        {"param_some_1", "param_other_1", "par_some_2_2", "par_other_2_2"},
        "param", "1")
        {"par_some_2_2", "par_other_2_2"},

    >>> EM(
        {"param_some_1", "param_other_1", "par_some_2_2", "par_other_2_2"},
        "par", "2")
        set()

    >>> EM(
        {"param_some_1", "param_other_1", "par_some_2_2", "par_other_2_2"},
        "por", "2")
        set()

    Parameters
    ----------
    user_config : dict, or dict.keys(), or similar
        The user configuration file - parameter keys should be parameter
        names.

    param_name : str
        The parameter name of the group to select. For example: `c3sym`.

    group_idx : str
        The number of the group: "1", "2", "3", etc. Only that group
        number will be selected.

    Returns
    -------
    set of strings
        A set with the selected parameters.
    """
    new: set[str] = set()
    for param in user_config:
        parts = param.split("_")
        try:
            pname, *_, molidx, pidx = parts
        except ValueError:
            continue
        if (
            pname + molidx == param_name
            and pidx == group_idx
            and belongs_to_multiple_index(parts)
        ):
            new.add(param)
    return new


# define public API by functionalizing the engines with the specific
# functions. The docstring below is the docstring of the function.
get_single_index_groups = partial(
    _get_groups,
    belongs_to_group=belongs_to_single_index,
    make_param_name=make_param_name_single_index,
    rejoin_parts=rejoin_parts_single_index,
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
blocks from the user configuration. See `reference` parameter.

When we execute this function, we want to know:

- param, the parameter name preffix
- how many parameters exist of the block (`something1`, `something2`, ...)
- _1, defines the number of the block

Returns a dictionary in the form of::

    {
        ("param", "1"): {
            "something1", "something2",
            "something_else", "something4"
            },
        ("otherparam", "1"): {"somename", "othername"},
        }

Parameters
----------
config : dictionary
    Where keys are the parameter names. The values are not important.

minimum : int
    Consider only the groups with at least `minimum` parameters.

reference : bool
    Whether reading a default config file or not. If true, only
    groups ending with "1" will be considered. This avoids capturing
    parameters for molecules, such as "mol_fix_origin_1",
    "mol_fix_origin_2", etc.

Returns
-------
dictionary - {tuple: dict}
"""

get_multiple_index_groups = partial(
    _get_groups,
    belongs_to_group=belongs_to_multiple_index,
    make_param_name=make_param_name_multiple_index,
    rejoin_parts=rejoin_parts_multiple_index,
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

Returns a dictionary in the form of::

    {("param1", "1"): {
        "something1", "something2", "something_else", "something4"}

Parameters
----------
config : dict
    Where keys are the parameter names. The values are not important.

minimum : int
    Consider only the groups with at least `minimum` parameters.

reference : bool
    Whether reading a default config file or not. If true, only
    groups ending with "1" will be considered. This avoids capturing
    parameters for molecules, such as "mol_fix_origin_1",
    "mol_fix_origin_2", etc.

Returns
-------
dictionary : {(tuple): dict}
"""

read_single_idx_groups_user_config = partial(
    _read_groups_in_user_config,
    get_user_groups=get_single_index_groups,
    extract_params=extract_single_index_params,
    _emsg_no_group=_emsg_no_group_single,
    _emsg_unexpected_params=_emsg_unexpected_params_single,
)
"""
Read single indexed groups in user config.

Parameters
----------
user_config : dict
    The user configuration dictionary.

default_groups : dict
    The groups present in the default configuration file for the
    specific module. This is the dictionary created by,
    `get_single_index_groups`.

Returns
-------
set
    A set with the new parameters in the user configuration file
    that are acceptable according to the expandable rules.
"""

read_multiple_idx_groups_user_config = partial(
    _read_groups_in_user_config,
    get_user_groups=get_multiple_index_groups,
    extract_params=extract_multiple_index_params,
    _emsg_no_group=_emsg_no_group_multiple,
    _emsg_unexpected_params=_emsg_unexpected_params_multiple,
)
"""
Read multiple indexed groups in user config.

Parameters
----------
user_config : dict
    The user configuration dictionary.

default_groups : dict
    The groups present in the default configuration file for the
    specific module. This is the dictionary created by
    `get_multiple_index_groups`.

Returns
-------
set
    A set with the new parameters in the user configuration file
    that are acceptable according to the expandable rules.
"""


def populate_mol_parameters_in_module(
    params: ParamMap, num_mols: int, defaults: ParamMap
) -> None:
    """
    Populate parameters dictionary with the needed molecule `mol_` parameters.

    The `mol_` prefixed parameters is a subclass of the expandable parameters.

    See :py:mod:`haddock.gear.expandable_parameters`.

    Modules require these parameters to be repeated for the number of input
    molecules.

    This function adds `mol_` parameters to the user input parameters,
    one per each `molecule` and for those which a values has not been
    added yet.

    Parameters
    ----------
    modules_params : dict
        A dictionary of parameters.

    Returns
    -------
    None
        Alter the dictionary in place.
    """
    # filter to get only the `mol_` parameters that are defaults
    # the get_trail_index is needed because defaults only have the '1'
    # suffixed parameters.
    mol_params = (
        p
        for p in list(params.keys())
        if is_mol_parameter(p) and get_trail_index(p) == "1"
    )

    for param, i in it.product(mol_params, range(1, num_mols + 1)):
        param_name = remove_trail_idx(param)
        params.setdefault(f"{param_name}_{i}", defaults[param])

    return
