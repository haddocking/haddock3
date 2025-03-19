"""Module to define specific validations in Haddock3."""
from haddock import config_expert_levels, _hidden_level
from haddock.core.defaults import RUNDIR, valid_run_dir_chars
from haddock.core.exceptions import ConfigurationError
from haddock.core.typing import FilePath, Any, ParamDict
from haddock.libs.libio import read_from_yaml, check_yaml_duplicated_parameters

_allowed_expert_levels = config_expert_levels + ("all", _hidden_level)
YAML_PARAMETER_BASE_SCHEME = (
    'title', 'short', 'long', 'group', 'explevel', 'type', 'default',
    )
TYPE_SPECIFIC_SCHEMES = {
    "integer": ("min", "max", ),
    "float": ("min", "max", "precision", ),
    "string": ("minchars", "maxchars", ),
    "list": ("minitems", "maxitems", ),
    "boolean": (),
    "file": (),
    "dict": (),
    }


def v_rundir(rundir: FilePath) -> None:
    """Validate string defining the run directory."""
    if set(str(rundir)) - set(valid_run_dir_chars):
        emsg = (
            f"The {RUNDIR!r} parameter can only have "
            r"[a-zA-Z0-9._-/\] characters."
            )
        raise ConfigurationError(emsg)


def validate_yaml_params_scheme(yaml_fpath: FilePath) -> None:
    """Validate a defaults.yaml file module parameters schemes.

    Parameters
    ----------
    yaml_fpath : str
        Path to the defaults.yaml file to check.
    """
    ycfg = read_from_yaml(yaml_fpath)
    # Loop over parameters
    for param_name, parameters in ycfg.items():
        try:
            validate_parameter_scheme(param_name, parameters)
        except AssertionError as error:
            raise AssertionError(f"{error} (in {yaml_fpath})")


def validate_parameter_scheme(
        param_name: str,
        parameters: ParamDict,
        ) -> None:
    """Validate a parameter scheme.

    Parameters
    ----------
    param_name : str
        Name of this parameter.
    parameters : dict[str, Any]
        Dictionary of param, value for this parameter.
    """
    for base_subparam in YAML_PARAMETER_BASE_SCHEME:
        # Check that this parameter contain the basic subparameters
        assert base_subparam in parameters.keys(), f"Parameter '{param_name}' is missing base subparam '{base_subparam}'"  # noqa : E501
        # Point avlue of this parameter
        param_value = parameters[base_subparam]
        # Evaluate expertize level
        if base_subparam == "explevel":
            # Check if expert level is ok
            assert parameters[base_subparam] in _allowed_expert_levels, f"Parameter '{param_name}' do not contain appropriate expert level"  # noqa : E501
        # Once we know the type, we can make parameter type specific checks
        elif base_subparam == "type":
            # Check if known type
            assert param_value in TYPE_SPECIFIC_SCHEMES.keys(), f"Parameter '{param_name}' contain an unknown parameter type: '{param_value}'"  # noqa : E501
            # Skip subsequent checks in case of type: dict
            if param_value == "dict":
                # Loop over nested parameters
                for dict_param_name, dict_param in parameters.items():
                    # Make sure they are parameters and not base ones
                    if dict_param_name not in YAML_PARAMETER_BASE_SCHEME:
                        # Recursive calls for nested parameters
                        validate_parameter_scheme(dict_param_name, dict_param)
                # Special case of dict, they do not have a default value...
                # So we skip the presence of 'default' in this type of params.
                # We break the loop before this check
                break

            # Check if all specific param for this type of parameter
            # are properly specified
            for specific_param in TYPE_SPECIFIC_SCHEMES[param_value]:
                assert specific_param in parameters.keys(), f"Parameter {param_name} is missing specific subparam '{specific_param}' for type '{base_subparam}'"  # noqa : E501


def validate_defaults_yaml(yaml_fpath: FilePath) -> None:
    """Validate a defaults.yaml file.

    Parameters
    ----------
    yaml_fpath : str
        Path to the defaults.yaml file to validate.
    """
    try:
        check_yaml_duplicated_parameters(yaml_fpath)
        validate_yaml_params_scheme(yaml_fpath)
    except AssertionError as assert_error:
        raise ConfigurationError(assert_error)
