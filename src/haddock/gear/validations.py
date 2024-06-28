"""Module to define specific validations in Haddock3."""
from haddock import config_expert_levels, _hidden_level
from haddock.core.defaults import RUNDIR, valid_run_dir_chars
from haddock.core.exceptions import ConfigurationError
from haddock.core.typing import FilePath
from haddock.libs.libio import read_from_yaml, check_yaml_duplicated_parameters


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


def validate_yaml_params_scheme(yaml_fpath: str) -> None:
    """Validate a default.yaml file parameters schemes.

    Parameters
    ----------
    yaml_fpath : str
        Path to the default.yaml file
    """
    ycfg = read_from_yaml(yaml_fpath)
    allowed_expert_levels = config_expert_levels + ("all", _hidden_level)
    # Loop over parameters
    for param, parameters in ycfg.items():
        for base_subparam in YAML_PARAMETER_BASE_SCHEME:
            # Check that this parameter contain the basic subparameters
            assert base_subparam in parameters.keys(), f"Parameter '{param}' is missing base subparam '{base_subparam}' (in {yaml_fpath})."  # noqa : E501
            param_value = parameters[base_subparam]
            if base_subparam == "explevel":
                # Check if expert level is ok
                assert parameters[base_subparam] in allowed_expert_levels, f"Parameter '{param}' do not contain appropriate expert level (in {yaml_fpath})"  # noqa : E501
            
            elif base_subparam == "type":
                # Check if known type
                assert param_value in TYPE_SPECIFIC_SCHEMES.keys(), f"Parameter {param} contain an unknown parameter type (in {yaml_fpath})"  # noqa : E501
                # Skip subsequent checks in case of type: dict
                if param_value == "dict":
                    break

                # Check if all specific param for this type of parameter
                # are properly specified
                for specific_param in TYPE_SPECIFIC_SCHEMES[param_value]:
                    assert specific_param in parameters.keys(), f"Parameter {param} is missing specific subparam '{specific_param}' for type '{base_subparam}' (in {yaml_fpath})"  # noqa : E501


def validate_defaults_yaml(yaml_fpath: str) -> None:
    """Validate a default.yaml file.

    Parameters
    ----------
    yaml_fpath : str
        Path to the default.yaml file
    """
    try:
        check_yaml_duplicated_parameters(yaml_fpath)
        validate_yaml_params_scheme(yaml_fpath)
    except AssertionError as assert_error:
        raise ConfigurationError(assert_error)
