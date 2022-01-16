"""
Parse data from HADDOCK3 to YAML and YAML to HADDOCK3 and related.

Accross this file you will see references to "yaml" as variables and
function names. In these cases, we always mean the HADDOCK3 YAML
configuration files which have specific keys.
"""
import os

from haddock import config_expert_levels
from haddock.libs.libio import read_from_yaml


def yaml2cfg_text_with_explevels(
        ycfg,
        module,
        expert_levels=config_expert_levels,
        ):
    """
    Convert HADDOCK3 YAML config to HADDOCK3 user config text.

    Adds commentaries with help strings.

    Parameters
    ----------
    ycfg : dict
        The dictionary representing the HADDOCK3 config file.

    module : str
        The module to which the config belongs to.

    expert_levels : list-list
        A list with the expert levels to consider. Defaults to all.
    """
    new_config = []
    new_config.append(f"[{module}]")

    for level in expert_levels:
        new_config.append(_yaml2cfg_text(ycfg[level], module, level))

    return os.linesep.join(new_config) + os.linesep


def _yaml2cfg_text(ycfg, module, subtitle=None):
    """
    Convert HADDOCK3 YAML config to HADDOCK3 user config text.

    Does not consider expert levels.
    See :func:`yaml2cfg_text_with_explevels` instead.

    Parameters
    ----------
    ycfg : dict
        The dictionary representing the HADDOCK3 YAML configuration.
        This configuration should NOT have the expertise levels. It
        expectes the first level of keys to be the parameter name.
    """
    params = []

    if subtitle:
        # empty strings represent new lines
        params.extend(["", f"# {subtitle}", ""])

    for param_name, param in ycfg.items():

        # "default" is not in param when the key points to a subdictionary
        # of parameters.
        if "default" not in param:

            params.append("")  # give extra space
            curr_module = f"{module}.{param_name}"
            params.append(f"[{curr_module}]")
            params.append(_yaml2cfg_text(param, module=curr_module))

        else:

            comment = []
            comment.append(param["hoover"])
            if "min" in param:
                comment.append(f"$min {param['min']} $max {param['max']}")
            if "length" in param:
                comment.append(f"$maxlen {param['length']}")

            params.append("{} = {!r}  # {}".format(
                param_name,
                param["default"],
                " ".join(comment),
                ))

            if param["type"] == "list":
                params.append(os.linesep)

    return os.linesep.join(params)


def read_from_yaml_config(cfg_file):
    """Read config from yaml by collapsing the expert levels."""
    ycfg = read_from_yaml(cfg_file)
    # there's no need to make a deep copy here, a shallow copy suffices.
    cfg = {}
    for level in config_expert_levels:
        cfg.update(flat_yaml_cfg(ycfg[level]))

    return cfg


def flat_yaml_cfg(cfg):
    """Flat a yaml config."""
    new = {}
    for param, values in cfg.items():
        try:
            new[param] = values["default"]
        except KeyError:
            new[param] = flat_yaml_cfg(values)
    return new
