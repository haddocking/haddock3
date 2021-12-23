"""Relates to logic or definition of parameters."""

config_mandatory_general_parameters = {
    'molecules',
    'run_dir',
    }
"""The mandatory general arguments of the configuration file."""

non_mandatory_general_parameters_defaults = {
    'self_contained': False,
    'ncores': 8,
    'cns_exec': "",
    'mode': 'local',
    'relative_envvars': False,
    }
