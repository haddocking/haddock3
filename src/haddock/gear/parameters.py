"""Relates to logic or definition of parameters."""
from haddock.core.defaults import cns_exec


config_mandatory_general_parameters = {
    'molecules',
    'run_dir',
    }
"""The mandatory general arguments of the configuration file."""

non_mandatory_general_parameters_defaults = {
    "self_contained": False,
    "ncores": 8,
    "cns_exec": cns_exec,
    "mode": "local",
    }
