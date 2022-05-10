"""Relates to logic or definition of parameters."""
from pathlib import Path

from haddock import haddock3_source_path
from haddock.gear.yaml2cfg import read_from_yaml_config


MANDATORY_YAML = Path(haddock3_source_path, "core", "mandatory.yaml")
"""The mandatory general arguments of the configuration file."""

OPTIONAL_YAML = Path(haddock3_source_path, "core", "optional.yaml")
"""The optional general arguments of the configuration file."""

_mandatory_parameters = read_from_yaml_config(MANDATORY_YAML)
config_mandatory_general_parameters = set(_mandatory_parameters)

_optional_parameters = read_from_yaml_config(OPTIONAL_YAML)
config_optional_general_parameters = set(_optional_parameters)
