"""Relates to logic or definition of parameters."""
from pathlib import Path

from haddock import haddock3_source_path
from haddock.gear.yaml2cfg import read_from_yaml_config


MANDATORY_YAML = Path(haddock3_source_path, "core", "mandatory.yaml")

mandatory_parameters = read_from_yaml_config(MANDATORY_YAML)

config_mandatory_general_parameters = set(mandatory_parameters)
