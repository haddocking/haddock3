"""General modules definitions."""
from pathlib import Path

from haddock import modules_defaults_path
from haddock.gear.yaml2cfg import read_from_yaml_config


modules_folder = Path(__file__).resolve().parent

_folder_match_regex = '[a-zA-Z]*/'
modules_category = {
    module.name: category.name
    for category in modules_folder.glob(_folder_match_regex)
    for module in category.glob(_folder_match_regex)
    }
"""Indexes each module in its specific category. Keys are Paths to the module,
values are their categories. Categories are the modules parent folders."""

# this dictionary defines non-mandatory general parameters that can be defined
# as global parameters thus affect all modules, or, instead, can be defined per
# module where the module definition overwrites global definition. Not all
# modules will use these parameters. It is the responsibility of the module to
# extract the parameters it needs.
# the config file is in modules/defaults.cfg
non_mandatory_general_parameters_defaults = \
    read_from_yaml_config(modules_defaults_path)
