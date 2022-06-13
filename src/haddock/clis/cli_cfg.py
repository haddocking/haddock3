"""
Extract the default configuration file for each module.

By default presents the parameters from the three expertise levels:
    `basic`, `intermediate`, and `guru`.

Usage::

    haddock3-cfg -m MODULE
    haddock3-cfg -m MODULE -l LEVEL
    haddock3-cfg -m topoaa -l all
    haddock3-cfg -m rigidbody -l all -g
"""
import argparse
import importlib
import os
import sys

from haddock import config_expert_levels
from haddock.modules import modules_category


ap = argparse.ArgumentParser(
    prog="HADDOCK3 config retriever",
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    "-m",
    dest="module",
    help="The module for which you want to retrieve the default configuration.",
    required=True,
    choices=sorted(modules_category.keys()),
    )

ap.add_argument(
    "-l",
    dest="explevel",
    required=False,
    help="The expertise level of the parameters. Defaults to \"all\".",
    default="all",
    choices=config_expert_levels + ("all",),
    )

ap.add_argument(
    '-g',
    dest='add_global',
    help="Add also the optional module's general parameters.",
    action="store_true",
    )


def _ap():
    return ap


# command-line client helper functions
# load_args, cli, maincli
def load_args(ap):
    """Load argument parser args."""
    return ap.parse_args()


def cli(ap, main):
    """Command-line interface entry point."""
    cmd = load_args(ap)
    main(**vars(cmd))


def maincli():
    """Execute main client."""
    cli(ap, main)


def main(module, explevel, add_global=True):
    """Extract the default configuration file for a given module."""
    from haddock import modules_defaults_path
    from haddock.gear.yaml2cfg import yaml2cfg_text
    from haddock.libs.libio import read_from_yaml

    module_name = ".".join((
        'haddock',
        'modules',
        modules_category[module],
        module,
        ))

    module_lib = importlib.import_module(module_name)
    cfg = module_lib.DEFAULT_CONFIG

    ycfg = read_from_yaml(cfg)
    new_config = yaml2cfg_text(ycfg, module, explevel)

    if add_global:
        general_cfg = read_from_yaml(modules_defaults_path)
        general_params_str = yaml2cfg_text(
            general_cfg,
            module=None,
            explevel="all",
            )
        comment = os.linesep.join((
            "# The parameters below are optional parameters. ",
            "# They can either be used as global parameters or as part ",
            "# of the module's parameters",
            ))

        new_config = os.linesep.join((
            comment,
            general_params_str,
            new_config,
            ))

    print(new_config, file=sys.stdout, flush=True)  # noqa: T201

    return 0


if __name__ == "__main__":
    sys.exit(maincli())
