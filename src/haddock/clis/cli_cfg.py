"""
Extract the default configuration file for each module.

By default presents the parameters from the three expertise levels:
    `basic`, `intermediate`, and `guru`.

Usage::

    haddock3-cfg -m MODULE
    haddock3-cfg -m MODULE -l LEVEL
    haddock3-cfg -m topoaa -l all
    haddock3-cfg -m rigidbody -l all -g
    haddock3-cfg -G
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

_group = ap.add_argument_group(
    "Mandatory parameters",
    "These parameters are mutually exclusive.",
    )
group = _group.add_mutually_exclusive_group(required=True)
groupm = ap.add_argument_group(
    title="Additional arguments when selecting a module",
    description=(
        "Use these optional arguments together with the `-m` option. "
        "If `-G` is given, these arguments will be ignored."
        ),
    )

group.add_argument(
    "-m",
    dest="module",
    help="The module for which you want to retrieve the default configuration.",
    choices=sorted(modules_category.keys()),
    default=None,
    )

group.add_argument(
    "-G",
    "--only-globals",
    dest="only_globals",
    help="Retrieve only the optional module's general parameters.",
    action="store_true",
    )

groupm.add_argument(
    "-l",
    dest="explevel",
    required=False,
    help="The expertise level of the parameters. Defaults to \"all\".",
    default="all",
    choices=config_expert_levels + ("all",),
    )

groupm.add_argument(
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


def main(module=None, explevel="all", add_globals=True):
    """
    Extract the defaults in the form of a run configuration file.

    Parameters
    ----------
    module : str or None
        The module name to extract the defaults from. ``None`` can be given
        if ``only_globals`` is ``True``.

    explevel : str
        Filter the parameters according to the expert level. Output all
        parameters that belong to the referred level or below. Choices
        are: all, easy, expert, and guru.

    add_globals : bool
        Whether to add the module's general global parameter. If this
        option is given and ``module`` is ``None``, outputs only the
        general parameters.
    """
    from haddock import modules_defaults_path
    from haddock.gear.yaml2cfg import yaml2cfg_text
    from haddock.libs.libio import read_from_yaml

    new_config = ''

    if add_globals:
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
            ))

    if module:

        module_name = ".".join((
            'haddock',
            'modules',
            modules_category[module],
            module,
            ))

        module_lib = importlib.import_module(module_name)
        cfg = module_lib.DEFAULT_CONFIG

        ycfg = read_from_yaml(cfg)
        module_config = yaml2cfg_text(ycfg, module, explevel)

        new_config = os.linesep.join((new_config, module_config))

    sys.stdout.write(new_config, flush=True)

    return 0


if __name__ == "__main__":
    sys.exit(maincli())
