"""
Extract the default configuration file for each module.

By default presents the parameters from the three expertise levels:
    `basic`, `intermediate`, and `guru`.

Usage::

    haddock3-cfg  # prints only the general parametes
    haddock3-cfg -g  # same as above
    haddock3-cfg -m MODULE  # prints only the module parameters
    haddock3-cfg -m MODULE -l LEVEL
    haddock3-cfg -m topoaa -l all
    haddock3-cfg -m rigidbody -l all -g  # prints the module and the globals
"""
import argparse
import importlib
import os
import sys

from pathlib import Path

from haddock import config_expert_levels, core_path
from haddock.core.typing import (
    ArgumentParser,
    Callable,
    ExpertLevel,
    Namespace,
    Optional,
)
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
    choices=sorted(modules_category.keys()),
    default=None,
)

ap.add_argument(
    "-l",
    dest="explevel",
    required=False,
    help='The expertise level of the parameters. Defaults to "all".',
    default="all",
    choices=config_expert_levels + ("all",),
)

ap.add_argument(
    "-g",
    "--globals",
    dest="global_params",
    help="Add also the optional module's general parameters.",
    action="store_true",
)

ap.add_argument(
    "-d",
    "--details",
    dest="details",
    help="Add detailed parameter description found in 'long'.",
    action="store_true",
    default=False,
)


def _ap() -> ArgumentParser:
    return ap


# command-line client helper functions
# load_args, cli, maincli
def load_args(ap: ArgumentParser) -> Namespace:
    """Load argument parser args."""
    return ap.parse_args()


def cli(ap: ArgumentParser, main: Callable[..., None]) -> None:
    """Command-line interface entry point."""
    cmd = load_args(ap)

    # I didn't want to to have the `--globals` param as a negative parameter
    # ('store_false'). However, `module` and `globals` can't both be false.
    # In order for the commands below performing the same:
    #
    # haddock3-cfg
    # haddock3-cfg -g
    #
    # we need this quick if statement. Which basically adjust the argparse to
    # the default values in the main function.
    # @joaomcteixeira
    if cmd.global_params is False and cmd.module is None:
        cmd.global_params = True

    main(**vars(cmd))


def maincli() -> None:
    """Execute main client."""
    cli(ap, main)


def main(
    module: Optional[str] = None,
    explevel: str = "all",
    global_params: bool = True,
    details: bool = False,
) -> None:
    """
    Extract the defaults in the form of a run configuration file.

    Parameters
    ----------
    module : str or None
        The module name to extract the defaults from.
        If ``None`` given, we expect ``global_params`` to be ``True``,
        otherwise nothing will be printed.

    explevel : str
        Filter the parameters according to the expert level. Output all
        parameters that belong to the referred level or below. Choices
        are: all, easy, expert, and guru.

    global_params : bool
        Whether to add the module's general global parameter. If
        ``True`` and ``module`` is ``None``, outputs only the general
        parameters.

    details : bool
        Whether to add the 'long' description of each parameter.
    """
    from haddock import modules_defaults_path
    from haddock.gear.yaml2cfg import yaml2cfg_text
    from haddock.libs.libio import read_from_yaml

    new_config = ""

    if global_params:
        # Read general parameters
        general_cfg = read_from_yaml(modules_defaults_path)
        general_params_str = yaml2cfg_text(
            general_cfg,
            module=None,
            explevel="all",
            details=details,
            )
        general_comment = os.linesep.join(
            (
                "#" * 50,
                "# The parameters below are optional parameters. ",
                "# They can either be used as global parameters or as part ",
                "# of the module's parameters",
                "#" * 50,
                )
            )

        # Read mandatory parmeters
        mandatory_cfg_path = Path(core_path, "mandatory.yaml")
        general_mandatory_cfg = read_from_yaml(mandatory_cfg_path)
        general_mandatory_params_str = yaml2cfg_text(
            general_mandatory_cfg,
            module=None,
            explevel="all",
            details=details,
            mandatory_param=True,
            )
        mandatory_comment = os.linesep.join(
            (
                "#" * 50,
                "# The parameters below are mandatory parameters.",
                "# They must be specified in your configuration file!",
                "#" * 50,
                )
            )

        # Read optional parameters
        optional_cfg_path = Path(core_path, "optional.yaml")
        general_optional_cfg = read_from_yaml(optional_cfg_path)
        general_optional_params_str = yaml2cfg_text(
            general_optional_cfg,
            module=None,
            explevel="all",
            details=details,
            )
        optional_comment = os.linesep.join(
            (
                "#" * 50,
                "# The parameters below are optional global parameters.",
                "#" * 50,
                )
            )
        # Concatenate all of them
        new_config = os.linesep.join(
            (
                mandatory_comment,
                general_mandatory_params_str,
                optional_comment,
                general_optional_params_str,
                general_comment,
                general_params_str,
                )
            )

    if module:
        module_name = ".".join(
            (
                "haddock",
                "modules",
                modules_category[module],
                module,
            )
        )

        module_lib = importlib.import_module(module_name)
        cfg = module_lib.DEFAULT_CONFIG

        ycfg = read_from_yaml(cfg)
        module_config = yaml2cfg_text(ycfg, module, explevel, details=details)

        new_config = os.linesep.join((new_config, module_config))

    sys.stdout.write(new_config)

    return


if __name__ == "__main__":
    sys.exit(maincli())  # type: ignore
