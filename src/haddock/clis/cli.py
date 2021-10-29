#!/usr/bin/env python3
import argparse
import logging
import sys
from argparse import ArgumentTypeError
from functools import partial

from haddock import current_version
from haddock.libs.libutil import file_exists
from haddock.gear.restart_run import add_restart_arg


# Command line interface parser
ap = argparse.ArgumentParser()

_arg_file_exist = partial(
    file_exists,
    exception=ArgumentTypeError,
    emsg="File {!r} does not exist or is not a file.")
ap.add_argument(
    "recipe",
    type=_arg_file_exist,
    help="The input recipe file path",
    )

add_restart_arg(ap)

ap.add_argument(
    "--setup",
    help="Only setup the run, do not execute",
    action="store_true",
    dest='setup_only',
    )

_log_levels = ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL")
ap.add_argument(
    "--log-level",
    default="INFO",
    choices=_log_levels,
    )

ap.add_argument(
    "-v",
    "--version",
    help="show version",
    action="version",
    version=f'{ap.prog} - {current_version}',
    )


def load_args(ap):
    """Load argument parser args."""
    return ap.parse_args()


def cli(ap, main):
    """Command-line interface entry point."""
    cmd = load_args(ap)
    main(**vars(cmd))


def maincli():
    """Main client execution."""
    cli(ap, main)


def main(
        recipe,
        restart=None,
        setup_only=False,
        log_level="INFO",
        ):
    """
    Execute HADDOCK3 client logic.

    Parameters
    ----------
    recipe : str or pathlib.Path
        The path to the recipe (config file).

    restart : int
        At which step to restart haddock3 run.

    setup_only : bool
        Whether to setup the run without running it.

    log_level : str
        The logging level: INFO, DEBUG, ERROR, WARNING, CRITICAL.
    """
    # anti-pattern to speed up CLI initiation
    from haddock.libs.libworkflow import WorkflowManager
    from haddock.gear.greetings import get_adieu, get_initial_greeting
    from haddock.gear.prepare_run import setup_run
    from haddock.core.exceptions import HaddockError, ConfigurationError

    # Configuring logging
    logging.basicConfig(
        level=log_level,
        format="[%(asctime)s] %(name)s:L%(lineno)d %(levelname)s - %(message)s",
        )

    # Special case only using print instead of logging
    logging.info(get_initial_greeting())

    try:
        params, other_params = setup_run(recipe, restart_from=restart)

    except HaddockError as err:
        logging.error(err)
        raise err

    if not setup_only:
        try:
            workflow = WorkflowManager(
                workflow_params=params,
                start=restart,
                **other_params,
                )

            # Main loop of execution
            workflow.run()

        except HaddockError as err:
            raise err
            logging.error(err)

    # Finish
    logging.info(get_adieu())


if __name__ == "__main__":
    sys.exit(maincli())
