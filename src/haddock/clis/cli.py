#!/usr/bin/env python3
import argparse
import sys
from argparse import ArgumentTypeError
from functools import partial

from haddock import current_version, has_terminal, log
from haddock.gear.restart_run import add_restart_arg
from haddock.libs.liblog import add_loglevel_arg, set_log_level_for_clis
from haddock.libs.libutil import file_exists


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

add_loglevel_arg(ap)

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

    set_log_level_for_clis(log, log_level, add_streamhandler=has_terminal)

    # Special case only using print instead of logging
    log.info(get_initial_greeting())

    try:
        params, other_params = setup_run(recipe, restart_from=restart)

    except ConfigurationError as err:
        log.error(err)
        sys.exit()

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
            log.error(err)

    # Finish
    log.info(get_adieu())


if __name__ == "__main__":
    sys.exit(maincli())
