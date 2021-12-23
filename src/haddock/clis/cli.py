#!/usr/bin/env python3
"""Main HADDOCK3 command-line client."""
import argparse
import sys
from argparse import ArgumentTypeError
from functools import partial
from pathlib import Path

from haddock import log, version
from haddock.gear.restart_run import add_restart_arg
from haddock.libs.libio import working_directory
from haddock.libs.liblog import (
    add_log_for_CLI,
    add_loglevel_arg,
    add_stringio_handler,
    log_file_name,
    log_formatters,
    )
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
    version=f'{ap.prog} - {version}',
    )


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
    from haddock.core.exceptions import HaddockError
    from haddock.gear.greetings import get_adieu, get_initial_greeting
    from haddock.gear.prepare_run import setup_run
    from haddock.libs.libworkflow import WorkflowManager

    # the io.StringIO handler is a trick to save the log while run_dir
    # is not read from the configuration file and the log can be saved
    # in the final file.
    #
    # Nonetheless, the log is saved in stdout and stderr in case it
    # breaks before the actual logfile is defined.
    # See lines further down on how we resolve this trick
    add_stringio_handler(
        log,
        log_level=log_level,
        formatter=log_formatters[log_level],
        )

    # Special case only using print instead of logging
    log.info(get_initial_greeting())

    try:
        params, other_params = setup_run(recipe, restart_from=restart)

    except HaddockError as err:
        log.error(err)
        raise err

    # here we the io.StringIO handler log information, and reset the log
    # handlers to fit the CLI and HADDOCK3 specifications.
    log_temporary = log.handlers[-1].stream.getvalue()
    _run_dir = other_params['run_dir']
    log_file = Path(_run_dir, log_file_name)
    add_log_for_CLI(log, log_level, log_file)

    # here we append the log information in the previous io.StringIO()
    # handler in the log file already in the run_dir
    with open(log_file, 'a') as fout:
        fout.write(log_temporary)

    if setup_only:
        log.info('We have setup the run only.')
        log.info(get_adieu())
        return

    with working_directory(other_params['run_dir']):
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
            log.error(err)
            sys.exit(1)

    # Finish
    log.info(get_adieu())


if __name__ == "__main__":
    sys.exit(maincli())
