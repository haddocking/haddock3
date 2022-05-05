#!/usr/bin/env python3
"""
Run HADDOCK3 docking simulation.

This is the main command-line client to run HADDOCK3 docking
simulations. To prepare a simulation, setup configuration file defining
a HADDOCK3 workflow and use this command-line client to execute that
workflow.

Usage::

    haddock3 -h
    haddock3 <CONFIG FILE>
"""
import argparse
import sys
from pathlib import Path

from haddock import log
from haddock.core.defaults import RUNDIR
from haddock.gear.restart_from_dir import (
    RESTART_FROM_DIR_DEFAULT,
    add_restart_from_dir,
    )
from haddock.gear.restart_run import add_restart_arg
from haddock.libs.libcli import add_version_arg, arg_file_exist
from haddock.libs.liblog import add_loglevel_arg


# Command line interface parser
ap = argparse.ArgumentParser()

ap.add_argument(
    "recipe",
    type=arg_file_exist,
    help="The input recipe file path",
    )

add_restart_arg(ap)
add_restart_from_dir(ap)

ap.add_argument(
    "--setup",
    help="Only setup the run, do not execute",
    action="store_true",
    dest='setup_only',
    )

add_loglevel_arg(ap)
add_version_arg(ap)


def _ap():
    return ap


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
        restart_from_dir=RESTART_FROM_DIR_DEFAULT,
        setup_only=False,
        log_level="INFO",
        ):
    """
    Run HADDOCK3 docking simulation.

    Parameters
    ----------
    recipe : str or pathlib.Path
        The path to the recipe (config file).

    restart : int, optional
        At which step to restart haddock3 run.

    setup_only : bool, optional
        Whether to setup the run without running it.

    log_level : str, optional
        The logging level: INFO, DEBUG, ERROR, WARNING, CRITICAL.
    """
    # anti-pattern to speed up CLI initiation
    from time import time

    from haddock.gear.greetings import get_adieu, get_initial_greeting
    from haddock.gear.prepare_run import setup_run
    from haddock.libs.libio import working_directory
    from haddock.libs.liblog import (
        add_log_for_CLI,
        add_stringio_handler,
        log_file_name,
        log_formatters,
        )
    from haddock.libs.libutil import (
        convert_seconds_to_min_sec,
        log_error_and_exit,
        )
    from haddock.libs.libworkflow import WorkflowManager

    start = time()
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

    log.info(get_initial_greeting())

    with log_error_and_exit():
        params, other_params = setup_run(
            recipe,
            restart_from=restart,
            restart_from_dir=restart_from_dir,
            )

    # here we the io.StringIO handler log information, and reset the log
    # handlers to fit the CLI and HADDOCK3 specifications.
    log_temporary = log.handlers[-1].stream.getvalue()
    _run_dir = other_params[RUNDIR]
    log_file = Path(_run_dir, log_file_name)
    add_log_for_CLI(log, log_level, log_file)

    # here we append the log information in the previous io.StringIO()
    # handler in the log file already in the run_dir
    with open(log_file, 'a') as fout:
        fout.write(log_temporary)

    if setup_only:
        log.info('We have setup the run, only.')
        log.info(get_adieu())
        return

    with (
            working_directory(other_params[RUNDIR]),
            log_error_and_exit(),
            ):

        workflow = WorkflowManager(
            workflow_params=params,
            start=restart,
            **other_params,
            )

        # Main loop of execution
        workflow.run()

    # Finish
    end = time()
    elapsed = convert_seconds_to_min_sec(end - start)
    log.info(f"This HADDOCK3 run took: {elapsed}")
    log.info(get_adieu())


if __name__ == "__main__":
    sys.exit(maincli())
