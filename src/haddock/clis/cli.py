#!/usr/bin/env python3
import argparse
import sys
import shutil
from argparse import ArgumentTypeError
from contextlib import contextmanager
from functools import partial
from pathlib import Path

from haddock import current_version, has_terminal, log
from haddock.gear.restart_run import add_restart_arg
from haddock.libs.liblog import add_loglevel_arg, set_log_for_cli
from haddock.libs.libutil import file_exists
from haddock.libs.libio import save_stream_to_file, save_streams_to_files


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


@contextmanager
def manage_config_errors(log_streams, log_names):
    try:
        yield
    except Exception as err:
        _msg = (
            'An error ocurred while reading the configuration file, '
            'hence the log files will be saved in the CWD and not in the '
            'specified run directory.'
            )
        log.info(_msg)
        log.error(err)
        save_streams_to_files(log_streams, log_names)
        sys.exit()


@contextmanager
def manage_run_errors():
    try:
        yield
    except (SystemExit, KeyboardInterrupt) as err:
        log.info('Something external to the code halted the execution.')
        log.exception(err)
    except Exception as err:
        log.info('An error ocurred while running the HADDOCK3 workflow.')
        log.exception(err)


def _run_haddock(restart, params, other_params):

    # anti-pattern to speed up CLI initiation
    from haddock.libs.libworkflow import WorkflowManager

    workflow = WorkflowManager(
        workflow_params=params,
        start=restart,
        **other_params,
        )

    # Main loop of execution
    workflow.run()

    return


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
    # small antipatterns to improve CLI startup speed
    from haddock.gear.greetings import get_adieu, get_initial_greeting
    from haddock.gear.prepare_run import setup_run

    log_streams, log_suggested_names = \
        set_log_for_cli(log, log_level, stream_to_stdout=has_terminal)

    log.info(get_initial_greeting())

    with manage_config_errors(log_streams, log_suggested_names):
        params, other_params = setup_run(recipe, restart_from=restart)

    if not setup_only:
        with manage_run_errors():
            _run_haddock(restart, params, other_params)

    _add_dir_run = lambda x: Path(other_params['run_dir'], Path(x).name)
    _log_names = list(map(_add_dir_run, log_suggested_names))
    save_streams_to_files(log_streams, _log_names)

    # Finish
    log.info(get_adieu())


if __name__ == "__main__":
    sys.exit(maincli())
