"""Module to define specific validations in Haddock3."""
import sys
from pathlib import Path

from haddock import log
from haddock.core.defaults import (
    RUNDIR,
    max_molecules_allowed,
    valid_run_dir_chars,
    )
from haddock.core.exceptions import ConfigurationError
from haddock.gear.greetings import get_goodbye_help


def v_rundir(rundir):
    """Validate string defining the run directory."""
    if set(str(rundir)) - set(valid_run_dir_chars):
        emsg = (
            f"The {RUNDIR!r} parameter can only have "
            r"[a-zA-Z0-9._-/\] characters."
            )
        raise ConfigurationError(emsg)


def v_maxmolecules(num_molecules):
    """Validate the num of molecules is allowed."""
    if len(num_molecules) > max_molecules_allowed:
        emsg = f"Too many molecules defined, max is {max_molecules_allowed}."
        raise ConfigurationError(emsg)


def v_rundir_exists(run_dir):
    """
    Validate if run directory already exists and is not empty.

    Exist execution if run directory exists.
    """
    _p = Path(run_dir)
    if _p.exists() and len(list(_p.iterdir())) > 0:
        log.info(
            f"The {RUNDIR!r} folder {str(_p)!r} exists and is not empty. "
            "We can't work on it unless you provide the `--restart` "
            "option. If you want to start a run from scratch, "
            "indicate a new folder, manually delete {str(_p)!r} first, "
            "or use the `--restart 0` option."
            )
        sys.exit(get_goodbye_help())
