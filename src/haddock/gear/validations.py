"""Module to define specific validations in Haddock3."""
from haddock.core.defaults import RUNDIR, valid_run_dir_chars
from haddock.core.exceptions import ConfigurationError


def v_rundir(rundir):
    """Validate string defining the run directory."""
    if set(str(rundir)) - set(valid_run_dir_chars):
        emsg = (
            f"The {RUNDIR!r} parameter can only have "
            r"[a-zA-Z0-9._-/\] characters."
            )
        raise ConfigurationError(emsg)
