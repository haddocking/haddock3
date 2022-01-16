"""Module to define specific validations in Haddock3."""
import string


def v_rundir(rundir):
    """Validate string defining the run directory."""
    valid = string.ascii_letters + string.digits + "._-/\\"
    if set(str(rundir)) - set(valid):
        emsg = \
            r"The `run_dir` parameter can only have [a-zA-Z0-9_-/\] characters."
        raise ValueError(emsg)
