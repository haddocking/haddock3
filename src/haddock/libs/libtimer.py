"""Tools related to timing functions."""
from contextlib import contextmanager
from time import time
from typing import Generator

from haddock import log


@contextmanager
def log_time(pre_msg: str) -> Generator[None, None, None]:
    """
    Log the time taken to execute the code under the context.

    Examples
    --------
    >>> with log_time("function took"):
        do_something()
        do_more()

    logged: "function tooked 10 minutes"

    Parameters
    ----------
    pre_msg : str
        String to log with appended time.

    See Also
    --------
    :py:func:`convert_seconds_to_min_sec`
    """
    start = time()
    yield
    end = time()
    elapsed = convert_seconds_to_min_sec(end - start)
    log.info(f"{pre_msg} {elapsed}")


def convert_seconds_to_min_sec(seconds: float) -> str:
    """
    Convert seconds to min&sec.

    Examples
    --------
    >>> convert_seconds_to_min_sec(60)
    1 minute

    >>> convert_seconds_to_min_sec(120)
    2 minutes and 0 seconds

    >>> convert_seconds_to_min_sec(179)
    2 minutes and 59 seconds

    Parameters
    ----------
    seconds : int
        The elapsed time in seconds. Seconds are round to integers.

    Returns
    -------
    str
    """
    seconds = int(round(seconds, 0))
    hours = seconds // 3600
    minutes = (seconds) // 60 % 60
    seconds = (seconds) % 60

    if hours:
        return f"{hours}h{minutes}m{seconds}s"

    if minutes:
        s = "" if minutes == 1 else "s"
        return f"{minutes} minute{s} and {seconds} seconds"

    s = "" if seconds == 1 else "s"
    return f"{seconds} seconds"
