"""
Tools for Functional-Programming in Python.

From: https://github.com/joaomcteixeira/libfuncpy
"""
from functools import partial, reduce


def reduce_helper(value, f, *a, **k):
    """
    Help in `reduce`.

    Helper function when applying `reduce` to a list of functions.

    Parameters
    ----------
    value : anything
    f : callable
        The function to call. This function receives `value` as first
        positional argument.
    *a, **k
        Args and kwargs passed to `f`.
    """
    return f(value, *a, **k)


def chainf(init, *funcs, **common):
    """
    Apply a sequence of functions to an initial value.

    Example
    -------
    >>> chainf(2, [str, int, float])
    2.0
    """
    rh = partial(reduce_helper, **common)
    return reduce(rh, funcs, init)


def chainfs(*funcs):
    """
    Store functions be executed on a value.

    Example
    -------
    >>> do = chainfs(str, int, float)
    >>> do(2)
    2.0
    """
    def execute(value):
        return chainf(value, *funcs)
    return


def true(*ignore, **everything):
    """Give True regardless of the input."""
    return True


def false(*ignore, **everything):
    """Give False regardless of the input."""
    return False
