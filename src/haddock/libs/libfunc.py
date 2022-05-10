"""
Tools for Functional-Programming in Python.

From: https://github.com/joaomcteixeira/libfuncpy
"""
from operator import is_not


def chainf(init, *funcs, **common):
    """
    Apply a sequence of functions to an initial value.

    Example
    -------
    >>> chainf(2, *[str, int, float])
    2.0

    Parameters
    ----------
    init : anything
        The initial value.

    **common : keyword arguments
        Common key word arguments to all functions.

    Returns
    -------
    anything
        The result of the chain of functions; this is, the return value
        of the last function.
    """
    for func in funcs:
        init = func(init, **common)
    return init


def chainfs(*funcs, **common):
    """
    Store functions be executed on a value.

    Example
    -------
    >>> do = chainfs(str, int, float)
    >>> do(2)
    2.0

    See Also
    --------
    :py:func:`chainf`
    """
    def execute(value):
        return chainf(value, *funcs, **common)

    return execute


def give_same(value):
    """Return what is given."""
    return value


def true(*ignore, **everything):
    """Give True regardless of the input."""
    return True


def false(*ignore, **everything):
    """Give False regardless of the input."""
    return False


def none(*ignore, **everything):
    """Give False regardless of the input."""
    return None


def nan(*ignore, **everything):
    """Give False regardless of the input."""
    return float('nan')


def not_none(value):
    """Give True if value is not None, or False otherwise."""
    return is_not(value, None)
