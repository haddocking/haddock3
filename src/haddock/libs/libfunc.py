"""
Tools for Functional-Programming in Python.

From: https://github.com/joaomcteixeira/libfuncpy
"""
from operator import is_not
from typing import Any, Callable, Literal


def chainf(init: Any, *funcs: Callable[..., Any], **common: Any) -> Any:
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


def chainfs(*funcs: Callable[..., Any], **common: Any) -> Callable[..., Any]:
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
    def execute(value: Any) -> Any:
        return chainf(value, *funcs, **common)

    return execute


def give_same(value: Any) -> Any:
    """Return what is given."""
    return value


def true(*ignore: Any, **everything: Any) -> Literal[True]:
    """Give True regardless of the input."""
    return True


def false(*ignore: Any, **everything: Any) -> Literal[False]:
    """Give False regardless of the input."""
    return False


def none(*ignore: Any, **everything: Any) -> Literal[None]:
    """Give None regardless of the input."""
    return None


def nan(*ignore: Any, **everything: Any) -> float:
    """Give nan regardless of the input."""
    return float('nan')


def not_none(value: Any) -> bool:
    """Give True if value is not None, or False otherwise."""
    return is_not(value, None)
