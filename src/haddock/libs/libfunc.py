"""
Tools for Functional-Programming in Python.

From: https://github.com/joaomcteixeira/libfuncpy
"""
from functools import reduce


def vartial(func, *args, **keywords):
    """
    Prepare a function with args and kwargs except for the first arg.
    Functions like `functools.partial` except that the resulting
    preprepared function expects the first positional argument.
    Example
    -------
    >>> pow2 = vartial(math.pow, 2)
    >>> pow2(3)
    9
    >>> pow2(4)
    16
    This is different from:
    >>> pow_base_3 = partial(math.pow, 3)
    >>> pow_base_3(2)
    9
    >>> pow_base_3(4)
    81
    """
    def newfunc(value, *fargs, **fkeywords):
        # newkeywords = {**keywords, **fkeywords}
        newkeywords = keywords.copy()
        newkeywords.update(fkeywords)
        return func(value, *args, *fargs, **newkeywords)
    newfunc.func = func
    newfunc.args = args
    newfunc.keywords = keywords
    return newfunc


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


def chainf(init, *funcs):
    """
    Apply a sequence of functions to an initial value.

    Example
    -------
    >>> chainf(2, [str, int, float])
    2.0
    """
    return reduce(reduce_helper, funcs, init)


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



