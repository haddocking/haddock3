"""
Tools for functional programming.

See more at https://github.com/joaomcteixeira/libfuncpy
"""


def true(*ignore, **everything):
    """Give True regardless of the input."""
    return True


def false(*ignore, **everything):
    """Give False regardless of the input."""
    return False


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
