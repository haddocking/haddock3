"""
Handles zero filling prefix for module folder names.

To facilitate folder sorting, HADDOCK3 module run folder names have a
integer prefix starting at 0.

Depending on the number of modules of the run, the integer prefix is a
single digit number or multi-digit number (00, 01, 02, 03...).

This gear module contains the "zero filling" logic in HADDOCK 3.

The main class of this module `zerofill` is already instantiated such
that its state can be used and edited throughout the HADDOCK3 library.

Examples
--------
>>> from haddock.gear import zerofill
>>> zerofill.read(modules)  # a dictionary containing 11 modules
>>> zerofill.zfnum
2
>>> zerofill.fill("topoaa", 0)
"00_topoaa"
"""
from math import ceil, log10


class _ZeroFill:
    """
    Register the zero fill number.

    Examples
    --------
    >>> _ZeroFill()
    >>> zerofill.read(modules)  # a dictionary containing five modules
    >>> zerofill.zfnum
    1
    >>> zerofill.fill("topoaa", 0)
    "0_topoaa"
    """

    def __init__(self, zfnum=None):
        self.zfnum = zfnum
        return

    def read(self, modules):
        """
        Register the zerofill number for current run.

        Zero fill number if available through the attribute: `zfnum`.

        Parameters
        ----------
        iterable
            Usually a dictionary or a list of the modules names.
            Any object implementing `len`.
        """
        self.zfnum = get_zerofill_for_modules(modules)

    def fill(self, name, num):
        """
        Fill a name with the corresponding zero filling.

        Examples
        --------
        >>> zerofill.zfum
        1

        >>> zerofill.fill("topoaa", 0)
        "0_topoaa"

        >>> zerofill.fill("rigidody", 1)
        "1_rigidbody"
        """
        return make_zero_fill(num, self.zfnum) + "_" + name


zero_fill = _ZeroFill()


def get_number_of_digits(num):
    """
    Get the number of digits of a number.

    10 has two digits.
    100 has three digits.
    """
    # also: return len(str(num)) :-)
    return max(ceil(log10(num + 1)), 1)


def get_zerofill_for_modules(modules):
    """
    Get a the prefix zerofill for modules.

    If there are 5 modules, zerofill digits is 1.

    If there are 10 modules, zerofill digits is 1 because counting
    starts at 0 (for topoaa).

    If there are 100 modules, zerofill digits is 2 because counting
    starts at 0 (for topoaa).

    This function is used in combination with `zero_fill`.
    """
    return get_number_of_digits(len(modules) - 1)


def make_zero_fill(number, digits):
    """Make a number string zero filled to the left."""
    return str(number).zfill(digits)
