"""Tools for functional programming."""
from operator import is_not


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


def is_str_int(value):
    """Return True if string is convertable to int."""
    try:
        int(value)
    except (ValueError, TypeError):
        return False
    else:
        return True


def is_str_float(value):
    """Return True if string is convertable to float."""
    try:
        float(value)
    except (ValueError, TypeError):
        return False
    else:
        return True
