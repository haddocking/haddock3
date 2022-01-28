"""Tools for functional programming."""


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
