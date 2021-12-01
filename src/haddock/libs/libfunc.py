"""Tools for functional programming."""


def true(*ignore, **everything):
    """Give True regardless of the input."""
    return True


def false(*ignore, **everything):
    """Give False regardless of the input."""
    return False
