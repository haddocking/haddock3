"""
In-house implementation of toml config files for HADDOCK3.

It does not implement all features of TOML files, but it does implement
the features needed for HADDOCK3. Reported differences with TOML and
specifications:

* Accepts repeated keys.
  * Repetitions get `.#` suffixes, where `#` is integer
* Does not allow comments after definitions (in the same line)
* Values in lists are parsed with `ast.literal_eval()`
* Parsing of values outside lists attempt:
  * int
  * float
  * datetime
  * bool
  * str
* lists can be defined in single lines or multi lines
   * multi-line lists must be followed by an empty line.
"""
import ast
from datetime import datetime


# string separator for internal reasons
_CRAZY = '_dnkljkqwdkjwql_'


def read_config(f):
    """Parse HADDOCK3 config file to a dictionary."""
    d = {}

    # fin can't be processed to a list at this stage
    # it must be kept as a generator
    fin = open(f)

    pure_lines = filter(_is_correct_line, map(_process_line, fin))
    key = None
    for line in pure_lines:

        # dictionary key found
        # gets keys
        if line.startswith('[') and line.endswith(']'):
            key = line.strip('[]')
            if key in d:
                key += _CRAZY + str(sum(1 for _k in d if _k.startswith(key)))

            d.setdefault(key, {})
            continue

        vk, v = (_.strip() for _ in line.split('='))

        # is the start of a multi-line list
        if line[0].isalpha() and '=' in line and not '[' in line:
            d1 = d.setdefault(key, {}) if key else d
            d1[vk] = _parse_value(v)

        # is a line with a list
        elif line[0].isalpha() and '=' in line and '[' in line and ']' in line:
            list_ = _eval_list_str(v)
            d1 = d.setdefault(key, {}) if key else d
            d1[vk] = list_

        elif line[0].isalpha() and '=' in line and '[' in line:
            idx = line.find('[')
            # need to send `fin` and not `pure_lines`
            block = line[idx:] + _get_list_block(fin)
            list_ = _eval_list_str(block)
            d1 = d.setdefault(key, {}) if key else d
            d1[vk] = list_

        else:
            raise ValueError(f'Can\'t process this line: {line!r}')

    fin.close()
    return _make_nested_keys(d)


def _is_correct_line(line):
    return not _is_comment(line)


def _is_comment(line):
    """Assert if line is a comment or a line to ignore."""
    is_comment = \
        line.startswith('#') \
        or not bool(line)

    return is_comment


def _process_line(line):
    return line.strip()


def _parse_value(raw):
    """Parse values as string to its value."""
    # ast.literal_eval is not used to allow special methods
    types = [
        int,
        float,
        _read_bool,
        datetime.fromisoformat,
        _clean_string,
        ]

    for method in types:
        try:
            return method(raw)
        except (ValueError, KeyError):
            continue


def _read_bool(raw):
    """
    Convert "true" and "false" related-strings to `True` and `False`.

    Parameters
    ----------
    raw : str
        Its lower() must be "true" or "false".

    Raises
    ------
    KeyError
        if string not valid.

    Returns
    -------
    bool
    """
    _ = {
        'true': True,
        'false': False,
        }
    return _[raw.lower()]


def _clean_string(s):
    """
    Clean string.

    Operations performed
    --------------------

    * removes quotation marks from strings
    """
    return s.strip("'\"")


def _eval_list_str(s):
    """Evaluates a string to a list."""
    s = '[' + s.strip(',[]') + ']'
    s = s.replace('true', 'True')
    s = s.replace('false', 'False')
    return ast.literal_eval(s)


def _get_list_block(fin):
    """
    Concatenates a list in a string from a multiline list.

    Considers that a list ends when an empty line is found.
    """
    block = []
    for line in fin:

        line = _process_line(line)
        if line.startswith('#'):
            continue

        if not line:
            break

        block.append(line)

    return ''.join(block).strip()


def _make_nested_keys(d):
    """
    Converts dotted nested keys into actual nested keys.

    Example
    -------
    {"one.two": {'foo': 'bar'}}
    {"one": {"two": {'foo': 'bar'}}}

    {"one.two": {...}, "one.foo": {...}}
    {"one": {"two": {...}}, {"foo": {...}}}
    """
    d1 = {}

    for key, value in d.items():

        if isinstance(value, dict):

            keys = [_k.replace(_CRAZY, '.') for _k in key.split('.')]
            dk = d1.setdefault(keys[0], {})

            for k in keys[1:]:
                dk = dk.setdefault(k, {})

            dk.update(value)

        else:
            d1[key] = value

    return d1


def get_module_name(name):
    """Gets the name according to the config parser."""
    return name.split('.')[0]
