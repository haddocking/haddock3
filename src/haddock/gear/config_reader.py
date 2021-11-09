"""
In-house implementation of enhanced toml-like config files for HADDOCK3.

(Likely) It does not implement all features of TOML files, but it does implement
the features needed for HADDOCK3. The config reader:

* Accepts repeated keys.
  * Repetitions get `.#` suffixes, where `#` is an integer

* Allows in-line comments for regex-defined values

* lines with special values not defined by regexes do not accept
  comments

* Regex defined values:
  * strings
  * numbers
  * null/nones
  * lists defined in one lines
  * lists defined in multiple lines (require empty line after the list
    definition)

* Values attempted without regex (do not accept comment in lines):
  * `datetime.fromisoformat`

To see examples of possible configuration files and generated
dictionaries see `test/test_config_reader.py`.
"""
import ast
import re
from contextlib import contextmanager
from datetime import datetime


# line regexes
# https://regex101.com/r/r4BlJf/1
_main_header_re = re.compile(r'^ *\[(\w+)\]')

# https://regex101.com/r/wKm07b/1
# thanks https://stackoverflow.com/questions/39158902
_sub_header_re = re.compile(r'^ *\[(\w+(?:\.\w+)+)\]')

# https://regex101.com/r/q2fuFl/1
_string_re = re.compile(r'''^ *(\w+) *= *("(.*?)"|'(.*?)')''')

# https://regex101.com/r/6X4j7n/1
_number_re = re.compile(r'^ *(\w+) *= *(\-?\d+\.?\d*|\-?\.\d+|\-?\.?\d+E\-?\d+|-?\d+\.?\d*E\d+)(?: |#|$)')

# https://regex101.com/r/K6yXbe/1
_none_re = re.compile(r'^ *(\w+) *= *([Nn]one|[Nn]ull)')

# https://regex101.com/r/YCZSAo/1
_list_one_liner_re = re.compile(r'^ *(\w+) *= *(\[.*\])')

# https://regex101.com/r/bWlaWB/1
_list_multiliner_re = re.compile(r" *(\w+) *= *\[\ *#?[^\]\n]*$")

# https://regex101.com/r/kY49lw/1
_true_re = re.compile(r'^ *(\w+) *= *([tT]rue)')
_false_re = re.compile(r'^ *(\w+) *= *([fF]alse)')


class NoGroupFoundError(Exception):
    pass


class ConfigFormatError(Exception):
    pass


class DuplicatedParameterError(Exception):
    pass


def read_config(f):
    """Parse HADDOCK3 config file to a dictionary."""
    with open(f, 'r') as fin:
        return _read_config(fin)


def _read_config(fin):
    """
    Read the config line by line.

    Parameters
    ----------
    fin : a openned file content generator
        A generator expressing the content of a file in lines. Cannot be
        a list.
    """
    # it must be kept as a generator
    d = {}

    # main keys always come before any subdictionary
    d1 = d

    pure_lines = filter(_is_correct_line, map(_process_line, fin))
    header = None
    for line in pure_lines:

        main_header_group = _main_header_re.match(line)
        sub_header_group = _sub_header_re.match(line)

        if main_header_group:
            header = _update_key_number(main_header_group[1], d)
            d1 = d.setdefault(header, {})

        elif sub_header_group:
            headers = sub_header_group[1].split('.')
            header_ = _update_key_number(headers[0], d, offset=-1)

            if header != header_:
                raise ConfigFormatError(
                    'A subheader cannot be defined before its main header.'
                    )

            d1 = d[header_]
            for head in headers[1:]:
                d1 = d1.setdefault(head, {})

        else:
            value_key, value = _read_value(line, fin)
            if value_key in d1:
                _msg = (
                    f'The parameter {value_key!r} is repeated. '
                    'Repeated parameter names are not allowed '
                    'for the same module.'
                    )
                raise DuplicatedParameterError(_msg)
            d1[value_key] = value


    return _remove_trailing_zeros_in_headers(d)


def _update_key_number(key, d, sep='.', offset=0):
    _num = str(sum(1 for k in d if k.startswith(key)) + offset)
    return key + sep + _num


def _read_value(line, fin):
        # evals if key:value are defined in a single line
        try:
            key, value = _get_one_line_group(line)
        except NoGroupFoundError:
            pass
        else:
            return key, value

        # evals if key:value is defined in multiple lines
        mll_group = _list_multiliner_re.match(line)

        if mll_group:
            key = mll_group[1]
            idx = line.find('[')
            # need to send `fin` and not `pure_lines`
            block = line[idx:] + _get_list_block(fin)
            list_ = _eval_list_str(block)
            return key, list_

        # if the flow reaches here...
        raise ValueError(f'Can\'t process this line: {line!r}')


def _is_correct_line(line):
    return not _is_comment(line)


def _is_comment(line):
    """
    Assert if line is a comment or a line to ignore.

    Expects striped lines.
    """
    is_comment = \
        line.startswith('#') \
        or not bool(line)

    return is_comment


def _process_line(line):
    return line.strip()


def _replace_bool(s):
    "Replace booleans and try to parse."""
    s = s.replace('true', 'True')
    s = s.replace('false', 'False')
    return s


def _eval_list_str(s):
    """
    Evaluates a string to a list.

    List string must be already enclosed in brackets `[]`.
    """
    s = _replace_bool(s)
    return ast.literal_eval(s)


def _get_one_line_group(line):
    """Attempt to identify a key:value pair in a single line."""
    for method, func in regex_single_line_methods:
        group = method.match(line)
        if group:
            return group[1], func(group[2])  # return first found

    # all regex-based methods have failed
    # now try methods not based on regex
    try:
        key, value = (s.strip() for s in line.split('='))
    except ValueError:
        raise NoGroupFoundError('Line does not match any group.')

    for method in regex_single_line_special_methods:
        try:
            return key, method(value)
        except Exception:
            continue

    raise NoGroupFoundError('Line does not match any group.')


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


def _remove_trailing_zeros_in_headers(d):
    """The suffix '.0' is removed from keys."""
    d1 = {}
    for key, value in d.items():
        if isinstance(value, dict) and key.endswith('.0'):
            d1[key[:-2]] = value
        else:
            d1[key] = value
    return d1


def get_module_name(name):
    """Gets the name according to the config parser."""
    return name.split('.')[0]


# methods to parse single line values
# (regex, func)
regex_single_line_methods = [
    (_string_re, ast.literal_eval),
    (_number_re, ast.literal_eval),
    (_none_re, lambda x: None),
    (_list_one_liner_re, _eval_list_str),
    (_true_re, lambda x: True),
    (_false_re, lambda x: False),
    ]

regex_single_line_special_methods = [
    datetime.fromisoformat,
    ]
