"""
Implementation of a TOML-like configuration file for HADDOCK3.

HADDOCK3 user configuration files follow TOML syntax plus additional
features required for HADDOCK3. Therefore, we have implemented our own
TOML-like parser to accommodate the extra features needed for HADDOCK3.

The most relevant features of HADDOCK3 user configuration files are:

**Accepts repeated keys:** blocks can have repeated names, for example:

.. code:: toml

        [topoaa]
        # parameters here...

        [rigidbody]
        # parameters here...

        [caprieval]
        # parameters here...

        [flexref]
        # parameters here...

        [caprieval]
        # parameters here...

**Allows in-line comments:**

.. code:: toml

    [module]
    parameter = value # some comment

**The following types are implemented**:

* strings
* numbers
* null/none
* ``nan``
* lists defined in one lines
* lists defined in multiple lines (require empty line after the list
  definition)
* `datetime.fromisoformat`

**To efficiently use this module, see functions:**

* :py:func:`haddock.gear.config_reader.read_config`
* :py:func:`haddock.gear.config_reader.get_module_name`
"""
import ast
import re
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path, PosixPath, WindowsPath

from haddock import EmptyPath
from haddock.core.defaults import RUNDIR
from haddock.core.exceptions import ConfigurationError
from haddock.libs.libfunc import false, give_same, nan, none, true


# The main value types are parsed from the configuration file using
# regular expressions. There is a regular expression for each type of
# parameter: str, numbers, null, etc. There are also regular expressions
# to capture lists. For each regex there is a link at regex101 that
# explains what is captured and what is not.
#
# Also, these regular expressions avoid introduction of strange
# characters that would be useless for the HADDOCK3 or may create
# problems.
#
# At the end of the file there is a list that maps each regular
# expression with the function that reads the captured value from the
# configuration file.
#
# the re.ASCII parameter makes sure non-ascii chars are rejected with \w key

# Captures the main headers.
# https://regex101.com/r/9urqti/1
_main_header_re = re.compile(r'^ *\[(\w+)\]', re.ASCII)

# Captures sub-headers
# https://regex101.com/r/6OpJJ8/1
# thanks https://stackoverflow.com/questions/39158902
_sub_header_re = re.compile(r'^ *\[(\w+(?:\.\w+)+)\]', re.ASCII)

# capture string parameters (paths are not strings, see later)
# https://regex101.com/r/0LYCAG/1
_string_re = re.compile(
    r'''^ *(\w+) *= *'''
    r'''("([\w\-]*?)"|'([\w\-]*?)')''',
    re.ASCII,
    )

# captures numbers
# https://regex101.com/r/Ur6TUK/1
_number_re = re.compile(
    r'^ *(\w+) *= *'
    r'(\-?\d+\.?\d*|\-?\.\d+|\-?\.?\d+[eE]\-?\d+|-?\d+\.?\d*[eE]\d+)(?: |#|$)',
    re.ASCII,
    )

# captures None or Null
# https://regex101.com/r/izatQ3/1
_none_re = re.compile(r'^ *(\w+) *= *([Nn]one|[Nn]ull)', re.ASCII)

# captures Nan
# https://regex101.com/r/3OVbqs/1
_nan_re = re.compile(r'^ *(\w+) *= *([nN][aA][nN])', re.ASCII)

# captures a list that is defined in a single line
# https://regex101.com/r/XtgU8I/1
_list_one_liner_re = re.compile(
    r'''^ *(\w+) *= *(\[[\w\ \,\-\"\'\[\]\.\\\/]*\])''',
    re.ASCII,
    )

# captures the beginning of a list defined in multiple lines. The rest
# of the list will be captured with a function later on.
# https://regex101.com/r/eakIzB/1
_list_multiliner_re = re.compile(
    r'''^ *(\w+) *= *\[(?:[\w\ \,\-\"\'\[\]\.\\\/]*'''
    r'''[\w\ \,\-\"\'\[\.\\\/](?:#+[\w\ ]*)?)?$''',
    re.ASCII,
    )

# captures the inside of a multiline list to make sure no strange chars
# are passed in. This is read within a special function so, there's no
# need for a regex101 link here.
_list_multiline_content_re = re.compile(
    r'''^ *([\w\ \,\-\"\'\[\]\.\\\/]*)''',
    re.ASCII,
    )

# captures true or false
# https://regex101.com/r/mkImqk/1
_true_re = re.compile(r'^ *(\w+) *= *([tT]rue)', re.ASCII)
# https://regex101.com/r/duMPt3/1
_false_re = re.compile(r'^ *(\w+) *= *([fF]alse)', re.ASCII)


# Some parameters are paths. The configuration reader reads those as
# pathlib.Paths so that they are injected properly in the rest of the
# workflow. In general, any parameter ending with `_fname` is a path,
# but there are also other parameters that are paths. Those need to be
# added to this list bellow:
_keys_that_accept_files = [
    "cns_exec",
    "executable",
    "molecules",
    RUNDIR,
    ]

# we join the list and inject it in the regex below
_keys_files_regex = "|".join(_keys_that_accept_files)

# regex accepting linux paths (PosixPaths)
# https://regex101.com/r/UQelVX/1
_file_linux_re = re.compile(
    r'''^ *(\w+_fname|'''
    f'''{_keys_files_regex}'''
    r''') *= *'''
    r'''("((\.{1,2}\/)*[\w\-\/]+[\w\-]\.?[\w\-]+)"'''
    r'''|'((\.{1,2}\/)*[\w\-\/]+[\w\-]\.?[\w\-]+)')''',
    re.ASCII,
    )

# regex accepting Windows paths (WindowsPaths)
# same as linux but with \.
# at the time of wrinting HADDOCK3 did not supported windows.
_file_windows_re = re.compile(
    r'''^ *(\w+_fname|'''
    f'''{_keys_files_regex}'''
    r''') *= *'''
    r'''("((\.{1,2}\\)*[\w\-\\]+[\w\-]\.?[\w\-]+)"'''
    r'''|'((\.{1,2}\\)*[\w\-\\]+[\w\-]\.?[\w\-]+)')''',
    re.ASCII,
    )

# Empty paths are defined path parameters with empty strings "".
# https://regex101.com/r/jtGbT2/1
_emptypath_re = re.compile(
    r'''^ *(\w+_fname|'''
    f'''{_keys_files_regex}'''
    r''') *= *(""|'')''',
    re.ASCII,
    )


# re.compile objects have a `.match` method. This _File_re class
# implements the logic to identify paths with a `.match` method so that
# all regex match methods can be chained, see below.
class _File_re:
    def __init__(self):
        """Map to file path regex."""
        self.re_groups = (
            (_file_linux_re, PosixPath),
            (_file_windows_re, WindowsPath),
            )

    def match(self, string):
        """
        Polymorphs re.match.

        Return the match group if the file exists.

        Else, return false.

        Made specifically to :func:`_get_one_line_group`.
        """
        for regex, PathType in self.re_groups:
            group = regex.match(string)
            if group:
                try:
                    PathType()
                except NotImplementedError:
                    emsg = (
                        "The path type you defined is not compatible with "
                        "your system. Maybe you define a Windows path and you "
                        "are running HADDOCK3 on Linux."
                        )
                    raise ConfigurationError(emsg) from NotImplementedError

                # strip the quotes so we don't need to use
                # ast.literal_eval
                value = group[2].strip("'\"")
                if value:
                    p = Path(value)
                    return group[0], group[1], p

        # everything else returns None by definition
        return None  # we want to return None to match the re.match()


# same as for _File_re but simulates _EmptyFilePath_re.
class _EmptyFilePath_re:
    def __init__(self):
        """Map to EmptyPath regex."""
        self.re = _emptypath_re

    def match(self, string):
        """
        Polymorphs re.match.

        Return the match group if the file exists.

        Else, return false.

        Made specifically to :func:`_get_one_line_group`.
        """
        group = self.re.match(string)
        if group:
            return group[0], group[1], EmptyPath()

        # everything else returns None by definition
        return None  # we want to return None to match the re.match()


# defines some specific exceptions
class NoGroupFoundError(Exception):
    """
    Exception if no group is found.

    Used for single line.
    """

    pass


class ConfigFormatError(Exception):
    """Exception if there is a format error."""

    pass


class DuplicatedParameterError(Exception):
    """Exception if duplicated parameters are found."""

    pass


class MultilineListDefinitionError(Exception):
    """Exception when a multiline list is wrongly formatted."""

    pass


# main public API
def read_config(fpath):
    """
    Read HADDOCK3 configure file to a dictionary.

    Parameters
    ----------
    fpath : str or :external:py:class:`pathlib.Path`
        Path to user configuration file.

    Returns
    -------
    dictionary
        Representing the user configuration file where first level of
        keys are the module names. Repeated modules will have a numeric
        suffix, for example: ``module.1``.
    """
    with open(fpath, 'r') as fin:
        return _read_config(fin)


def _read_config(fin):
    """
    Read the config line by line.

    Parameters
    ----------
    fin : a openned file content generator
        A generator expressing the content of a file in lines. Cannot be
        a list.

    Returns
    -------
    dictionary
        Representing the user configuration file where first level of
        keys are the module names. Repeated modules will have a numeric
        suffix, for example: ``module.1``.
    """
    # main dictionary to store the configuration
    d = {}

    # main keys always come before any subdictionary
    d1 = d

    pure_lines = filter(_is_correct_line, map(_process_line, fin))
    header = None
    # this loops only over the lines that have meaningful information
    for line in pure_lines:

        # collects header and subheader
        main_header_group = _main_header_re.match(line)
        sub_header_group = _sub_header_re.match(line)

        # prepares header
        if main_header_group:
            header = _update_key_number(main_header_group[1], d)
            d1 = d.setdefault(header, {})

        # prepares subheader
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

        # is a value line
        else:
            # reads the parameter line
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
    """
    Update key number pairs.

    Despite the original config can accept repated header keys,
    a python dictionary cannot. So a integer suffix is added.
    """
    _num = str(sum(1 for k in d if k.startswith(key)) + offset)
    return key + sep + _num


def _read_value(line, fin):
    """Read parameter:value pairs."""
    # attempts to read a parameter defined in a single line
    try:
        key, value = _get_one_line_group(line)
    except NoGroupFoundError:
        # if the group is not found, try the multiline list
        pass
    else:
        return key, value

    # evals if key:value is defined in multiple lines
    mll_group = _list_multiliner_re.match(line)

    # tries to read the multiline list.
    if mll_group:

        # define the error message to be used later
        emsg = (
            f"Can't process this line {line!r}. "
            "The multiline list is not properly formatted."
            )

        key = mll_group[1]
        idx = line.find('[')
        # need to send `fin` and not `pure_lines`
        with _multiline_block_error(emsg, ValueError):
            # defines a block to pass to _eval_list_str
            block = line[idx:] + _get_list_block(fin)

        with _multiline_block_error(emsg, SyntaxError):
            list_ = _eval_list_str(block)

        return key, list_

        # if the flow reaches here...
        raise MultilineListDefinitionError(emsg)

    raise NoGroupFoundError(f"Could not read line: {line!r}.")


@contextmanager
def _multiline_block_error(emsg, exception):
    """Error context for multiline lists."""
    try:
        yield
    except exception as err:
        raise MultilineListDefinitionError(emsg) from err


def _is_correct_line(line):
    """
    Define what is a correct line.

    Can be expanded in the future if needed.
    """
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
    """Replace booleans and before passing to ast.literal_eval."""
    s = s.replace('true', 'True')
    s = s.replace('false', 'False')
    return s


def _eval_list_str(s):
    """
    Evaluate a string to a list.

    List string must be already enclosed in brackets `[]`.
    """
    s = _replace_bool(s)
    return ast.literal_eval(s)


def _get_one_line_group(line):
    """Attempt to identify a key:value pair in a single line."""
    # attempts to read `key = value` pairs in a line.
    # the order by which the reading is performed is important
    # that is, the order of `regex_single_line_methods`.
    for method, func in regex_single_line_methods:
        group = method.match(line)
        if group:
            # return the first found
            return group[1], func(group[2])

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

    raise NoGroupFoundError(f"Line does not match any group: {line!r}")


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

        group = _list_multiline_content_re.match(line)
        if group:
            block.append(line)
        else:
            raise ValueError("Bad group in multiline list.")

    return ''.join(block).strip()


def _remove_trailing_zeros_in_headers(d):
    """Remove the suffix '.0' from keys."""
    d1 = {}
    for key, value in d.items():
        if isinstance(value, dict) and key.endswith('.0'):
            d1[key[:-2]] = value
        else:
            d1[key] = value
    return d1


def get_module_name(name):
    """Get the name according to the config parser."""
    return name.split('.')[0]


# methods to parse single line values
# the order of this list matters
# the idea is simple: one regex compile object (or polymorphism) maps to
# a function that should process the captured value or return some
# predifined value accordingly.
regex_single_line_methods = (
    (_File_re(), give_same),
    (_EmptyFilePath_re(), give_same),
    (_string_re, ast.literal_eval),
    (_number_re, ast.literal_eval),
    (_nan_re, nan),
    (_none_re, none),
    (_list_one_liner_re, _eval_list_str),
    (_true_re, true),
    (_false_re, false),
    )

# methods to parse other values not defined by regex
# at the time of writing there are no parameters that require a date
# but, for a matter of demonstration, I added that functionality here
regex_single_line_special_methods = [
    datetime.fromisoformat,
    ]
