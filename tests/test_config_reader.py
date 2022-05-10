"""Test config reader."""
import os
from datetime import datetime
from pathlib import Path, PosixPath, WindowsPath

import pytest

from haddock import EmptyPath
from haddock.gear import config_reader


def on_windows():
    """Check we if are on windows."""
    try:
        WindowsPath()
    except NotImplementedError:
        return False
    else:
        return True


def on_linux():
    """Check we if are on windows."""
    try:
        PosixPath()
    except NotImplementedError:
        return False
    else:
        return True


@pytest.mark.parametrize(
    'line,expected',
    [
        ('[header]', 'header'),
        ('[under_scores_work]', 'under_scores_work'),
        ],
    )
def test_main_header_re(line, expected):
    """Test header regex."""
    result = config_reader._main_header_re.match(line)
    assert result[1] == expected


@pytest.mark.parametrize(
    'line,expected',
    [
        ('[another.header]', 'another.header'),
        ('[another.header_2]', 'another.header_2'),
        ],
    )
def test_sub_header_re(line, expected):
    """Test header regex."""
    result = config_reader._sub_header_re.match(line)
    assert result[1] == expected


@pytest.mark.parametrize(
    'line',
    [
        '[[header]]',
        '[héàder]',
        'value = "some_string"',
        '[header with spaces]',
        '[not.valid]',
        ],
    )
def test_main_header_re_wrong(line):
    """Test main header wrong."""
    assert config_reader._main_header_re.match(line) is None


@pytest.mark.parametrize(
    'line',
    [
        '[[header]]',
        'value = "some_string"',
        '[header with spaces]',
        '[not.valid with spaces]',
        '[single]',
        ],
    )
def test_sub_header_re_wrong(line):
    """Test sub header wrong."""
    assert config_reader._sub_header_re.match(line) is None


@pytest.mark.parametrize(
    'line,name,value',
    [
        ('value = "some_string"', 'value', '"some_string"'),
        ("value = 'some'", 'value', "'some'"),
        (
            "value = 'someSOME123790'",
            'value',
            "'someSOME123790'",
            ),
        ("value = ''", 'value', "''"),
        ('value = ""', 'value', '""'),
        ("var2='other'", 'var2', "'other'"),
        ("var2='other'#somecomment", 'var2', "'other'"),
        ("var2='other'    #somecomment", 'var2', "'other'"),
        ("var_2='other'    #somecomment", 'var_2', "'other'"),
        ("v_ar_2='other'    #somecomment", 'v_ar_2', "'other'"),
        ],
    )
def test_string_re(line, name, value):
    """Test string regex."""
    result = config_reader._string_re.match(line)
    assert result[1] == name
    assert result[2] == value


@pytest.mark.parametrize(
    'line',
    [
        "value='s*ràngë'",
        "vàlüé = 'some'",
        "value='.'",
        'value=1',
        'value=other',
        'value=true',
        'value=string with spaces',
        'value = "string_with-other/chars"',
        'value = ["list"]',
        ],
    )
def test_string_re_wrong(line):
    """Test string wrong examples."""
    assert config_reader._string_re.match(line) is None


@pytest.mark.skipif(on_windows(), reason="We are on Windows")
@pytest.mark.parametrize(
    'line,name,value_fname',
    [
        (r'value_fname = "file"', 'value_fname', Path('file')),
        (
            r'molecules = "../../data/2oob.pdb"',
            'molecules',
            Path('../../data/2oob.pdb'),
            ),
        (
            'value_fname = "file_that_does_not_exist"#with comments',
            "value_fname",
            Path("file_that_does_not_exist"),
            ),
        (
            "value_fname = './path/to/file'",
            "value_fname",
            Path('path', 'to', 'file')),
        (
            "value_fname = '../path/to/file'",
            "value_fname",
            Path('..', 'path', 'to', 'file')),
        ],
    )
def test_File_re_linux(line, name, value_fname):
    """Test string regex."""
    path_re = config_reader._File_re()
    result = path_re.match(line)
    assert result[1] == name
    assert result[2] == value_fname


@pytest.mark.skipif(on_linux(), reason="We are on Linux")
@pytest.mark.parametrize(
    'line,name,value_fname',
    [
        (r'value_fname = "file"', 'value_fname', Path('file')),
        (
            r'value_fname = ".\file_that_does_not_exist"#with comments',
            "value_fname",
            Path("file_that_does_not_exist"),
            ),
        (
            r"value_fname = '.\path\to\file'",
            "value_fname",
            Path('path', 'to', 'file')),
        (
            r"value_fname = '..\path\to\file'",
            "value_fname",
            Path('..', 'path', 'to', 'file')),
        ],
    )
def test_File_re_windows(line, name, value_fname):
    """Test string regex."""
    path_re = config_reader._File_re()
    result = path_re.match(line)
    assert result[1] == name
    assert result[2] == value_fname


@pytest.mark.parametrize(
    'line',
    [
        "value = 4",
        "value=''",
        "value='some_file'",
        "value='./path with spaces/file'",
        "value='./path.with.dots/file'",
        r"value='.\path with spaces\file'",
        r"value='./path/with/strângë/ch*rs/",
        "value_fname = 4",
        "vàlüé_fname = 'somepath'",
        "value_fname=''",
        "value_fname='./'",
        "value_fname='./path with spaces/file'",
        "value_fname='./path.with.dots/file'",
        r"value_fname='.\path with spaces\file'",
        r"value_fname='./path/with/strângë/ch*rs/",
        ],
    )
def test_File_re_wrong(line):
    """Test string regex."""
    path_re = config_reader._File_re()
    result = path_re.match(line)
    assert result is None


@pytest.mark.parametrize(
    "string, name",
    [
        ("value_fname = ''", "value_fname"),
        ('value_fname = ""', "value_fname"),
        ]
    )
def test_EmptyPath_re(string, name):
    """Test empty path."""
    emptypath_re = config_reader._EmptyFilePath_re()
    result = emptypath_re.match(string)
    assert result[1] == name
    assert isinstance(result[2], EmptyPath)


@pytest.mark.parametrize(
    "string",
    [
        "value_fname = 'something'",
        'value_fname = "./path/file.ext"',
        'value = 4',
        'vàlue = ""',
        ]
    )
def test_EmptyPath_re_wrong(string):
    """Test empty path."""
    emptypath_re = config_reader._EmptyFilePath_re()
    result = emptypath_re.match(string)
    assert result is None


@pytest.mark.parametrize(
    'line,name,value',
    [
        ('value = null', 'value', 'null'),
        ("value = Null", 'value', 'Null'),
        ("var2=None", 'var2', 'None'),
        ("var2=none#somecomment", 'var2', 'none'),
        ("var2=None    #somecomment", 'var2', 'None'),
        ],
    )
def test_none_re(line, name, value):
    """Test none regex."""
    result = config_reader._none_re.match(line)
    assert result[1] == name
    assert result[2] == value


@pytest.mark.parametrize(
    'line',
    [
        'value=1',
        'valué = None',
        'value=other',
        'value=true',
        'value="some_string"',
        'value = ["list"]',
        ],
    )
def test_none_re_wrong(line):
    """Test None wrong examples."""
    assert config_reader._none_re.match(line) is None


@pytest.mark.parametrize(
    'line,name,value',
    [
        ('value = nan', 'value', 'nan'),
        ('value = nAn', 'value', 'nAn'),
        ('value = naN', 'value', 'naN'),
        ('value = NaN', 'value', 'NaN'),
        ("var2=nan", 'var2', 'nan'),
        ("var2=NaN#somecomment", 'var2', 'NaN'),
        ("var2=NaN    #somecomment", 'var2', 'NaN'),
        ],
    )
def test_nan_re(line, name, value):
    """Test none regex."""
    result = config_reader._nan_re.match(line)
    assert result[1] == name
    assert result[2] == value


@pytest.mark.parametrize(
    'line,number',
    [
        ('value = 00', "00"),
        ('value = 00\n', "00"),
        ('value = 0', "0"),
        ('value = -0', "-0"),
        ('value = 0.', "0."),
        ('value = .0', ".0"),
        ('value = 12.3', "12.3"),
        ('w_vdw_0 = 0.01', "0.01"),
        ('value = 12.34', "12.34"),
        ('value = -12.34', "-12.34"),
        ('value = 12.34 # with comment', "12.34"),
        ('value = 10E40', "10E40"),
        ('value = 1E4', "1E4"),
        ('value = 1e4', "1e4"),
        ('value = 1E-4', "1E-4"),
        ('value = -10E40', "-10E40"),
        ('value = -10E-40', "-10E-40"),
        ('value = -.1E-4', "-.1E-4"),
        ('value = 10.2E30', "10.2E30"),
        ('value = 10.2e30', "10.2e30"),
        ('value = .10E30', ".10E30"),
        ('value = -10.2E30', "-10.2E30"),
        ('value = -10.2e30', "-10.2e30"),
        ('value = -.10E30', "-.10E30"),
        ('value = -.10E-30', "-.10E-30"),
        ('value = 123#with comment', '123'),
        ('value = -.10E-30#withcomment', "-.10E-30"),
        ('value = -.10e-30#withcomment', "-.10e-30"),
        ('value = -.10#E-30', "-.10"),
        ('va_lue_0 = -.10#E-30', "-.10"),
        ('va_lue_0 = -.10#e-30', "-.10"),
        ],
    )
def test_number_re(line, number):
    """Test number regex."""
    result = config_reader._number_re.match(line)
    assert result[1] == line.split('=')[0].strip()
    assert result[2] == number


@pytest.mark.parametrize(
    'line',
    [
        "válúe = 12",
        "value = 12.34wrong",
        "value = 12.34.12",
        "value = 1E4.4",
        "value = 1e4.4",
        "value = .10E30E",
        "value = 10.2.E30",
        "value = 10.2.e30",
        "value = E10",
        "value = e10",
        "value = .E10",
        "value = -e",
        "value = -E10",
        "value = 10-10E19",
        "value = 10-10e19",
        ],
    )
def test_number_re_wrong(line):
    """Test number regex wrong examples."""
    assert config_reader._number_re.match(line) is None


@pytest.mark.parametrize(
    'line,name,value',
    [
        ("ports = [ 8000, 8001, 8002 ]", "ports", "[ 8000, 8001, 8002 ]"),
        (
            '_data = [ ["gamma", "delta"], [1, 2] ] # valid comment',
            "_data",
            '[ ["gamma", "delta"], [1, 2] ]',
            ),
        (
            '_data = [ ["gamma", "delta"], [1, 2] ]# valid comment\n',
            "_data",
            '[ ["gamma", "delta"], [1, 2] ]',
            ),
        ("ports = [8000]]]", "ports", "[8000]]]"),
        ("ports=[8000]", "ports", "[8000]"),
        # wrong bracket closing are accepted at this stage,
        # It will then raise error when evaluating
        ],
    )
def test_list_one_liner_re(line, name, value):
    """Test regex capturing lists defined in one line."""
    result = config_reader._list_one_liner_re.match(line)
    assert result[1] == name
    assert result[2] == value


@pytest.mark.parametrize(
    'line',
    [
        'value = ][8000]',
        'valué = [8000]',
        'value = 8000',
        'value = [8000',
        'value ="somestring',
        ],
    )
def test_list_one_liner_re_wrong(line):
    """Test regex capturing lists defined in one line."""
    assert config_reader._list_one_liner_re.match(line) is None


@pytest.mark.parametrize(
    'line',
    [
        'value=[',
        'value = [ ## some comment',
        'value = [\n',
        ],
    )
def test_list_multi_liner(line):
    """Test regex captures the first line of a multiline list."""
    assert config_reader._list_multiliner_re.match(line)


@pytest.mark.parametrize(
    'line',
    [
        'value=',
        'valué=[',
        'value = "some" ## some comment',
        'value = [1]',
        ],
    )
def test_list_multi_liner_wrong(line):
    """Test regex captures the first line of a multiline list."""
    assert config_reader._list_multiliner_re.match(line) is None


@pytest.mark.parametrize(
    'line,name,value',
    [
        ('value = 1', 'value', 1),
        ('value = "some"', 'value', "some"),
        ('list = [12, 13]', 'list', [12, 13]),
        ('path_fname = "./here"', 'path_fname', Path("./here")),
        (
            'date = 1979-05-27T07:32:00-08:00',
            'date',
            datetime.fromisoformat('1979-05-27T07:32:00-08:00'),
            ),
        ],
    )
def test_get_one_line_group(line, name, value):
    """Test get one line values."""
    n, v = config_reader._get_one_line_group(line)
    assert n == name
    assert v == value


@pytest.mark.parametrize(
    'lines,expected',
    [
        (
            ('"some value", 15\n', ']\n', '\n'),
            ('"some value", 15]')
            ),
        ],
    )
def test_get_block_list(lines, expected):
    """Test captures multiline list."""
    r = config_reader._get_list_block(lines)
    assert r == expected


@pytest.mark.parametrize(
    'header,name',
    [
        ('header', 'header'),
        ('header.1', 'header'),
        ('header.1.1..1', 'header'),
        ],
    )
def test_get_module_name(header, name):
    """Test get module name."""
    assert config_reader.get_module_name(header) == name


@pytest.mark.parametrize(
    's,expected',
    [
        ("[12, 13]", [12, 13]),
        (
            "[true, True, false, False, 'all']",
            [True, True, False, False, "all"]
            ),
        ],
    )
def test_eval_list_string(s, expected):
    """Test conver string to list."""
    r = config_reader._eval_list_str(s)
    assert r == expected


@pytest.mark.parametrize(
    'line,expected',
    [
        ('# this is a comment', True),
        ('value = 123', False),
        ('', True),
        ('    # other comment', False),  # expects striped lines
        ],
    )
def test_is_comment(line, expected):
    """Test if line if comment."""
    r = config_reader._is_comment(line)
    assert r == expected


@pytest.mark.parametrize(
    'line,expected',
    [
        ('# this is a comment', False),
        ('value = 123', True),
        ('', False),
        ('    # other comment', True),  # expects striped lines
        ],
    )
def test_is_correct_line(line, expected):
    """Test if line if comment."""
    r = config_reader._is_correct_line(line)
    assert r == expected


@pytest.mark.parametrize(
    'line,expected',
    [
        (' ', ''),
        ('  value = 123', 'value = 123'),
        ('', ''),
        ('    # other      \n', '# other'),
        ],
    )
def test_process_line(line, expected):
    """Test processes lines."""
    r = config_reader._process_line(line)
    assert r == expected


_config_example_1 = """
# some comment
num1 = 10
name="somestring"
null_value = None#some comment
w_vdw_1 = 1.0
w_vdw_0 = 0.01
w_vdw_2 = 1.0
[headerone]
name = "theotherstring"
_list = [
    12,
    "foo",
    [56, 86],
    ]

[headerone.weights]
val = 1
bsa = 0.1
elec = -1
list1 = [1, 2,"3"]

[headerone] #some comment
[headerone.weights]
other = 50
"""

_config_example_dict_1 = {
    "num1": 10,
    "name": "somestring",
    "null_value": None,
    "w_vdw_0": 0.01,
    "w_vdw_1": 1.0,
    "w_vdw_2": 1.0,
    "headerone": {
        "name": "theotherstring",
        "_list": [12, "foo", [56, 86]],
        "weights": {
            "val": 1,
            "bsa": 0.1,
            "elec": -1,
            "list1": [1, 2, "3"],
            },
        },
    "headerone.1": {'weights': {'other': 50}},
    }

_config_example_dict_2 = {
    "num1": 10,
    "some_path_fname": Path("pointing", "to", "some", "path"),
    "module": {
        "name": ["../some/file", "../some/otherfile"],
        "d1": {
            "var1": 1,
            "var2": None,
            "d2": {
                "var3": True,
                "list_": [1, 2, 3, 4],
                },
            },
        },
    }

_config_example_2 = """num1 = 10
some_path_fname = "./pointing/to/some/path"
[module]
name = ["../some/file", "../some/otherfile"]
[module.d1]
var1 = 1
var2 = None
[module.d1.d2]
var3 = True
list_ = [
    1,
    2,
    3, 4
    ]
"""

# this examples shows the behaviour of subkey repetition
_config_example_3 = """
val = 1
[header]
[header.d1]
val2 = 10
[header.d1]
val3 = 20
"""

_config_example_dict_3 = {
    "val": 1,
    "header": {
        "d1": {
            "val2": 10,
            "val3": 20,
            }
        }
    }


@pytest.mark.parametrize(
    'config,expected',
    [
        (_config_example_1, _config_example_dict_1),
        (_config_example_2, _config_example_dict_2),
        (_config_example_3, _config_example_dict_3),
        ],
    )
def test_read_config(config, expected):
    """Test read config."""
    r = config_reader._read_config((i for i in config.split(os.linesep)))
    assert r == expected


_config_broken_1 = """
name1 = "good"
[header.subdict]
number = 15
"""

_config_broken_2 = """
name1 = "good"
[module]
s = None
[header.subdict]  # this should give an error
number = 15
"""

_config_broken_3 = """
name1 = "good"
name1 = "bad"
"""

_config_broken_4 = """
val = 1
[header]
[header.d1]
val2 = 10
[header.d1]
val2 = 20
"""

_config_broken_5 = """
broken_list = [
    1,
    2,3
    ]
name = 1
"""

_config_broken_6 = """
broken_list = [ 1,
    2,3]
name = 1
"""

_config_broken_7 = """
broken_list = [
    1,
    2,3,
    strängé,
    +
    ]

"""

_config_broken_8 = """
nämé1 = "good"
"""


@pytest.mark.parametrize(
    'config, error',
    [
        (_config_broken_1, config_reader.ConfigFormatError),
        (_config_broken_2, config_reader.ConfigFormatError),
        (_config_broken_3, config_reader.DuplicatedParameterError),
        (_config_broken_4, config_reader.DuplicatedParameterError),
        (_config_broken_5, config_reader.MultilineListDefinitionError),
        (_config_broken_6, config_reader.MultilineListDefinitionError),
        (_config_broken_7, config_reader.MultilineListDefinitionError),
        (_config_broken_8, config_reader.NoGroupFoundError),
        ],
    )
def test_config_format_errors(config, error):
    """Test broken configs and errors."""
    with pytest.raises(error):
        config_reader._read_config((i for i in config.split(os.linesep)))
