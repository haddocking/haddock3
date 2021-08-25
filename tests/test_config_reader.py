"""Test config reader."""
from datetime import datetime


import pytest


from haddock.gear import config_reader


@pytest.mark.parametrize(
    'line,expected',
    [
        ('[header]', 'header'),
        ('[another.header]', 'another.header'),
        ],
    )
def test_header_re(line, expected):
    """Test header regex."""
    result = config_reader._header_re.match(line)
    assert result[1] == expected


@pytest.mark.parametrize(
    'line',
    [
        '[[header]]',
        'value = "some string"',
        ],
    )
def test_header_re_wrong(line):
    assert config_reader._header_re.match(line) is None


@pytest.mark.parametrize(
    'line,name,value',
    [
        ('value = "some string"', 'value', '"some string"'),
        ("value = 'some'", 'value', "'some'"),
        ("var2='other'", 'var2', "'other'"),
        ("var2='other'#somecomment", 'var2', "'other'"),
        ("var2='other'    #somecomment", 'var2', "'other'"),
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
        'value=1',
        'value=other',
        'value=true',
        'value = ["list"]',
        ],
    )
def test_string_re_wrong(line):
    assert config_reader._string_re.match(line) is None


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
        ('value = 12.34', "12.34"),
        ('value = -12.34', "-12.34"),
        ('value = 12.34 # with comment', "12.34"),
        ('value = 10E40', "10E40"),
        ('value = 1E4', "1E4"),
        ('value = 1E-4', "1E-4"),
        ('value = -10E40', "-10E40"),
        ('value = -10E-40', "-10E-40"),
        ('value = -.1E-4', "-.1E-4"),
        ('value = 10.2E30', "10.2E30"),
        ('value = .10E30', ".10E30"),
        ('value = -10.2E30', "-10.2E30"),
        ('value = -.10E30', "-.10E30"),
        ('value = -.10E-30', "-.10E-30"),
        ],
    )
def test_number_re(line, number):
    result = config_reader._number_re.match(line)
    assert result[1] == "value"
    assert result[2] == number


@pytest.mark.parametrize(
    'line',
    [
        "value = 12.34wrong",
        "value = 12.34.12",
        "value = 1E4.4",
        "value = .10E30E",
        "value = 10.2.E30",
        "value = 123#with wrong comment",
        "value = E10",
        "value = .E10",
        "value = -e",
        "value = -E10",
        "value = 10-10E19",
        ],
    )
def test_number_re_wront(line):
    assert config_reader._number_re.match(line) is None


@pytest.mark.parametrize(
    'line,name,value',
    [
        ("ports = [ 8000, 8001, 8002 ]", "ports", "[ 8000, 8001, 8002 ]"),
        ('_data = [ ["gamma", "delta"], [1, 2] ] # valid comment', "_data", '[ ["gamma", "delta"], [1, 2] ]'),
        ('_data = [ ["gamma", "delta"], [1, 2] ]# valid comment\n', "_data", '[ ["gamma", "delta"], [1, 2] ]'),
        ("ports = [8000]]]", "ports", "[8000]]]"),
        ("ports=[8000]", "ports", "[8000]"),
        # wrong bracket formation are accepted. It will then raise error when evaluating
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
    "Test regex captures the first line of a multiline list."""
    assert config_reader._list_multiliner_re.match(line)


@pytest.mark.parametrize(
    'line',
    [
        'value=',
        'value = "some" ## some comment',
        'value = [1]',
        ],
    )
def test_list_multi_liner_wrong(line):
    "Test regex captures the first line of a multiline list."""
    assert config_reader._list_multiliner_re.match(line) is None


@pytest.mark.parametrize(
    'line,name,value',
    [
        ('value = 1', 'value', 1),
        ('value = "some"', 'value', "some"),
        ('list = [12, 13]', 'list', [12, 13]),
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
    'd,expected',
    [
        (
            {"one.two": {'foo': 'bar'}},
            {"one": {"two": {'foo': 'bar'}}},
            ),
        (
            {"one": {'foo': 'bar'}},
            {"one": {'foo': 'bar'}},
            ),
        (
            {"one": {'foo': 'bar'}, "two": 2},
            {"one": {'foo': 'bar'}, "two": 2},
            ),
        (
            {
                "one.two": {'foo': 'bar'},
                "one.foo": {'zoo': 'animal'},
                "other": 15,
                },
            {
                "one": {
                    "two": {'foo': 'bar'},
                    "foo": {'zoo': 'animal'},
                    },
                "other": 15,
                },
            ),
        (
            {
                "one.two": {'foo': 'bar'},
                "one.two.loo": {'zoo': 'animal'},
                "other": 15,
                },
            {
                "one": {
                    "two": {
                        'foo': 'bar',
                        "loo": {"zoo": "animal"}
                        },
                    },
                "other": 15,
                },
            ),
        ],
    )
def test_make_nested_keys(d, expected):
    """Test making nested keys from numbered headers."""
    r = config_reader._make_nested_keys(d)
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
    'lines,expected',
    [
        (
            """
            # some comment
            num1 = 10
            name="some string"
            [headerone]
            name = "the other string"
            _list = [
                12,
                "foo",
                [56, 86],
                ]

            [headerone]
            """.split('\n'),
            {
                "num1": 10,
                "name": "some string",
                "headerone": {
                    "name": "the other string",
                    "_list": [12, "foo", [56, 86]],
                    },
                "headerone.1": {},
                },
            ),
        ],
    )
def test_read_config(lines,expected):
    """Test read config."""
    print(lines)
    r = config_reader._read_config((i for i in lines))
    assert r == expected


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
        ('    # other comment', False), # expects striped lines
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
        ('    # other comment', True), # expects striped lines
        ],
    )
def test_is_comment(line, expected):
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
