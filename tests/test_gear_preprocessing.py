"""Test preprocessing operations."""
import os
from itertools import zip_longest
from pathlib import Path

import pytest

from haddock.gear import preprocessing as pp

from . import broken_pdb, corrected_pdb, residues_top


def test_open_or_give():
    in1 = open(broken_pdb)
    input_ = [
        broken_pdb,  # path inside list
        str(broken_pdb),  # path inside list
        in1,
        broken_pdb.read_text().split(os.linesep),  # list
        tuple(broken_pdb.read_text().split(os.linesep)),  # tuple
        ]

    result = pp._open_or_give(input_)
    assert len(result) == 5
    for r in result:
        assert len(r) == 192
    in1.close()


@pytest.mark.parametrize(
    "value",
    [1, 1.1, {1: None}, None],
    )
def test_open_or_give_wrong(value):
    with pytest.raises(TypeError):
        pp._open_or_give(value)


def test_read_additional_residues():
    result = pp.read_additional_residues(residues_top)
    assert isinstance(result, tuple)
    assert result == ("DA2", "DE3", "DI", "DO", "DU1")


def test_read_additional_residues_lines():
    lines = Path(residues_top).read_text().split(os.linesep)
    result = pp.read_additional_residues.original(lines)
    assert isinstance(result, tuple)
    assert result == ("DA2", "DE3", "DI", "DO", "DU1")


def test_wrep_pdb_tidy():
    pp.wrep_pdb_tidy(pp._open_or_give([broken_pdb])[0], strict=False)


def test_wrep_pdb_tidy_strict():
    pp.wrep_pdb_tidy(pp._open_or_give([broken_pdb])[0], strict=True)


models_okay = """MODEL        1
ATOM      3  CA  ARG B   4      37.080  43.455  -3.421  1.00  0.00           C
ATOM      3  CA  GLU B   6      33.861  45.127  -2.233  1.00  0.00           C
ATOM      3  CA  ALA B   7      35.081  45.036   1.305  1.00  0.00           C
TER
ENDMDL
MODEL        2
ATOM      3  CA  ARG B   4      37.080  43.455  -3.421  1.00  0.00           C
ATOM      3  CA  GLU B   6      33.861  45.127  -2.233  1.00  0.00           C
ATOM      3  CA  ALA B   7      35.081  45.036   1.305  1.00  0.00           C
TER
ENDMDL
MODEL        3
ATOM      3  CA  ARG B   4      37.080  43.455  -3.421  1.00  0.00           C
ATOM      3  CA  GLU B   6      33.861  45.127  -2.233  1.00  0.00           C
ATOM      3  CA  ALA B   7      35.081  45.036   1.305  1.00  0.00           C
TER
ENDMDL
"""


models_wrong_2 = """MODEL        1
ATOM      3  CA  ARG B   4      37.080  43.455  -3.421  1.00  0.00           C
ATOM      3  CA  GLU B   6      33.861  45.127  -2.233  1.00  0.00           C
ATOM      3  CA  ALA B   7      35.081  45.036   1.305  1.00  0.00           C
TER
ENDMDL
MODEL        2
ATOM      3  CA  GLU B   6      33.861  45.127  -2.233  1.00  0.00           C
ATOM      3  CA  ALA B   7      35.081  45.036   1.305  1.00  0.00           C
TER
ENDMDL
MODEL        3
ATOM      3  CA  ARG B   4      37.080  43.455  -3.421  1.00  0.00           C
ATOM      3  CA  GLU B   6      33.861  45.127  -2.233  1.00  0.00           C
ATOM      3  CA  ALA B   7      35.081  45.036   1.305  1.00  0.00           C
TER
ENDMDL
"""

models_wrong_3 = """MODEL        1
ATOM      3  CA  ARG B   4      37.080  43.455  -3.421  1.00  0.00           C
ATOM      3  CA  GLU B   6      33.861  45.127  -2.233  1.00  0.00           C
ATOM      3  CA  ALA B   7      35.081  45.036   1.305  1.00  0.00           C
TER
ENDMDL
MODEL        2
ATOM      3  CA  ARG B   4      37.080  43.455  -3.421  1.00  0.00           C
ATOM      3  CA  GLU B   6      33.861  45.127  -2.233  1.00  0.00           C
ATOM      3  CA  ALA B   7      35.081  45.036   1.305  1.00  0.00           C
TER
ENDMDL
MODEL        3
ATOM      3  CA  GLU B   6      33.861  45.127  -2.233  1.00  0.00           C
ATOM      3  CA  ALA B   7      35.081  45.036   1.305  1.00  0.00           C
TER
ENDMDL
"""


def test_check_models():
    pp.models_should_have_the_same_labels(models_okay.split(os.linesep))


def test_check_models_wrong_2():
    _lines = models_wrong_2.split(os.linesep)
    with pytest.raises(pp.ModelsDifferError) as error:
        pp.models_should_have_the_same_labels(_lines)
    assert str(error.value) == "Labels in MODEL 2 differ from MODEL 1."


def test_check_models_wrong_3():
    _lines = models_wrong_3.split(os.linesep)
    with pytest.raises(pp.ModelsDifferError) as error:
        pp.models_should_have_the_same_labels(_lines)
    assert str(error.value) == "Labels in MODEL 3 differ from MODEL 1."


rep_chains_1 = [
    [
        'ATOM      3  CA  ARG B   4      37.080  43.455  -3.421  1.00  0.00           C  ',
        'ATOM      3  CA  GLU B   6      33.861  45.127  -2.233  1.00  0.00           C  ',
        'ATOM      3  CA  ALA B   7      35.081  45.036   1.305  1.00  0.00           C  ',
        ],
    [
        'ATOM      3  CA  ARG B   4      37.080  43.455  -3.421  1.00  0.00           C  ',
        'ATOM      3  CA  GLU B   6      33.861  45.127  -2.233  1.00  0.00           C  ',
        'ATOM      3  CA  ALA B   7      35.081  45.036   1.305  1.00  0.00           C  ',
        ],
    ]

expected_rep_chains_1 = [
    [
        'ATOM      3  CA  ARG B   4      37.080  43.455  -3.421  1.00  0.00           C  ',
        'ATOM      3  CA  GLU B   6      33.861  45.127  -2.233  1.00  0.00           C  ',
        'ATOM      3  CA  ALA B   7      35.081  45.036   1.305  1.00  0.00           C  ',
        ],
    [
        'ATOM      3  CA  ARG A   4      37.080  43.455  -3.421  1.00  0.00      A    C  ',
        'ATOM      3  CA  GLU A   6      33.861  45.127  -2.233  1.00  0.00      A    C  ',
        'ATOM      3  CA  ALA A   7      35.081  45.036   1.305  1.00  0.00      A    C  ',
        ],
    ]

rep_chains_no_rep = [
    [
        'ATOM      3  CA  ARG A   4      37.080  43.455  -3.421  1.00  0.00           C  ',
        'ATOM      3  CA  GLU A   6      33.861  45.127  -2.233  1.00  0.00           C  ',
        'ATOM      3  CA  ALA A   7      35.081  45.036   1.305  1.00  0.00           C  ',
        ],
    [
        'ATOM      3  CA  ARG B   4      37.080  43.455  -3.421  1.00  0.00           C  ',
        'ATOM      3  CA  GLU B   6      33.861  45.127  -2.233  1.00  0.00           C  ',
        'ATOM      3  CA  ALA B   7      35.081  45.036   1.305  1.00  0.00           C  ',
        ],
    ]


@pytest.mark.parametrize(
    'in_, expected',
    [
        (rep_chains_1, expected_rep_chains_1),
        (rep_chains_no_rep, rep_chains_no_rep),
        ],
    )
def test_correct_equal_chain_segids(in_, expected):
    result = pp.correct_equal_chain_segids(in_)

    # made for loop to facilitate visualization of errors
    for r, e in zip(result, expected):
        assert r == e


multiple_chainIDs_1 = [
    'ATOM      3  CA  ARG A   4      37.080  43.455  -3.421  1.00  0.00           C  ',
    'ATOM      3  CA  GLU B   6      33.861  45.127  -2.233  1.00  0.00           C  ',
    'TER',
    'ATOM      3  CA  ALA C   7      35.081  45.036   1.305  1.00  0.00           C  ',
    ]

multiple_chainIDs_2 = [
    'ATOM      3  CA  ARG C   4      37.080  43.455  -3.421  1.00  0.00           C  ',
    'ATOM      3  CA  GLU B   6      33.861  45.127  -2.233  1.00  0.00           C  ',
    'ATOM      3  CA  ALA A   7      35.081  45.036   1.305  1.00  0.00           C  ',
    ]

# mind the added segID
expected_multiple_chainIDs_1 = [
    'ATOM      3  CA  ARG A   4      37.080  43.455  -3.421  1.00  0.00      A    C  ',
    'ATOM      3  CA  GLU A   6      33.861  45.127  -2.233  1.00  0.00      A    C  ',
    'TER                  A                                                           ',  # why this extra space?
    'ATOM      3  CA  ALA A   7      35.081  45.036   1.305  1.00  0.00      A    C  ',
    ]

expected_multiple_chainIDs_2 = [
    'ATOM      3  CA  ARG C   4      37.080  43.455  -3.421  1.00  0.00      C    C  ',
    'ATOM      3  CA  GLU C   6      33.861  45.127  -2.233  1.00  0.00      C    C  ',
    'ATOM      3  CA  ALA C   7      35.081  45.036   1.305  1.00  0.00      C    C  ',
    ]


@pytest.mark.parametrize(
    "in_,expected",
    [
        (multiple_chainIDs_1, expected_multiple_chainIDs_1),
        (multiple_chainIDs_2, expected_multiple_chainIDs_2),
        # returns the same thing
        (expected_multiple_chainIDs_2, expected_multiple_chainIDs_2),
        ]
    )
def test_homogenize_chainIDs(in_, expected):
    """Test homogenize chain IDs."""
    result = pp.homogenize_chains(in_)

    # assert line by line to facilitate finding errors
    for r_line, e_line in zip(result, expected):
        assert len(r_line) == len(e_line)
    for r_line, e_line in zip(result, expected):
        assert r_line == e_line


nochain_chainIDs_1 = [
    'ATOM      3  CA  ARG     4      37.080  43.455  -3.421  1.00  0.00              ',
    'ATOM      3  CA  GLU     6      33.861  45.127  -2.233  1.00  0.00              ',
    'ATOM      3  CA  ALA     7      35.081  45.036   1.305  1.00  0.00              ',
    ]

nochain_chainIDs_2 = [
    'ATOM      3  CA  ARG A   4      37.080  43.455  -3.421  1.00  0.00              ',
    'ATOM      3  CA  GLU A   6      33.861  45.127  -2.233  1.00  0.00              ',
    'ATOM      3  CA  ALA A   7      35.081  45.036   1.305  1.00  0.00              ',
    ]

nochain_chainIDs_3 = [
    'ATOM      3  CA  ARG     4      37.080  43.455  -3.421  1.00  0.00      A       ',
    'ATOM      3  CA  GLU     6      33.861  45.127  -2.233  1.00  0.00      A       ',
    'ATOM      3  CA  ALA     7      35.081  45.036   1.305  1.00  0.00      A       ',
    ]

nochain_chainIDs_4 = [
    'ATOM      3  CA  ARG B   4      37.080  43.455  -3.421  1.00  0.00      A       ',
    'ATOM      3  CA  GLU B   6      33.861  45.127  -2.233  1.00  0.00      A       ',
    'ATOM      3  CA  ALA B   7      35.081  45.036   1.305  1.00  0.00      A       ',
    ]

expected_nochain_chainIDs_1 = [
    'ATOM      3  CA  ARG A   4      37.080  43.455  -3.421  1.00  0.00      A       ',
    'ATOM      3  CA  GLU A   6      33.861  45.127  -2.233  1.00  0.00      A       ',
    'ATOM      3  CA  ALA A   7      35.081  45.036   1.305  1.00  0.00      A       ',
    ]

expected_nochain_chainIDs_4 = [
    'ATOM      3  CA  ARG B   4      37.080  43.455  -3.421  1.00  0.00      B       ',
    'ATOM      3  CA  GLU B   6      33.861  45.127  -2.233  1.00  0.00      B       ',
    'ATOM      3  CA  ALA B   7      35.081  45.036   1.305  1.00  0.00      B       ',
    ]

chain_and_seg_IDs = [
    'ATOM      3  CA  ARG A   4      37.080  43.455  -3.421  1.00  0.00      A       ',
    'ATOM      3  CA  GLU A   6      33.861  45.127  -2.233  1.00  0.00      A       ',
    'ATOM      3  CA  ALA A   7      35.081  45.036   1.305  1.00  0.00      A       ',
    ]


@pytest.mark.parametrize(
    'in_,expected',
    [
        (nochain_chainIDs_1, expected_nochain_chainIDs_1),
        (nochain_chainIDs_2, expected_nochain_chainIDs_1),
        (nochain_chainIDs_3, expected_nochain_chainIDs_1),
        (nochain_chainIDs_4, expected_nochain_chainIDs_4),
        (chain_and_seg_IDs, chain_and_seg_IDs),
        ]
    )
def test_solve_nochainID(in_, expected):
    """Test solve nochainID function."""
    result = pp.solve_no_chainID_no_segID(in_)
    for r, e in zip(result, expected):
        assert r == e


# test cases to test adding the ion charges
# even indexes are the input, odd indexes are the expected results
ion_cases = [
    # 0
    'HETATM 3833 ZN+2 ZN2 A  42      21.391  -8.794  33.944  1.00 24.37          ZN  ',
    'HETATM 3833 ZN+2 ZN2 A  42      21.391  -8.794  33.944  1.00 24.37          ZN+2',

    # 1
    'HETATM 3833 ZN   ZN2 A  42      21.391  -8.794  33.944  1.00 24.37          ZN  ',
    'HETATM 3833 ZN+2 ZN2 A  42      21.391  -8.794  33.944  1.00 24.37          ZN+2',

    # 2
    'HETATM 3833 ZN+2 ZN  A  42      21.391  -8.794  33.944  1.00 24.37          ZN  ',
    'HETATM 3833 ZN+2 ZN2 A  42      21.391  -8.794  33.944  1.00 24.37          ZN+2',

    # 3
    'HETATM 3833 ZN   ZN  A  42      21.391  -8.794  33.944  1.00 24.37          ZN  ',
    'HETATM 3833 ZN    ZN A  42      21.391  -8.794  33.944  1.00 24.37          ZN  ',

    # 4
    'HETATM 3833 ZN    ZN A  42      21.391  -8.794  33.944  1.00 24.37          ZN+2',
    'HETATM 3833 ZN+2 ZN2 A  42      21.391  -8.794  33.944  1.00 24.37          ZN+2',

    # 5
    'HETATM 3834  K    K1 A  42      21.391  -8.794  33.944  1.00 24.37              ',
    'HETATM 3834  K+1  K1 A  42      21.391  -8.794  33.944  1.00 24.37           K+1',

    # 6
    'HETATM 3835 NI    NI A  42      21.391  -8.794  33.944  1.00 24.37              ',
    'HETATM 3835 NI    NI A  42      21.391  -8.794  33.944  1.00 24.37          NI  ',

    # 7
    'HETATM 3834  F-1   F A  42      21.391  -8.794  33.944  1.00 24.37              ',
    'HETATM 3834  F-1  F1 A  42      21.391  -8.794  33.944  1.00 24.37           F-1',

    # 8
    'HETATM    3 CA    CA B   4      37.080  43.455  -3.421  1.00  0.00      B   CA+2',
    'HETATM    3 CA+2 CA2 B   4      37.080  43.455  -3.421  1.00  0.00      B   CA+2',

    # 9 lines not concerning ions - note element is not defined
    # carbon alpha conflicting with calcium
    'ATOM      3  CA  ARG B   4      37.080  43.455  -3.421  1.00  0.00      B       ',
    'ATOM      3  CA  ARG B   4      37.080  43.455  -3.421  1.00  0.00      B       ',

    # 10 other lines not concerning ions - atom not CA (carbon alpha)
    'ATOM      1  N   ALA A   1       0.000   0.000   0.000  0.00  0.00           N  ',
    'ATOM      1  N   ALA A   1       0.000   0.000   0.000  0.00  0.00           N  ',

    # 11 lines not concerning ions - here the element is defined
    'ATOM      3  CA  ARG B   4      37.080  43.455  -3.421  1.00  0.00      B    C  ',
    'ATOM      3  CA  ARG B   4      37.080  43.455  -3.421  1.00  0.00      B    C  ',

    # 12 exception
    'HETATM 3834  F    FA A  42      21.391  -8.794  33.944  1.00 24.37              ',
    'HETATM 3834  F     F A  42      21.391  -8.794  33.944  1.00 24.37           F  ',
    ]


# made this way to make each line a separate tests
@pytest.fixture(params=list(zip(ion_cases[::2], ion_cases[1::2])))
def ion_cases_fixture(request):
    """Give ion test cases."""
    return request.param[0], request.param[1]


def test_correct_ion_charges(ion_cases_fixture):
    input_, expected = ion_cases_fixture

    # 1. puts input_ inside a list because the function only accepts lists
    # 2. the function is wrapped around list() because it is a generator
    result = list(pp.add_charges_to_ions([input_]))

    assert len(result) == 1
    # compares two strings
    assert result[0] == expected


def test_process_ion_case_atom():
    """Test processing an exception case."""
    inp = 'HETATM 3834  F    FA A  42      21.391  -8.794  33.944  1.00 24.37              '  # noqa: E501
    exp = 'HETATM 3834  F     F A  42      21.391  -8.794  33.944  1.00 24.37           F  '  # noqa: E501
    result = pp._process_ion_case_atom(inp)
    assert result == exp


@pytest.mark.parametrize(
    "lines, expected",
    [
        (
            'ATOM      3  CA  GLU A   6      33.861  45.127  -2.233  1.00  0.00      A       ',
            'ATOM      3  CA  GLU B   6      33.861  45.127  -2.233  1.00  0.00      A       ',
            ),
        ]
    )
def test_wrep_chain(lines, expected):
    """Test if dry report works."""
    result = pp.wrep_pdb_chain([lines], "B", report=True)
    assert result[0] == expected


def test_process_pdbs():
    """."""
    result = pp.process_pdbs(broken_pdb)
    assert len(result) == 1

    expected = corrected_pdb.read_text().rstrip(os.linesep).split(os.linesep)
    Path('testpreprocessing.pdb').write_text(os.linesep.join(result[0]))

    for i, (rline, eline) in enumerate(zip_longest(result[0], expected)):
        assert rline == eline, i
