"""Test preprocessing operations."""
import logging
import os
from pathlib import Path

import pytest

from haddock.gear import preprocessing as pp

from . import broken_pdb, data_folder, good_pdb


def test_open_or_give_3():
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
        assert len(r) == 241
    in1.close()


def test_wrep_pdb_tidy():
    pp.wrep_pdb_tidy(pp._open_or_give([broken_pdb])[0], strict=False)


def test_wrep_pdb_tidy_strict():
    pp.wrep_pdb_tidy(pp._open_or_give([broken_pdb])[0], strict=True)


#@pytest.mark.skip
def test_all(caplog):
    """."""
    caplog.set_level(logging.WARNING)

    result = pp.process_pdbs(broken_pdb)

    output = Path(
        broken_pdb.parent,
        broken_pdb.stem + '_processed').with_suffix(broken_pdb.suffix)

    output.write_text(os.linesep.join(result[0]))

    #lines_results = output.read_text().split('\n')
    #lines_expected = good_pdb.read_text().split('\n')

    # asserts the number of WARNING messages is appropriate
    # see: https://docs.pytest.org/en/6.2.x/logging.html
    #assert len(caplog.text.split(os.linesep)) > 0

    #assert lines_results[-1] == lines_expected[-1]
    #assert len(lines_results) == len(lines_expected)
    #for i in range(len(lines_results)):
    #    assert lines_results[i] == lines_expected[i]



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
    assert str(error.value) == "MODEL 2 differs from MODEL 1."


def test_check_models_wrong_3():
    _lines = models_wrong_3.split(os.linesep)
    with pytest.raises(pp.ModelsDifferError) as error:
        pp.models_should_have_the_same_labels(_lines)
    assert str(error.value) == "MODEL 3 differs from MODEL 1."


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
    'TER                  A                                                           ',  # wtf this extra space?
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
        ]
    )
def test_homogenize_chainIDs(in_, expected):
    result = pp.homogenize_chains(in_)
    print(result)
    for r, e in zip(result, expected):
        assert len(result) == len(expected)
    for r, e in zip(result, expected):
        assert result == expected


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


@pytest.mark.parametrize(
    'in_,expected',
    [
        (nochain_chainIDs_1, expected_nochain_chainIDs_1),
        (nochain_chainIDs_2, expected_nochain_chainIDs_1),
        (nochain_chainIDs_3, expected_nochain_chainIDs_1),
        (nochain_chainIDs_4, expected_nochain_chainIDs_4),
        ]
    )
def test_solve_nochainID(in_, expected):
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
    ]

# made this way to make each line a separate tests
@pytest.fixture(params=list(zip(ion_cases[::2], ion_cases[1::2])))
def ion_cases_fixture(request):
    return request.param[0], request.param[1]


def tests_correct_ion_charges(ion_cases_fixture):
    input_, expected = ion_cases_fixture

    # 1. puts input_ inside a list because the function only accepts lists
    # 2. the function is wrapped around list() because it is a generator
    result = list(pp.add_charges_to_ions([input_]))

    assert len(result) == 1
    # compares two strings
    assert result[0] == expected
