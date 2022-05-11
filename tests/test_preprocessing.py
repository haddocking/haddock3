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
        assert len(r) == 240
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
