"""Test supported molecules."""
import os

import pytest

from haddock.core import supported_molecules as SM


ZN = """RESIdue ZN2 {zinc 2+}
  GROUp
    ATOM ZN+2 TYPE=ZN+2 CHARge=+0.96 END ! comment
END {ZN2}"""

CU = """RESIdue CU2 {copper 2+}
  GROUp
    ATOM CU+2 TYPE=CU+2 CHARge=+2.0 END
END {CU2}"""


@pytest.mark.parametrize(
    "text, expected",
    [
        (ZN, ("ZN2", ("ZN+2",), 0.96)),
        (CU, ("CU2", ("CU+2",), 2)),
        ]
    )
def test_read_top_file(text, expected):
    """Test read top file to CNSTopologyResidue."""
    result = SM._read_residues_from_top_file(text.split(os.linesep))
    assert isinstance(result, list)
    assert isinstance(result[0], SM.CNSTopologyResidue)

    resname, atom, charge = expected
    assert result[0].resname == resname
    assert result[0].atoms == atom
    assert result[0].charge == charge


def test_read_top_file_2():
    text = ZN + os.linesep + CU
    result = SM._read_residues_from_top_file(text.split(os.linesep))

    assert isinstance(result, list)
    assert len(result) == 2
    assert all(isinstance(r, SM.CNSTopologyResidue) for r in result)

    assert result[0].resname == "ZN2"
    assert result[0].atoms == ("ZN+2",)
    assert result[0].charge == 0.96

    assert result[1].resname == "CU2"
    assert result[1].atoms == ("CU+2",)
    assert result[1].charge == 2
