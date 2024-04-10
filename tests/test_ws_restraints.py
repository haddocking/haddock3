from base64 import b64encode
import gzip
from pathlib import Path
from textwrap import dedent

from fastapi.testclient import TestClient
import pytest

from haddock.clis.restraints.webservice import app

from . import golden_data

@pytest.fixture
def example_pdb_file_gzipped_base64():
    file = Path(golden_data, "protprot_onechain.pdb")
    content = file.read_bytes()
    gzipped = gzip.compress(content)
    return b64encode(gzipped).decode()

@pytest.fixture
def client():
    return TestClient(app)

def test_calc_accessibility(client: TestClient, example_pdb_file_gzipped_base64: str):
    body = {
        "structure": example_pdb_file_gzipped_base64,
        "cutoff": 0.1,
    }
    response = client.post("/calc_accessibility", json=body)
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    result = response.json()
    assert all(isinstance(i, int) for i in result["A"])
    assert len(result["A"]) == 118

def test_actpass_to_ambig(client: TestClient):
    body = {
        "active1": [
            1,
        ],
        "active2": [2],
        "passive1": [3],
        "passive2": [4],
        "segid1": "A",
        "segid2": "B",
    }
    response = client.post("/actpass_to_ambig", json=body)
    assert response.status_code == 200
    assert response.headers["content-type"].startswith("text/plain")
    expected = dedent(
        """\
        assign (resi 1 and segid A)
        (
               (resi 2 and segid B)
                or
               (resi 4 and segid B)
        ) 2.0 2.0 0.0

        assign (resi 2 and segid B)
        (
               (resi 1 and segid A)
                or
               (resi 3 and segid A)
        ) 2.0 2.0 0.0"""
    )
    assert response.text == expected
