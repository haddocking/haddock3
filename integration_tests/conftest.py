from pathlib import Path

import pytest

from haddock.modules.analysis.caprieval.capri import load_contacts
from haddock.libs.libontology import PDBFile


def calc_fnat_with_caprieval(model: Path, native: Path) -> float:
    model_pdb = PDBFile(model)
    native_pdb = PDBFile(native)

    model_contacts = load_contacts(model_pdb)
    native_contacts = load_contacts(native_pdb)

    intersection = native_contacts & model_contacts

    fnat = len(intersection) / float(len(model_contacts))

    return fnat


@pytest.fixture
def calc_fnat():
    return calc_fnat_with_caprieval
