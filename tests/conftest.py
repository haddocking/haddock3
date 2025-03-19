import shutil
import tempfile
from pathlib import Path
from typing import Any, Generator

import pytest

from haddock.libs.libontology import PDBFile

from . import golden_data


@pytest.fixture(name="protprot_input_list")
def fixture_protprot_input_list():
    """Prot-prot input."""

    with (
        tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_1,
        tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_2,
    ):
        dst_prot_1 = tmp_1.name
        dst_prot_2 = tmp_2.name

        src_prot_1 = Path(golden_data, "protprot_complex_1.pdb")
        src_prot_2 = Path(golden_data, "protprot_complex_2.pdb")

        shutil.copy(src_prot_1, dst_prot_1)
        shutil.copy(src_prot_2, dst_prot_2)

        pdb_obj_1 = PDBFile(file_name=dst_prot_1, path=Path(dst_prot_1).parent)
        pdb_obj_2 = PDBFile(file_name=dst_prot_2, path=Path(dst_prot_2).parent)

        yield [pdb_obj_1, pdb_obj_2]


@pytest.fixture(name="protprot_1bkd_input_list")
def fixture_protprot_1bkd_input_list():
    """
    Prot-prot input for target 1bkd.

    Heterogeneous ensemble and big protein.
    """
    with (
        tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_1,
        tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_2,
    ):
        dst_prot_1 = tmp_1.name
        dst_prot_2 = tmp_2.name

        src_prot_1 = Path(golden_data, "protprot_1bkd_1.pdb")
        src_prot_2 = Path(golden_data, "protprot_1bkd_2.pdb")

        shutil.copy(src_prot_1, dst_prot_1)
        shutil.copy(src_prot_2, dst_prot_2)

        pdb_obj_1 = PDBFile(file_name=dst_prot_1, path=Path(dst_prot_1).parent)
        pdb_obj_2 = PDBFile(file_name=dst_prot_2, path=Path(dst_prot_2).parent)

        yield [pdb_obj_1, pdb_obj_2]


@pytest.fixture(name="protdna_input_list")
def fixture_protdna_input_list():
    """Prot-DNA input."""
    with (
        tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_1,
        tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_2,
    ):
        dst_prot_1 = tmp_1.name
        dst_prot_2 = tmp_2.name

        src_prot_1 = Path(golden_data, "protdna_complex_1.pdb")
        src_prot_2 = Path(golden_data, "protdna_complex_2.pdb")

        shutil.copy(src_prot_1, dst_prot_1)
        shutil.copy(src_prot_2, dst_prot_2)

        pdb_obj_1 = PDBFile(file_name=dst_prot_1, path=Path(dst_prot_1).parent)
        pdb_obj_2 = PDBFile(file_name=dst_prot_2, path=Path(dst_prot_2).parent)

        yield [pdb_obj_1, pdb_obj_2]


@pytest.fixture(name="protlig_input_list")
def fixture_protlig_input_list():
    """Protein-Ligand input."""
    with (
        tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_1,
        tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_2,
    ):
        dst_prot_1 = tmp_1.name
        dst_prot_2 = tmp_2.name

        src_prot_1 = Path(golden_data, "protlig_complex_1.pdb")
        src_prot_2 = Path(golden_data, "protlig_complex_2.pdb")

        shutil.copy(src_prot_1, dst_prot_1)
        shutil.copy(src_prot_2, dst_prot_2)

        pdb_obj_1 = PDBFile(file_name=dst_prot_1, path=Path(dst_prot_1).parent)
        pdb_obj_2 = PDBFile(file_name=dst_prot_2, path=Path(dst_prot_2).parent)

        yield [pdb_obj_1, pdb_obj_2]


@pytest.fixture(name="protprot_onechain_list")
def fixture_protprot_onechain_list() -> Generator[list[PDBFile], Any, Any]:
    """Protein-Protein complex with a single chain ID."""
    with (
        tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_1,
        tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_2,
    ):
        dst_prot_1 = tmp_1.name
        dst_prot_2 = tmp_2.name

        src_prot_1 = Path(golden_data, "protprot_complex_1.pdb")
        src_prot_2 = Path(golden_data, "protprot_onechain.pdb")

        shutil.copy(src_prot_1, dst_prot_1)
        shutil.copy(src_prot_2, dst_prot_2)

        pdb_obj_1 = PDBFile(file_name=dst_prot_1, path=Path(dst_prot_1).parent)
        pdb_obj_2 = PDBFile(file_name=dst_prot_2, path=Path(dst_prot_2).parent)

        yield [pdb_obj_1, pdb_obj_2]
