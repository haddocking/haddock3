import os
import shutil
import tempfile
from pathlib import Path
from typing import Any, Generator

import pytest

from haddock.libs.libcapri import CAPRI
from haddock.libs.libontology import PDBFile

from . import golden_data


def remove_aln_files(class_name):
    """Helper function to remove intermediary alignment files."""
    file_l = [
        Path(class_name.path, "blosum62.izone"),
        Path(class_name.path, "blosum62_A.aln"),
        Path(class_name.path, "blosum62_B.aln"),
    ]
    for f in file_l:
        if f.exists():
            os.unlink(f)


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

@pytest.fixture(name="protprot_input_list_cg")
def fixture_protprot_input_list_cg():
    """Prot-prot input."""

    with (
        tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_1,
        tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_2,
    ):
        dst_prot_1 = tmp_1.name
        dst_prot_2 = tmp_2.name

        src_prot_1 = Path(golden_data, "protprot_complex_1_cg.pdb")
        src_prot_2 = Path(golden_data, "protprot_complex_2_cg.pdb")

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


@pytest.fixture(name="protprot_1bkd_input_list_cg")
def fixture_protprot_1bkd_input_list_cg():
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

        src_prot_1 = Path(golden_data, "protprot_1bkd_1_cg.pdb")
        src_prot_2 = Path(golden_data, "protprot_1bkd_2_cg.pdb")

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


@pytest.fixture(name="protdna_input_list_cg")
def fixture_protdna_input_list_cg():
    """Prot-DNA input."""
    with (
        tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_1,
        tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_2,
    ):
        dst_prot_1 = tmp_1.name
        dst_prot_2 = tmp_2.name

        src_prot_1 = Path(golden_data, "protdna_complex_1_cg.pdb")
        src_prot_2 = Path(golden_data, "protdna_complex_2_cg.pdb")

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


@pytest.fixture(name="prot_HETATMlig_input_list")
def fixture_prot_HETATMlig_input_list():
    """Protein-Ligand with HETATM input."""
    with (
        tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_1,
        tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp_2,
    ):
        dst_prot_1 = tmp_1.name
        dst_prot_2 = tmp_2.name

        src_prot_1 = Path(golden_data, "protlig_complex_1_HETATM.pdb")
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


# ---------------------------------------------------------------------------
# CAPRI object params fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(name="params")
def fixture_params():
    """???"""
    return {
        "receptor_chain": "A",
        "ligand_chains": ["B"],
        "aln_method": "sequence",
        "allatoms": False,
        "keep_hetatm": False,
    }


@pytest.fixture(name="params_all")
def fixture_params_all():
    """???"""
    return {
        "receptor_chain": "A",
        "ligand_chains": ["B"],
        "aln_method": "sequence",
        "allatoms": True,
        "keep_hetatm": False,
    }


@pytest.fixture(name="params_hetatm")
def fixture_params_hetatm():
    """Parameters when HETATM must be kept."""
    return {
        "receptor_chain": "A",
        "ligand_chains": ["B"],
        "aln_method": "sequence",
        "allatoms": False,
        "keep_hetatm": True,
    }


# ---------------------------------------------------------------------------
# CAPRI object fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(name="protdna_caprimodule")
def fixture_protdna_caprimodule(protdna_input_list, params):
    """Protein-DNA CAPRI module."""
    reference = protdna_input_list[0]
    model = protdna_input_list[1]
    _path = protdna_input_list[0].path
    capri = CAPRI(
        identificator=42,
        reference=reference,
        model=model,
        path=_path,
        params=params,
    )

    yield capri

    remove_aln_files(capri)


@pytest.fixture(name="protlig_caprimodule")
def fixture_protlig_caprimodule(protlig_input_list, params):
    """Protein-Ligand CAPRI module."""
    reference = protlig_input_list[0]
    model = protlig_input_list[1]
    _path = protlig_input_list[0].path
    capri = CAPRI(
        identificator=42,
        reference=reference,
        model=model,
        path=_path,
        params=params,
    )

    yield capri

    remove_aln_files(capri)


@pytest.fixture(name="protprot_caprimodule")
def fixture_protprot_caprimodule(protprot_input_list, params):
    """Protein-Protein CAPRI module."""
    reference = protprot_input_list[0]
    model = protprot_input_list[1]
    _path = protprot_input_list[0].path
    capri = CAPRI(
        identificator=str(42),
        reference=reference,
        model=model,
        path=_path,
        params=params,
    )

    yield capri

    remove_aln_files(capri)


@pytest.fixture(name="protprot_allatm_caprimodule")
def fixture_protprot_allatm_caprimodule(protprot_input_list, params_all):
    """Protein-Protein CAPRI module."""
    reference = protprot_input_list[0]
    model = protprot_input_list[1]
    _path = protprot_input_list[0].path
    capri = CAPRI(
        identificator=42,
        reference=reference,
        model=model,
        path=_path,
        params=params_all,
    )

    yield capri

    remove_aln_files(capri)


@pytest.fixture(name="protprot_1bkd_caprimodule")
def fixture_protprot_1bkd_caprimodule(protprot_1bkd_input_list, params):
    """Protein-Protein CAPRI module for target 1BKD."""
    reference = protprot_1bkd_input_list[0]
    model = protprot_1bkd_input_list[1]
    _path = protprot_1bkd_input_list[0].path
    capri = CAPRI(
        identificator=str(42),
        reference=reference,
        model=model,
        path=_path,
        params=params,
    )

    yield capri

    remove_aln_files(capri)


@pytest.fixture(name="protprot_onechain_ref_caprimodule")
def fixture_protprot_onechain_ref_caprimodule(protprot_onechain_list, params):
    """Protein-Protein CAPRI module with a single chain structure as ref."""

    # reference, model = protprot_onechain_list
    model, reference = protprot_onechain_list
    _path = reference.path

    capri = CAPRI(
        identificator=str(42),
        reference=reference,
        model=model,
        path=_path,
        params=params,
    )

    yield capri

    remove_aln_files(capri)


@pytest.fixture(name="protprot_onechain_mod_caprimodule")
def fixture_protprot_onechain_mod_caprimodule(protprot_onechain_list, params):
    """Protein-Protein CAPRI module with a single chain structure as model."""

    reference, model = protprot_onechain_list
    _path = reference.path

    capri = CAPRI(
        identificator=str(42),
        reference=reference,
        model=model,
        path=_path,
        params=params,
    )

    yield capri

    remove_aln_files(capri)


@pytest.fixture(name="protdna_caprimodule_cg")
def fixture_protdna_caprimodule_cg(protdna_input_list_cg, params):
    """Protein-DNA CAPRI module (coarse-grained)."""
    reference = protdna_input_list_cg[0]
    model = protdna_input_list_cg[1]
    _path = protdna_input_list_cg[0].path
    capri = CAPRI(
        identificator=42,
        reference=reference,
        model=model,
        path=_path,
        params=params,
        ff="martini2",
    )

    yield capri

    remove_aln_files(capri)


@pytest.fixture(name="protprot_caprimodule_cg")
def fixture_protprot_caprimodule_cg(protprot_input_list_cg, params):
    """Protein-Protein CAPRI module (coarse-grained)."""
    reference = protprot_input_list_cg[0]
    model = protprot_input_list_cg[1]
    _path = protprot_input_list_cg[0].path
    capri = CAPRI(
        identificator=str(42),
        reference=reference,
        model=model,
        path=_path,
        params=params,
        ff="martini2",
    )

    yield capri

    remove_aln_files(capri)


@pytest.fixture(name="protprot_allatm_caprimodule_cg")
def fixture_protprot_allatm_caprimodule_cg(protprot_input_list_cg, params_all):
    """Protein-Protein CAPRI module, all atoms (coarse-grained)."""
    reference = protprot_input_list_cg[0]
    model = protprot_input_list_cg[1]
    _path = protprot_input_list_cg[0].path
    capri = CAPRI(
        identificator=42,
        reference=reference,
        model=model,
        path=_path,
        params=params_all,
        ff="martini2",
    )

    yield capri

    remove_aln_files(capri)


@pytest.fixture(name="prot_hetatmlig_caprimodule")
def fixture_prot_hetatmlig_caprimodule(prot_HETATMlig_input_list, params_hetatm):
    """Protein-Ligand CAPRI module (HETATM ligand)."""
    reference = prot_HETATMlig_input_list[0]
    model = prot_HETATMlig_input_list[1]
    _path = prot_HETATMlig_input_list[1].path
    capri = CAPRI(
        identificator=42,
        reference=reference,
        model=model,
        path=_path,
        params=params_hetatm,
    )

    yield capri

    remove_aln_files(capri)
