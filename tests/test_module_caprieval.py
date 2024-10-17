"""Test the CAPRI module."""

import os
import random
import shutil
import tempfile
import uuid
from pathlib import Path

import numpy as np
import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.caprieval.capri import (
    CAPRI,
    calc_stats,
    capri_cluster_analysis,
    extract_data_from_capri_class,
    get_previous_cns_step,
    load_contacts,
    rank_according_to_score,
)

from . import golden_data


def remove_aln_files(class_name):
    """Remove intermediary alignment files."""
    file_l = [
        Path(class_name.path, "blosum62.izone"),
        Path(class_name.path, "blosum62_A.aln"),
        Path(class_name.path, "blosum62_B.aln"),
    ]
    for f in file_l:
        if f.exists():
            os.unlink(f)


def read_capri_file(fname):
    """???"""
    file_content = [
        e.split()[1:] for e in open(fname).readlines() if not e.startswith("#")
    ]
    return file_content


@pytest.fixture(name="params")
def fixture_params():
    """???"""
    return {
        "receptor_chain": "A",
        "ligand_chains": ["B"],
        "aln_method": "sequence",
        "allatoms": False,
    }


@pytest.fixture(name="params_all")
def fixture_params_all():
    """???"""
    return {
        "receptor_chain": "A",
        "ligand_chains": ["B"],
        "aln_method": "sequence",
        "allatoms": True,
    }


@pytest.fixture(name="protdna_caprimodule")
def fixture_protdna_caprimodule(protdna_input_list, params):
    """Protein-DNA CAPRI module."""
    reference = protdna_input_list[0]
    model = protdna_input_list[1]
    _path = protdna_input_list[0].path
    capri = CAPRI(
        identificator=str(42),
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
        identificator=str(42),
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
        identificator=str(42),
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


# Q: Not used anywhere?
# @pytest.fixture
# def protprot_caprimodule_parallel(protprot_input_list):
#     """Protein-Protein CAPRI module."""
#     reference = protprot_input_list[0].rel_path
#     model = protprot_input_list[1].rel_path
#     capri = CAPRI(
#         reference=reference,
#         model=model,
#         receptor_chain="A",
#         ligand_chains=["B"],
#         aln_method="sequence",
#         path=golden_data,
#         identificator=0,
#         core_model_idx=0,
#     )
#     yield capri
#     remove_aln_files(capri)


def test_protprot_irmsd(protprot_caprimodule):
    """Test protein-protein i-rmsd calculation."""
    # using standard cutoff
    protprot_caprimodule.calc_irmsd(cutoff=10.0)
    assert np.isclose(protprot_caprimodule.irmsd, 8.33, atol=0.01)
    # using default cutoff = 5.0
    protprot_caprimodule.calc_irmsd()
    assert np.isclose(protprot_caprimodule.irmsd, 7.38, atol=0.01)


def test_protprot_lrmsd(protprot_caprimodule):
    """Test protein-protein l-rmsd calculation."""
    protprot_caprimodule.calc_lrmsd()
    assert np.isclose(protprot_caprimodule.lrmsd, 20.94, atol=0.01)


def test_protprot_ilrmsd(protprot_caprimodule):
    """Test protein-protein i-l-rmsd calculation."""
    protprot_caprimodule.calc_ilrmsd()
    assert np.isclose(protprot_caprimodule.ilrmsd, 18.25, atol=0.01)


def test_protprot_fnat(protprot_caprimodule):
    """Test protein-protein fnat calculation."""
    protprot_caprimodule.calc_fnat()
    assert np.isclose(protprot_caprimodule.fnat, 0.05, atol=0.01)


def test_protprot_dockq(protprot_caprimodule):
    """Test protein-protein dockq calculation."""
    protprot_caprimodule.irmsd = 7.38
    protprot_caprimodule.fnat = 0.05
    protprot_caprimodule.lrmsd = 15.9
    protprot_caprimodule.calc_dockq()
    assert np.isclose(protprot_caprimodule.dockq, 0.10, atol=0.01)


def test_protprot_global_rmsd(protprot_caprimodule):
    """Test protein-protein l-rmsd calculation."""
    protprot_caprimodule.calc_global_rmsd()
    assert np.isclose(protprot_caprimodule.rmsd, 8.62, atol=0.01)


def test_protprot_bb_vs_all_atoms(
    protprot_caprimodule,
    protprot_allatm_caprimodule,
):
    """Test difference between all and backbone atoms selection."""
    assert protprot_caprimodule.allatoms is False
    assert protprot_allatm_caprimodule.allatoms is True
    assert protprot_caprimodule.atoms != protprot_allatm_caprimodule.atoms
    only_bb_atoms = [a for atn in protprot_caprimodule.atoms.values() for a in atn]
    all_atoms = [a for atn in protprot_allatm_caprimodule.atoms.values() for a in atn]
    # Check length
    assert len(only_bb_atoms) < len(all_atoms)
    # Make sure all bb are also included in all atoms
    for rname, atoms in protprot_caprimodule.atoms.items():
        for atn in atoms:
            assert atn in protprot_allatm_caprimodule.atoms[rname]


def test_protprot_allatoms(protprot_allatm_caprimodule):
    """Test protein-protein all atoms calculations."""
    # Compute values with all atoms calculations
    protprot_allatm_caprimodule.calc_fnat()
    protprot_allatm_caprimodule.calc_lrmsd()
    protprot_allatm_caprimodule.calc_irmsd()
    protprot_allatm_caprimodule.calc_ilrmsd()
    # process checks
    assert np.isclose(protprot_allatm_caprimodule.lrmsd, 21.10, atol=0.01)


def test_protprot_1bkd_irmsd(protprot_1bkd_caprimodule):
    """Test protein-protein i-rmsd calculation."""
    protprot_1bkd_caprimodule.calc_irmsd(cutoff=10.0)
    assert np.isclose(protprot_1bkd_caprimodule.irmsd, 8.16, atol=0.01)


def test_protprot_1bkd_lrmsd(protprot_1bkd_caprimodule):
    """Test protein-protein l-rmsd calculation."""
    protprot_1bkd_caprimodule.calc_lrmsd()
    assert np.isclose(protprot_1bkd_caprimodule.lrmsd, 28.47, atol=0.01)


def test_protprot_1bkd_ilrmsd(protprot_1bkd_caprimodule):
    """Test protein-protein i-l-rmsd calculation."""
    protprot_1bkd_caprimodule.calc_ilrmsd()
    assert np.isclose(protprot_1bkd_caprimodule.ilrmsd, 21.71, atol=0.01)


def test_protprot_1bkd_fnat(protprot_1bkd_caprimodule):
    """Test protein-protein fnat calculation."""
    protprot_1bkd_caprimodule.calc_fnat()
    assert np.isclose(protprot_1bkd_caprimodule.fnat, 0.07, atol=0.01)


def test_protprot_1bkd_dockq(protprot_1bkd_caprimodule):
    """Test protein-protein dockq calculation."""
    protprot_1bkd_caprimodule.irmsd = 8.16
    protprot_1bkd_caprimodule.fnat = 0.07
    protprot_1bkd_caprimodule.lrmsd = 28.47
    protprot_1bkd_caprimodule.calc_dockq()
    assert np.isclose(protprot_1bkd_caprimodule.dockq, 0.06, atol=0.01)


def test_protprot_1bkd_allatoms(protprot_1bkd_caprimodule):
    """Test protein-protein all atoms calculations."""
    # modify module allatoms parameter
    protprot_1bkd_caprimodule.allatoms = True
    protprot_1bkd_caprimodule.atoms = protprot_1bkd_caprimodule._load_atoms(
        protprot_1bkd_caprimodule.model,
        protprot_1bkd_caprimodule.reference,
        full=protprot_1bkd_caprimodule.allatoms,
    )
    # Compute values with all atoms calculations
    protprot_1bkd_caprimodule.calc_fnat()
    protprot_1bkd_caprimodule.calc_lrmsd()
    protprot_1bkd_caprimodule.calc_irmsd()
    protprot_1bkd_caprimodule.calc_ilrmsd()
    # Check l-rmsd with Pymol cmd
    # align protprot_1bkd_1 and chain A and not hydrogen, protprot_1bkd_2 and chain A and not hydrogen, cycles=0
    # rms_cur protprot_1bkd_1 and chain R and not hydrogen, protprot_1bkd_2 and chain R and not hydrogen
    assert np.isclose(protprot_1bkd_caprimodule.lrmsd, 28.530, atol=0.01)
    assert np.isclose(protprot_1bkd_caprimodule.irmsd, 8.181, atol=0.01)
    assert np.isclose(protprot_1bkd_caprimodule.ilrmsd, 21.5878, atol=0.01)
    assert np.isclose(protprot_1bkd_caprimodule.fnat, 0.07, atol=0.01)


def test_protlig_irmsd(protlig_caprimodule):
    """Test protein-ligand i-rmsd calculation."""
    protlig_caprimodule.calc_irmsd()
    assert np.isclose(protlig_caprimodule.irmsd, 0.22, atol=0.01)


def test_protlig_lrmsd(protlig_caprimodule):
    """Test protein-ligand l-rmsd calculation."""
    protlig_caprimodule.calc_lrmsd()
    assert np.isclose(protlig_caprimodule.lrmsd, 0.49, atol=0.01)


def test_protlig_ilrmsd(protlig_caprimodule):
    """Test protein-ligand i-l-rmsd calculation."""
    protlig_caprimodule.calc_ilrmsd()
    assert np.isclose(protlig_caprimodule.ilrmsd, 0.49, atol=0.01)


def test_protlig_fnat(protlig_caprimodule):
    """Test protein-ligand fnat calculation."""
    protlig_caprimodule.calc_fnat()
    assert np.isclose(protlig_caprimodule.fnat, 1.0, atol=0.01)


def test_protlig_allatoms(protlig_caprimodule):
    """Test protein-protein all atoms calculations."""
    # modify module allatoms parameter
    protlig_caprimodule.allatoms = True
    protlig_caprimodule.atoms = protlig_caprimodule._load_atoms(
        protlig_caprimodule.model,
        protlig_caprimodule.reference,
        full=protlig_caprimodule.allatoms,
    )
    # Compute values with all atoms calculations
    protlig_caprimodule.calc_fnat()
    protlig_caprimodule.calc_lrmsd()
    protlig_caprimodule.calc_irmsd()
    protlig_caprimodule.calc_ilrmsd()
    assert np.isclose(protlig_caprimodule.irmsd, 0.16, atol=0.01)
    assert np.isclose(protlig_caprimodule.lrmsd, 0.49, atol=0.01)
    assert np.isclose(protlig_caprimodule.ilrmsd, 0.49, atol=0.01)
    assert np.isclose(protlig_caprimodule.fnat, 1.0, atol=0.01)


def test_protdna_irmsd(protdna_caprimodule):
    """Test protein-dna i-rmsd calculation."""
    protdna_caprimodule.calc_irmsd()
    assert np.isclose(protdna_caprimodule.irmsd, 2.05, atol=0.01)


def test_protdna_lrmsd(protdna_caprimodule):
    """Test protein-dna l-rmsd calculation."""
    protdna_caprimodule.calc_lrmsd()
    assert np.isclose(protdna_caprimodule.lrmsd, 6.13, atol=0.01)


def test_protdna_ilrmsd(protdna_caprimodule):
    """Test protein-dna i-l-rmsd calculation."""
    protdna_caprimodule.calc_ilrmsd()
    assert np.isclose(protdna_caprimodule.ilrmsd, 5.97, atol=0.01)


def test_protdna_fnat(protdna_caprimodule):
    """Test protein-dna fnat calculation."""
    protdna_caprimodule.calc_fnat()
    assert np.isclose(protdna_caprimodule.fnat, 0.49, atol=0.01)


def test_protdna_allatoms(protdna_caprimodule):
    """Test protein-protein all atoms calculations."""
    # modify module allatoms parameter
    protdna_caprimodule.allatoms = True
    protdna_caprimodule.atoms = protdna_caprimodule._load_atoms(
        protdna_caprimodule.model,
        protdna_caprimodule.reference,
        full=protdna_caprimodule.allatoms,
    )
    # Compute values with all atoms calculations
    protdna_caprimodule.calc_fnat()
    protdna_caprimodule.calc_lrmsd()
    protdna_caprimodule.calc_irmsd()
    protdna_caprimodule.calc_ilrmsd()
    # Check l-rmsd with Pymol cmd
    # align protdna_complex_1 and chain A and not hydrogen, protdna_complex_2 and chain A and not hydrogen, cycles=0
    # rms_cur protdna_complex_1 and chain B and not hydrogen, protdna_complex_2 and chain B and not hydrogen
    assert np.isclose(protdna_caprimodule.lrmsd, 6.194, atol=0.01)
    assert np.isclose(protdna_caprimodule.irmsd, 2.321, atol=0.01)
    assert np.isclose(protdna_caprimodule.ilrmsd, 5.955, atol=0.01)
    assert np.isclose(protdna_caprimodule.fnat, 0.4878, atol=0.01)




def test_identify_protprotinterface(protprot_caprimodule, protprot_input_list):
    """Test the interface identification."""
    protprot_complex = protprot_input_list[0]

    observed_interface = protprot_caprimodule.identify_interface(
        protprot_complex, cutoff=5.0
    )
    expected_interface = {
        "A": [37, 38, 39, 43, 45, 71, 75, 90, 94, 96, 132],
        "B": [52, 51, 16, 54, 53, 56, 11, 12, 17, 48],
    }

    for ch in expected_interface.keys():
        assert sorted(observed_interface[ch]) == sorted(expected_interface[ch])


def test_identify_protdnainterface(protdna_caprimodule, protdna_input_list):
    """Test the interface identification."""
    protdna_complex = protdna_input_list[0]

    observed_interface = protdna_caprimodule.identify_interface(
        protdna_complex, cutoff=5.0
    )
    expected_interface = {
        "A": [10, 16, 17, 18, 27, 28, 29, 30, 32, 33, 38, 39, 41, 42, 43, 44],
        "B": [4, 3, 2, 33, 32, 5, 6, 34, 35, 31, 7, 30],
    }

    for ch in expected_interface.keys():
        assert sorted(observed_interface[ch]) == sorted(expected_interface[ch])


def test_identify_protliginterface(protlig_caprimodule, protlig_input_list):
    """Test the interface identification."""
    protlig_complex = protlig_input_list[0]

    observed_interface = protlig_caprimodule.identify_interface(
        protlig_complex, cutoff=5.0
    )
    expected_interface = {
        "A": [
            118,
            119,
            151,
            152,
            178,
            179,
            222,
            224,
            227,
            246,
            276,
            277,
            292,
            294,
            348,
            371,
            406,
        ],
        "B": [500],
    }

    for ch in expected_interface.keys():
        assert sorted(observed_interface[ch]) == sorted(expected_interface[ch])


def test_load_contacts(protprot_input_list):
    """Test loading contacts."""
    protprot_complex = protprot_input_list[0]
    observed_con_set = load_contacts(protprot_complex, cutoff=5.0)
    expected_con_set = {
        ("A", 39, "B", 52),
        ("A", 43, "B", 54),
        ("A", 45, "B", 56),
        ("A", 38, "B", 16),
        ("A", 75, "B", 17),
        ("A", 94, "B", 16),
        ("A", 39, "B", 51),
        ("A", 39, "B", 54),
        ("A", 90, "B", 17),
        ("A", 96, "B", 17),
        ("A", 45, "B", 12),
        ("A", 39, "B", 53),
        ("A", 38, "B", 51),
        ("A", 132, "B", 48),
        ("A", 71, "B", 17),
        ("A", 132, "B", 51),
        ("A", 90, "B", 16),
        ("A", 94, "B", 51),
        ("A", 37, "B", 52),
        ("A", 45, "B", 11),
    }

    assert observed_con_set == expected_con_set


def test_add_chain_from_segid(protprot_caprimodule):
    """Test replacing the chainID with segID."""
    tmp = tempfile.NamedTemporaryFile(delete=True)
    pdb_f = Path(golden_data, "protein_segid.pdb")
    shutil.copy(pdb_f, tmp.name)
    # this will replace-in-place
    protprot_caprimodule.add_chain_from_segid(tmp.name)

    with open(tmp.name) as fh:
        for line in fh:
            if line.startswith("ATOM"):
                assert line[21] == "A"


def test_calc_stats():
    """Test the calculation of statistics."""
    observed_mean, observed_std = calc_stats([2, 2, 4, 5])
    assert np.isclose(observed_mean, 3.25, atol=0.01)
    assert np.isclose(observed_std, 1.3, atol=0.01)


def test_capri_cluster_analysis(protprot_caprimodule, protprot_input_list):
    """Test the cluster analysis."""
    model1, model2 = protprot_input_list[0], protprot_input_list[1]
    model1.clt_rank, model2.clt_rank = 1, 2
    model1.clt_id, model2.clt_id = 1, 2
    model1.score, model2.score = 42.0, 50.0
    protprot_caprimodule.irmsd = 0.1
    protprot_caprimodule.fnat = 1.0
    protprot_caprimodule.lrmsd = 1.2
    protprot_caprimodule.ilrmsd = 4.3
    protprot_caprimodule.rmsd = (0.01,)
    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)
        capri_cluster_analysis(
            capri_list=[protprot_caprimodule, protprot_caprimodule],
            model_list=[model1, model2],
            output_fname="capri_clt.txt",
            clt_threshold=5,
            sort_key="score",
            sort_ascending=True,
            path=Path("."),
        )

        assert Path("capri_clt.txt").stat().st_size != 0

        observed_outf_l = read_capri_file("capri_clt.txt")
        expected_outf_l = [
            [
                "cluster_id",
                "n",
                "under_eval",
                "score",
                "score_std",
                "irmsd",
                "irmsd_std",
                "fnat",
                "fnat_std",
                "lrmsd",
                "lrmsd_std",
                "dockq",
                "dockq_std",
                "ilrmsd",
                "ilrmsd_std",
                "rmsd",
                "rmsd_std",
                "caprieval_rank",
            ],
            [
                "1",
                "1",
                "yes",
                "42.000",
                "0.000",
                "0.100",
                "0.000",
                "1.000",
                "0.000",
                "1.200",
                "0.000",
                "nan",
                "nan",
                "4.300",
                "0.000",
                "0.010",
                "0.000",
                "1",
            ],
            [
                "2",
                "1",
                "yes",
                "50.000",
                "0.000",
                "0.100",
                "0.000",
                "1.000",
                "0.000",
                "1.200",
                "0.000",
                "nan",
                "nan",
                "4.300",
                "0.000",
                "0.010",
                "0.000",
                "2",
            ],
        ]
        assert observed_outf_l == expected_outf_l

        # test sorting
        capri_cluster_analysis(
            capri_list=[protprot_caprimodule, protprot_caprimodule],
            model_list=[model1, model2],
            output_fname="capri_clt.txt",
            clt_threshold=5,
            sort_key="cluster_rank",
            sort_ascending=False,
            path=Path("."),
        )

        observed_outf_l = read_capri_file("capri_clt.txt")
        expected_outf_l = [
            [
                "cluster_id",
                "n",
                "under_eval",
                "score",
                "score_std",
                "irmsd",
                "irmsd_std",
                "fnat",
                "fnat_std",
                "lrmsd",
                "lrmsd_std",
                "dockq",
                "dockq_std",
                "rmsd",
                "rmsd_std",
                "caprieval_rank",
            ],
            [
                "2",
                "1",
                "yes",
                "50.000",
                "0.000",
                "0.100",
                "0.000",
                "1.000",
                "0.000",
                "1.200",
                "0.000",
                "nan",
                "nan",
                "0.010",
                "0.000",
                "2",
            ],
            [
                "1",
                "1",
                "yes",
                "42.000",
                "0.000",
                "0.100",
                "0.000",
                "1.000",
                "0.000",
                "1.200",
                "0.000",
                "nan",
                "nan",
                "0.010",
                "0.000",
                "1",
            ],
        ]


def test_single_chain_reference(protprot_onechain_ref_caprimodule, params):
    """Test correct values if reference has a single chain."""
    # fnat
    protprot_onechain_ref_caprimodule.calc_fnat()
    assert np.isnan(protprot_onechain_ref_caprimodule.fnat)
    # irmsd
    protprot_onechain_ref_caprimodule.calc_irmsd()
    assert np.isnan(protprot_onechain_ref_caprimodule.irmsd)
    # lrmsd
    protprot_onechain_ref_caprimodule.calc_lrmsd()
    assert np.isnan(protprot_onechain_ref_caprimodule.lrmsd)
    # ilrmsd
    protprot_onechain_ref_caprimodule.calc_ilrmsd()
    assert np.isnan(protprot_onechain_ref_caprimodule.ilrmsd)
    # dockq
    protprot_onechain_ref_caprimodule.calc_dockq()
    assert np.isnan(protprot_onechain_ref_caprimodule.dockq)


def test_single_chain_model(protprot_onechain_mod_caprimodule, params):
    """Test correct values if model has a single chain."""
    # fnat
    protprot_onechain_mod_caprimodule.calc_fnat()
    assert np.isclose(protprot_onechain_mod_caprimodule.fnat, 0.00, atol=0.01)
    # irmsd might be different than zero
    protprot_onechain_mod_caprimodule.calc_irmsd()
    assert np.isclose(protprot_onechain_mod_caprimodule.irmsd, 12.85, atol=0.01)
    # lrmsd
    protprot_onechain_mod_caprimodule.calc_lrmsd()
    assert np.isnan(protprot_onechain_mod_caprimodule.lrmsd)
    # ilrmsd
    protprot_onechain_mod_caprimodule.calc_ilrmsd()
    assert np.isnan(protprot_onechain_mod_caprimodule.ilrmsd)
    # dockq
    protprot_onechain_mod_caprimodule.calc_dockq()
    assert np.isnan(protprot_onechain_mod_caprimodule.dockq)


def test_get_previous_cns_step():
    "Test getting the previous CNS step."
    mock_steps = ["0_topoaa", "1_emscoring", "1_emscoring_fake", "2_clustfcc"]
    assert get_previous_cns_step(mock_steps, 2) == "1_emscoring"
    mock_steps_2 = ["bla", "test"]
    assert get_previous_cns_step(mock_steps_2, 1) is None


# TODO: Move to integration tests
def test_protprot_1bkd_swapped_chains(protprot_1bkd_caprimodule):
    """Test protein-protein l-rmsd and ilrmsd calculation with swapped chains."""
    protprot_1bkd_caprimodule.r_chain = "S"
    protprot_1bkd_caprimodule.l_chain = "R"
    protprot_1bkd_caprimodule.calc_lrmsd()
    assert np.isclose(protprot_1bkd_caprimodule.lrmsd, 19.21, atol=0.01)
    protprot_1bkd_caprimodule.calc_ilrmsd()
    assert np.isclose(protprot_1bkd_caprimodule.ilrmsd, 16.33, atol=0.01)


# TODO: Move to integration tests
def test_protdna_swapped_chains(protdna_caprimodule):
    """Test protein-dna l-rmsd and ilrmsd calculation with swapped chains."""
    protdna_caprimodule.r_chain = "B"
    protdna_caprimodule.l_chain = "A"
    protdna_caprimodule.calc_lrmsd()
    assert np.isclose(protdna_caprimodule.lrmsd, 4.81, atol=0.01)
    protdna_caprimodule.calc_ilrmsd()
    assert np.isclose(protdna_caprimodule.ilrmsd, 4.91, atol=0.01)


def test_capri_run(mocker):
    """???"""

    mock_get_align_func = mocker.Mock(
        return_value={"numbering": {1: 1, 2: 2}, "chain_dict": {1: 2}}
    )
    mocker.patch(
        "haddock.modules.analysis.caprieval.capri.get_align",
        return_value=mock_get_align_func,
    )

    mocker.patch.object(CAPRI, "_load_atoms", return_value=None)
    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)
        capri = CAPRI(
            identificator="test",
            model=Path("some-file"),
            path=Path("."),
            reference=Path("some-file"),
            params={
                "fnat": True,
                "irmsd": True,
                "lrmsd": True,
                "ilrmsd": True,
                "dockq": True,
                "global_rmsd": True,
                "allatoms": True,
                "receptor_chain": "A",
                "ligand_chains": ["B"],
                "alignment_method": None,
                "lovoalign_exec": None,
                "fnat_cutoff": 10,
                "irmsd_cutoff": 10,
            },
        )
        rand_fnat = random.random()
        mocker.patch.object(
            capri,
            "calc_fnat",
            side_effect=lambda cutoff: setattr(capri, "fnat", rand_fnat),
        )
        rand_irmsd = random.random()
        mocker.patch.object(
            capri,
            "calc_irmsd",
            side_effect=lambda cutoff: setattr(capri, "irmsd", rand_irmsd),
        )
        rand_lrmsd = random.random()
        mocker.patch.object(
            capri, "calc_lrmsd", side_effect=lambda: setattr(capri, "lrmsd", rand_lrmsd)
        )
        rand_ilrmsd = random.random()
        mocker.patch.object(
            capri,
            "calc_ilrmsd",
            side_effect=lambda cutoff: setattr(capri, "ilrmsd", rand_ilrmsd),
        )
        rand_dockq = random.random()
        mocker.patch.object(
            capri, "calc_dockq", side_effect=lambda: setattr(capri, "dockq", rand_dockq)
        )

        rand_global_rmsd = random.random()
        mocker.patch.object(
            capri,
            "calc_global_rmsd",
            side_effect=lambda: setattr(capri, "rmsd", rand_global_rmsd),
        )


        capri.run()

        # The only logic to be tested is if the methods are called
        assert capri.fnat == pytest.approx(rand_fnat)
        assert capri.irmsd == pytest.approx(rand_irmsd)
        assert capri.lrmsd == pytest.approx(rand_lrmsd)
        assert capri.ilrmsd == pytest.approx(rand_ilrmsd)
        assert capri.dockq == pytest.approx(rand_dockq)
        assert capri.rmsd == pytest.approx(rand_global_rmsd)


def test_rank_according_to_score():
    """test ranking according to score."""
    data = {
        1: {
            "score": 3.0,
            "caprieval_rank": 99999,
        },
        2: {
            "score": 2.0,
            "caprieval_rank": 99999,
        },
        3: {
            "score": 1.0,
            "caprieval_rank": 99999,
        },
    }

    ranked_data = rank_according_to_score(
        data=data, sort_key="score", sort_ascending=True
    )

    assert ranked_data[1]["caprieval_rank"] == 1
    assert ranked_data[2]["caprieval_rank"] == 2
    assert ranked_data[3]["caprieval_rank"] == 3


def test_extract_data_from_capri_class(mocker):
    """???"""

    mocker.patch(
        "haddock.modules.analysis.caprieval.capri.write_nested_dic_to_file",
        return_value=None,
    )
    mocker.patch.object(CAPRI, "_load_atoms", return_value=None)

    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)

        c = CAPRI(
            path=Path("."),
            identificator="42",
            model=Path("anything"),
            reference=Path("anything"),
            params={
                "allatoms": True,
                "receptor_chain": "X",
                "ligand_chains": ["X"],
            },
        )

        # Create random entries
        random_energy = random.random()
        random_model = PDBFile(
            file_name=str(uuid.uuid4()),
            score=42,
            unw_energies={"energy": random_energy},
        )

        random_clt_id = random.randint(0, 100)
        random_clt_rank = random.randint(0, 100)
        random_clt_model_rank = random.randint(0, 100)
        random_md5 = str(uuid.uuid4())
        random_score = random.random()
        random_irmsd = random.random()
        random_fnat = random.random()
        random_lrmsd = random.random()
        random_ilrmsd = random.random()
        random_dockq = random.random()
        random_rmsd = random.random()

        c.model = random_model
        c.model.clt_id = random_clt_id
        c.model.clt_rank = random_clt_rank
        c.model.clt_model_rank = random_clt_model_rank
        c.md5 = random_md5
        c.score = random_score
        c.irmsd = random_irmsd
        c.fnat = random_fnat
        c.lrmsd = random_lrmsd
        c.ilrmsd = random_ilrmsd
        c.dockq = random_dockq
        c.rmsd = random_rmsd

        observed_data = extract_data_from_capri_class(
            capri_objects=[c],
            sort_key="score",
            sort_ascending=True,
            output_fname=Path(""),
        )

        assert observed_data is not None

        assert observed_data[1]["model"] == random_model
        assert observed_data[1]["md5"] == random_md5
        assert observed_data[1]["score"] == random_score
        assert observed_data[1]["irmsd"] == random_irmsd
        assert observed_data[1]["fnat"] == random_fnat
        assert observed_data[1]["lrmsd"] == random_lrmsd
        assert observed_data[1]["ilrmsd"] == random_ilrmsd
        assert observed_data[1]["dockq"] == random_dockq
        assert observed_data[1]["rmsd"] == random_rmsd
        assert observed_data[1]["energy"] == random_energy
        assert observed_data[1]["cluster_id"] == random_clt_id
        assert observed_data[1]["cluster_ranking"] == random_clt_rank
        assert observed_data[1]["model-cluster_ranking"] == random_clt_model_rank
