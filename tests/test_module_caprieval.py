"""Test the CAPRI module."""
import os
import shutil
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.caprieval.capri import CAPRI

from . import golden_data


def round_two_dec(dic):
    """Round a rms dictionary to two digits."""
    return dict((k, round(dic[k], 2)) for k in dic)


def remove_aln_files(class_name):
    """Remove intermediary alignment files."""
    file_l = [Path(class_name.path, 'blosum62.izone'),
              Path(class_name.path, 'blosum62_A.aln'),
              Path(class_name.path, 'blosum62_B.aln')]
    for f in file_l:
        if f.exists():
            os.unlink(f)


@pytest.fixture
def protprot_input_list():
    """Prot-prot input."""
    return [
        PDBFile(Path(golden_data, "protprot_complex_1.pdb"), path=golden_data),
        PDBFile(Path(golden_data, "protprot_complex_2.pdb"), path=golden_data)
        ]


@pytest.fixture
def protdna_input_list():
    """Prot-DNA input."""
    return [
        PDBFile(Path(golden_data, "protdna_complex_1.pdb"), path=golden_data),
        PDBFile(Path(golden_data, "protdna_complex_2.pdb"), path=golden_data)
        ]


@pytest.fixture
def protlig_input_list():
    """Protein-Ligand input."""
    return [
        PDBFile(Path(golden_data, "protlig_complex_1.pdb"), path=golden_data),
        PDBFile(Path(golden_data, "protlig_complex_2.pdb"), path=golden_data),
        ]


@pytest.fixture
def protdna_caprimodule(protdna_input_list):
    """Protein-DNA CAPRI module."""
    ref = protdna_input_list[0].rel_path
    capri = CAPRI(
        reference=ref,
        model_list=protdna_input_list,
        receptor_chain="A",
        ligand_chain="B",
        aln_method="sequence",
        path=golden_data,
        )

    yield capri

    remove_aln_files(capri)


@pytest.fixture
def protlig_caprimodule(protlig_input_list):
    """Protein-Ligand CAPRI module."""
    ref = protlig_input_list[0].rel_path
    capri = CAPRI(
        reference=ref,
        model_list=protlig_input_list,
        receptor_chain="A",
        ligand_chain="B",
        aln_method="sequence",
        path=golden_data,
        )

    yield capri

    remove_aln_files(capri)


@pytest.fixture
def protprot_caprimodule(protprot_input_list):
    """Protein-Protein CAPRI module."""
    ref = protprot_input_list[0].rel_path
    capri = CAPRI(
        reference=ref,
        model_list=protprot_input_list,
        receptor_chain="A",
        ligand_chain="B",
        aln_method="sequence",
        path=golden_data,
        )

    yield capri

    remove_aln_files(capri)


@pytest.fixture
def protprot_caprimodule_parallel(protprot_input_list):
    """Protein-Protein CAPRI module."""
    ref = protprot_input_list[0].rel_path
    capri = CAPRI(
        reference=ref,
        model_list=protprot_input_list,
        receptor_chain="A",
        ligand_chain="B",
        aln_method="sequence",
        path=golden_data,
        core=0,
        core_model_idx=0
        )

    yield capri

    remove_aln_files(capri)


def test_protprot_irmsd(protprot_caprimodule, protprot_input_list):
    """Test protein-protein i-rmsd calculation."""
    observed_irmsd_dic = protprot_caprimodule.irmsd()
    # round the irmsd to 2 digits so we can compare
    observed_irmsd_dic = round_two_dec(observed_irmsd_dic)

    expected_irmsd_dic = {
        protprot_input_list[0]: 0.0,
        protprot_input_list[1]: 7.92,
        }

    assert observed_irmsd_dic == expected_irmsd_dic


def test_protprot_lrmsd(protprot_caprimodule, protprot_input_list):
    """Test protein-protein l-rmsd calculation."""
    observed_lrmsd_dic = protprot_caprimodule.lrmsd()
    observed_lrmsd_dic = round_two_dec(observed_lrmsd_dic)

    expected_lrmsd_dic = {
        protprot_input_list[0]: 0.0,
        protprot_input_list[1]: 15.85,
        }

    assert observed_lrmsd_dic == expected_lrmsd_dic


def test_protprot_ilrmsd(protprot_caprimodule, protprot_input_list):
    """Test protein-protein i-l-rmsd calculation."""
    observed_ilrmsd_dic = protprot_caprimodule.ilrmsd()
    observed_ilrmsd_dic = round_two_dec(observed_ilrmsd_dic)

    expected_ilrmsd_dic = {
        protprot_input_list[0]: 0.0,
        protprot_input_list[1]: 9.66,
        }

    assert observed_ilrmsd_dic == expected_ilrmsd_dic


def test_protprot_fnat(protprot_caprimodule, protprot_input_list):
    """Test protein-protein fnat calculation."""
    observed_fnat_dic = protprot_caprimodule.fnat()
    observed_fnat_dic = round_two_dec(observed_fnat_dic)

    expected_fnat_dic = {
        protprot_input_list[0]: 1.0,
        protprot_input_list[1]: 0.05,
        }

    assert observed_fnat_dic == expected_fnat_dic


def test_protlig_irmsd(protlig_caprimodule, protlig_input_list):
    """Test protein-ligand i-rmsd calculation."""
    observed_irmsd_dic = protlig_caprimodule.irmsd()
    # round the irmsd to 2 digits so we can compare
    observed_irmsd_dic = round_two_dec(observed_irmsd_dic)

    expected_irmsd_dic = {
        protlig_input_list[0]: 0.0,
        protlig_input_list[1]: 0.22,
        }

    assert observed_irmsd_dic == expected_irmsd_dic


def test_protlig_lrmsd(protlig_caprimodule, protlig_input_list):
    """Test protein-ligand l-rmsd calculation."""
    observed_lrmsd_dic = protlig_caprimodule.lrmsd()
    # round the irmsd to 2 digits so we can compare
    observed_lrmsd_dic = round_two_dec(observed_lrmsd_dic)

    expected_lrmsd_dic = {
        protlig_input_list[0]: 0.0,
        protlig_input_list[1]: 0.51,
        }

    assert observed_lrmsd_dic == expected_lrmsd_dic


def test_protlig_ilrmsd(protlig_caprimodule, protlig_input_list):
    """Test protein-ligand i-l-rmsd calculation."""
    observed_ilrmsd_dic = protlig_caprimodule.ilrmsd()
    # round the irmsd to 2 digits so we can compare
    observed_ilrmsd_dic = round_two_dec(observed_ilrmsd_dic)

    expected_ilrmsd_dic = {
        protlig_input_list[0]: 0.0,
        protlig_input_list[1]: 0.5,
        }

    assert observed_ilrmsd_dic == expected_ilrmsd_dic


def test_protlig_fnat(protlig_caprimodule, protlig_input_list):
    """Test protein-ligand fnat calculation."""
    observed_fnat_dic = protlig_caprimodule.fnat()
    # round the irmsd to 2 digits so we can compare
    observed_fnat_dic = round_two_dec(observed_fnat_dic)

    expected_fnat_dic = {protlig_input_list[0]: 1.0, protlig_input_list[1]: 1.0}

    assert observed_fnat_dic == expected_fnat_dic


def test_protdna_irmsd(protdna_caprimodule, protdna_input_list):
    """Test protein-dna i-rmsd calculation."""
    observed_irmsd_dic = protdna_caprimodule.irmsd()
    # round the irmsd to 2 digits so we can compare
    observed_irmsd_dic = round_two_dec(observed_irmsd_dic)

    expected_irmsd_dic = {
        protdna_input_list[0]: 0.0,
        protdna_input_list[1]: 2.05,
        }

    assert observed_irmsd_dic == expected_irmsd_dic


def test_protdna_lrmsd(protdna_caprimodule, protdna_input_list):
    """Test protein-dna l-rmsd calculation."""
    observed_lrmsd_dic = protdna_caprimodule.lrmsd()
    # round the irmsd to 2 digits so we can compare
    observed_lrmsd_dic = round_two_dec(observed_lrmsd_dic)

    expected_lrmsd_dic = {
        protdna_input_list[0]: 0.0,
        protdna_input_list[1]: 4.19,
        }

    assert observed_lrmsd_dic == expected_lrmsd_dic


def test_protdna_ilrmsd(protdna_caprimodule, protdna_input_list):
    """Test protein-dna i-l-rmsd calculation."""
    observed_ilrmsd_dic = protdna_caprimodule.ilrmsd()
    # round the irmsd to 2 digits so we can compare
    observed_ilrmsd_dic = round_two_dec(observed_ilrmsd_dic)

    expected_ilrmsd_dic = {
        protdna_input_list[0]: 0.0,
        protdna_input_list[1]: 1.89,
        }

    assert observed_ilrmsd_dic == expected_ilrmsd_dic


def test_protdna_fnat(protdna_caprimodule, protdna_input_list):
    """Test protein-dna fnat calculation."""
    observed_fnat_dic = protdna_caprimodule.fnat()
    # round the irmsd to 2 digits so we can compare
    observed_fnat_dic = round_two_dec(observed_fnat_dic)

    expected_fnat_dic = {
        protdna_input_list[0]: 1.0,
        protdna_input_list[1]: 0.49,
        }

    assert observed_fnat_dic == expected_fnat_dic


def test_output(protprot_caprimodule):
    """Test the writing of capri.tsv file."""
    factor = 1
    clt_id = 1
    for m in protprot_caprimodule.model_list:
        protprot_caprimodule.irmsd_dic[m] = 0.111 * factor
        protprot_caprimodule.fnat_dic[m] = 0.333 * factor
        protprot_caprimodule.lrmsd_dic[m] = 0.444 * factor
        protprot_caprimodule.ilrmsd_dic[m] = 0.555 * factor
        m.clt_id = clt_id
        clt_id += 1
        factor += 1.5

    sortby_key = "fnat"
    sort_ascending = True
    clt_threshold = 1
    protprot_caprimodule.output(
        clt_threshold,
        sortby_key=sortby_key,
        sort_ascending=sort_ascending,
        )

    ss_fname = Path(protprot_caprimodule.path, "capri_ss.tsv")
    clt_fname = Path(protprot_caprimodule.path, "capri_clt.tsv")

    assert ss_fname.stat().st_size != 0
    assert clt_fname.stat().st_size != 0

    # remove the model column since its name will depend on where we are running
    #  the test
    observed_outf_l = [e.split()[1:] for e in open(
        ss_fname).readlines() if not e.startswith('#')]
    expected_outf_l = [
        ['caprieval_rank', 'score', 'irmsd', 'fnat', 'lrmsd', 'ilrmsd',
         'cluster-id', 'cluster-ranking', 'model-cluster-ranking'],
        ['1', 'nan', '0.111', '0.333', '0.444', '0.555', '1', '-', '-'],
        ['2', 'nan', '0.278', '0.833', '1.110', '1.388', '2', '-', '-']]

    assert observed_outf_l == expected_outf_l

    observed_outf_l = [e.split() for e in open(
        clt_fname).readlines() if not e.startswith('#')]
    expected_outf_l = [
        ['cluster_rank', 'cluster_id', 'n', 'under_eval', 'score', 'score_std',
         'irmsd', 'irmsd_std', 'fnat', 'fnat_std', 'lrmsd', 'lrmsd_std',
         'dockqn', 'dockq_std', 'caprieval_rank'],
        ['-', '1', '1', '-', 'nan', 'nan', '0.111', '0.000', '0.333', '0.000',
         '0.444', '0.000', 'nan', 'nan', '1'],
        ['-', '2', '1', '-', 'nan', 'nan', '0.278', '0.000', '0.833', '0.000',
         '1.110', '0.000', 'nan', 'nan', '2']]

    assert observed_outf_l == expected_outf_l

    os.unlink(ss_fname)
    os.unlink(clt_fname)


def test_output_parallel(protprot_caprimodule_parallel):
    """Test the writing of capri.tsv file."""
    factor = 1
    clt_id = 1
    for m in protprot_caprimodule_parallel.model_list:
        protprot_caprimodule_parallel.irmsd_dic[m] = 0.111 * factor
        protprot_caprimodule_parallel.fnat_dic[m] = 0.333 * factor
        protprot_caprimodule_parallel.lrmsd_dic[m] = 0.444 * factor
        protprot_caprimodule_parallel.ilrmsd_dic[m] = 0.555 * factor
        m.clt_id = clt_id
        clt_id += 1
        factor += 1.5

    sortby_key = "fnat"
    sort_ascending = True
    clt_threshold = 1
    protprot_caprimodule_parallel.output(
        clt_threshold,
        sortby_key=sortby_key,
        sort_ascending=sort_ascending,
        )

    # check that the parallel files are present
    ss_fname_0 = Path(protprot_caprimodule_parallel.path, "capri_ss_0.tsv")
    clt_fname_0 = Path(protprot_caprimodule_parallel.path, "capri_clt_0.tsv")

    assert ss_fname_0.stat().st_size != 0
    assert clt_fname_0.stat().st_size != 0

    # replicate the previous task, on files *_0.tsv, to ensure consistency
    observed_outf_l = [e.split()[1:] for e in open(
        ss_fname_0).readlines() if not e.startswith('#')]
    expected_outf_l = [
        ['caprieval_rank', 'score', 'irmsd', 'fnat', 'lrmsd', 'ilrmsd',
         'cluster-id', 'cluster-ranking', 'model-cluster-ranking'],
        ['1', 'nan', '0.111', '0.333', '0.444', '0.555', '1', '-', '-'],
        ['2', 'nan', '0.278', '0.833', '1.110', '1.388', '2', '-', '-']]

    assert observed_outf_l == expected_outf_l

    observed_outf_l = [e.split() for e in open(
        clt_fname_0).readlines() if not e.startswith('#')]
    expected_outf_l = [
        ['cluster_rank', 'cluster_id', 'n', 'under_eval', 'score', 'score_std',
         'irmsd', 'irmsd_std', 'fnat', 'fnat_std', 'lrmsd', 'lrmsd_std',
         'dockqn', 'dockq_std', 'caprieval_rank'],
        ['-', '1', '1', '-', 'nan', 'nan', '0.111', '0.000', '0.333', '0.000',
         '0.444', '0.000', 'nan', 'nan', '1'],
        ['-', '2', '1', '-', 'nan', 'nan', '0.278', '0.000', '0.833', '0.000',
         '1.110', '0.000', 'nan', 'nan', '2']]

    assert observed_outf_l == expected_outf_l

    os.unlink(ss_fname_0)
    os.unlink(clt_fname_0)


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

    assert observed_interface == expected_interface


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

    assert observed_interface == expected_interface


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

    assert observed_interface == expected_interface


def test_load_contacts(protprot_caprimodule, protprot_input_list):
    """Test loading contacts."""
    protprot_complex = protprot_input_list[0]
    observed_con_set = protprot_caprimodule.load_contacts(
        protprot_complex, cutoff=5.0
        )
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


# def test_write_coord_dic():
#     pass


# def test_write_coords():
#     pass


# def test_write_pymol_viz():
#     pass
