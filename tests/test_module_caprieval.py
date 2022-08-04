"""Test the CAPRI module."""
import os
import shutil
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.caprieval.capri import (
    CAPRI,
    calc_stats,
    capri_cluster_analysis,
    rearrange_ss_capri_output,
    )

from . import golden_data


def round_two_dec(value):
    """Round a value to two digits."""
    return round(value, 2)


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
def params():
    return {
        "receptor_chain": "A", "ligand_chain": "B", "aln_method": "sequence"
        }


@pytest.fixture
def protdna_caprimodule(protdna_input_list, params):
    """Protein-DNA CAPRI module."""
    reference = protdna_input_list[0].rel_path
    model = protdna_input_list[1].rel_path
    capri = CAPRI(
        identificator=42,
        reference=reference,
        model=model,
        path=golden_data,
        params=params,
        )

    yield capri

    remove_aln_files(capri)


@pytest.fixture
def protlig_caprimodule(protlig_input_list, params):
    """Protein-Ligand CAPRI module."""
    reference = protlig_input_list[0].rel_path
    model = protlig_input_list[1].rel_path
    capri = CAPRI(
        identificator=42,
        reference=reference,
        model=model,
        path=golden_data,
        params=params,
        )

    yield capri

    remove_aln_files(capri)


@pytest.fixture
def protprot_caprimodule(protprot_input_list, params):
    """Protein-Protein CAPRI module."""
    reference = protprot_input_list[0].rel_path
    model = protprot_input_list[1]
    capri = CAPRI(
        identificator=42,
        reference=reference,
        model=model,
        path=golden_data,
        params=params,
        )

    yield capri

    remove_aln_files(capri)


@pytest.fixture
def protprot_caprimodule_parallel(protprot_input_list):
    """Protein-Protein CAPRI module."""
    reference = protprot_input_list[0].rel_path
    model = protprot_input_list[1].rel_path
    capri = CAPRI(
        reference=reference,
        model=model,
        receptor_chain="A",
        ligand_chain="B",
        aln_method="sequence",
        path=golden_data,
        identificator=0,
        core_model_idx=0
        )

    yield capri

    remove_aln_files(capri)


def test_protprot_irmsd(protprot_caprimodule):
    """Test protein-protein i-rmsd calculation."""
    protprot_caprimodule.calc_irmsd()
    assert round_two_dec(protprot_caprimodule.irmsd) == 7.38


def test_protprot_lrmsd(protprot_caprimodule):
    """Test protein-protein l-rmsd calculation."""
    protprot_caprimodule.calc_lrmsd()
    assert round_two_dec(protprot_caprimodule.lrmsd) == 15.9


def test_protprot_ilrmsd(protprot_caprimodule):
    """Test protein-protein i-l-rmsd calculation."""
    protprot_caprimodule.calc_ilrmsd()
    assert round_two_dec(protprot_caprimodule.ilrmsd) == 9.67


def test_protprot_fnat(protprot_caprimodule):
    """Test protein-protein fnat calculation."""
    protprot_caprimodule.calc_fnat()
    assert round_two_dec(protprot_caprimodule.fnat) == 0.05


def test_protlig_irmsd(protlig_caprimodule):
    """Test protein-ligand i-rmsd calculation."""
    protlig_caprimodule.calc_irmsd()
    assert round_two_dec(protlig_caprimodule.irmsd) == 0.22


def test_protlig_lrmsd(protlig_caprimodule):
    """Test protein-ligand l-rmsd calculation."""
    protlig_caprimodule.calc_lrmsd()
    assert round_two_dec(protlig_caprimodule.lrmsd) == 0.51


def test_protlig_ilrmsd(protlig_caprimodule):
    """Test protein-ligand i-l-rmsd calculation."""
    protlig_caprimodule.calc_ilrmsd()
    assert round_two_dec(protlig_caprimodule.ilrmsd) == 0.5


def test_protlig_fnat(protlig_caprimodule):
    """Test protein-ligand fnat calculation."""
    protlig_caprimodule.calc_fnat()
    assert round_two_dec(protlig_caprimodule.fnat) == 1.0


def test_protdna_irmsd(protdna_caprimodule):
    """Test protein-dna i-rmsd calculation."""
    protdna_caprimodule.calc_irmsd()
    assert round_two_dec(protdna_caprimodule.irmsd) == 2.05


def test_protdna_lrmsd(protdna_caprimodule):
    """Test protein-dna l-rmsd calculation."""
    protdna_caprimodule.calc_lrmsd()
    assert round_two_dec(protdna_caprimodule.lrmsd) == 4.19


def test_protdna_ilrmsd(protdna_caprimodule):
    """Test protein-dna i-l-rmsd calculation."""
    protdna_caprimodule.calc_ilrmsd()
    assert round_two_dec(protdna_caprimodule.ilrmsd) == 1.89


def test_protdna_fnat(protdna_caprimodule, protdna_input_list):
    """Test protein-dna fnat calculation."""
    protdna_caprimodule.calc_fnat()
    assert round_two_dec(protdna_caprimodule.fnat) == 0.49


def test_make_output(protprot_caprimodule):
    """Test the writing of capri.tsv file."""
    protprot_caprimodule.model.clt_id = 1
    protprot_caprimodule.model.clt_rank = 1
    protprot_caprimodule.model.clt_model_rank = 10

    protprot_caprimodule.make_output()

    ss_fname = Path(
        protprot_caprimodule.path,
        f"capri_ss_{protprot_caprimodule.identificator}.tsv"
        )

    assert ss_fname.stat().st_size != 0

    # remove the model column since its name will depend on where we are running
    #  the test
    observed_outf_l = [e.split()[1:] for e in open(
        ss_fname).readlines() if not e.startswith('#')]
    expected_outf_l = [
        ['md5', 'caprieval_rank', 'score', 'irmsd', 'fnat', 'lrmsd', 'ilrmsd',
         'dockq', 'cluster-id', 'cluster-ranking', 'model-cluster-ranking'],
        ['-', '-', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', '1', '1', '10'], ]

    assert observed_outf_l == expected_outf_l

    os.unlink(ss_fname)


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


def test_rearrange_ss_capri_output():
    """Test rearranging the capri output."""
    with open(f"{golden_data}/capri_ss_1.tsv", 'w') as fh:
        fh.write(
            "model	caprieval_rank	score	irmsd	fnat	lrmsd	ilrmsd	"
            "dockq	cluster-id	cluster-ranking	"
            "model-cluster-ranking" + os.linesep)
        fh.write(
            "../1_emscoring/emscoring_909.pdb	1	-424.751	0.000	"
            "1.000	0.000	0.000	1.000	-	-	-" + os.linesep)
    rearrange_ss_capri_output(
        'capri_ss.txt',
        output_count=1,
        sort_key="score",
        sort_ascending=True,
        path=golden_data)

    assert Path('capri_ss.txt').stat().st_size != 0
    Path('capri_ss.txt').unlink()


def test_calc_stats():
    """Test the calculation of statistics."""
    observed_mean, observed_std = calc_stats([2, 2, 4, 5])
    assert round_two_dec(observed_mean) == 3.25
    assert round_two_dec(observed_std) == 1.3


def test_capri_cluster_analysis(protprot_caprimodule, protprot_input_list):
    """Test the cluster analysis."""
    model = protprot_input_list[0]
    model.clt_rank = 1
    model.clt_id = 1
    model.score = 42.0
    protprot_caprimodule.irmsd = 0.1
    protprot_caprimodule.fnat = 1.0
    protprot_caprimodule.lrmsd = 1.2
    protprot_caprimodule.ilrmsd = 4.3
    capri_cluster_analysis(
        capri_list=[protprot_caprimodule],
        model_list=[model],
        output_fname="capri_clt.txt",
        clt_threshold=5,
        sort_key="score",
        sort_ascending=True,
        path=Path("."))

    assert Path('capri_clt.txt').stat().st_size != 0

    Path('capri_clt.txt').unlink()


def test_check_chains(protprot_caprimodule):
    """Test correct checking of chains."""
    obs_ch = [["A", "C"],
              ["A", "B"],
              ["S", "E", "B", "A"],
              ["C", "D"]]
    
    # assuming exp chains are A and B
    exp_ch = [["A", "C"],
              ["A", "B"],
              ["A", "B"],
              ["C", "D"]]

    for n in range(len(obs_ch)):
        obs_r_chain, obs_l_chain = protprot_caprimodule.check_chains(obs_ch[n])
        exp_r_chain, exp_l_chain = exp_ch[n][0], exp_ch[n][1]
        assert obs_r_chain == exp_r_chain
        assert obs_l_chain == exp_l_chain
