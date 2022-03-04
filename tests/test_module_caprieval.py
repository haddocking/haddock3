"""Test the CAPRI module."""
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.caprieval.capri import (
    CAPRI,
    align_seq,
    dump_as_izone,
    get_align,
    get_atoms,
    make_range,
    pdb2fastadic,
    )

from . import golden_data


def array_to_list(np_array):
    """Transform a numpy array in a nested list."""
    return [list(e) for e in np_array]


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
    rankby_key = "irmsd"
    rank_ascending = True
    clt_threshold = 1
    protprot_caprimodule.output(
        clt_threshold,
        sortby_key=sortby_key,
        sort_ascending=sort_ascending,
        rankby_key=rankby_key,
        rank_ascending=rank_ascending,
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
        ['caprieval_rank', 'cluster_rank', 'cluster_id', 'n', 'under_eval',
         'score', 'score_std', 'irmsd', 'irmsd_std', 'fnat', 'fnat_std',
         'lrmsd', 'lrmsd_std', 'dockqn', 'dockq_std'],
        ['1', '-', '1', '1', '-', 'nan', 'nan', '0.111', '0.000', '0.333',
         '0.000', '0.444', '0.000', 'nan', 'nan'],
        ['2', '-', '2', '1', '-', 'nan', 'nan', '0.278', '0.000', '0.833',
         '0.000', '1.110', '0.000', 'nan', 'nan']]

    assert observed_outf_l == expected_outf_l

    os.unlink(ss_fname)
    os.unlink(clt_fname)


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


def test_calc_rmsd(protprot_caprimodule):
    """Test the RMSD calculation."""
    V = [
        [0.4975639180542826, 2.4634453459517913, 13.502915985158635],
        [11.227311917422782, 16.269820879072906, -3.962511851240976],
        [-7.0306634127649135, 4.226046643850031, -8.654935533437529],
        [-1.592886678968348, -10.874055594404654, -5.820287656251305],
        [-4.372921399498893, 2.094055858210638, -8.93135311766362],
        [-1.2876205748591463, -1.4125224802187486, 15.515992917804834],
        [6.43227627252424, 17.468089367358402, 11.04719510106607],
        [5.88420286284558, 9.382362089324166, -5.708004059876668],
        [4.238076687025332, 3.8305358597829944, -5.857279410346993],
        [5.145346893131023, -2.7233805318567046, 7.910771966100926],
        ]
    W = [
        [-2.2927341772151664, 1.3226734177215107, 12.252610126582304],
        [8.353265822784834, 10.995673417721513, -7.888389873417699],
        [-12.224734177215169, -8.316326582278489, 6.479610126582301],
        [-6.9257341772151655, 6.939673417721508, -1.281389873417698],
        [-11.438734177215167, -5.325326582278489, 8.2086101265823],
        [-3.1317341772151686, -1.6973265822784924, 15.408610126582303],
        [-0.2497341772151671, 15.53367341772151, 4.4986101265822995],
        [5.965265822784836, 2.219673417721509, -7.501389873417699],
        [5.648265822784836, -3.0373265822784887, -4.8293898734177],
        [5.240265822784835, -3.3833265822784924, 9.9936101265823],
        ]

    rmsd = protprot_caprimodule.calc_rmsd(V, W)

    assert round(rmsd, 2) == 12.02


def test_kabsch(protprot_caprimodule):
    """Test the Kabsch algorithm."""
    P = [
        [-3.811237974683542, -0.12069367088607574, 7.200868354430379],
        [-5.913237974683542, 8.536306329113923, 6.286868354430379],
        [12.839762025316457, 8.507306329113923, 4.5878683544303795],
        [4.856762025316458, -3.7666936708860748, -10.57713164556962],
        [-1.9222379746835419, -4.171693670886076, -6.596131645569622],
        [-12.186237974683543, -4.807693670886075, -10.77713164556962],
        [0.3087620253164576, -7.500693670886075, 8.632868354430379],
        [-0.08023797468354221, -7.218693670886075, 2.9678683544303794],
        [9.232762025316458, 18.419306329113926, 8.949868354430379],
        [-1.7482379746835424, -8.395693670886075, 7.604868354430379],
        ]

    Q = [
        [6.0952658227848495, -5.630326582278489, 8.033610126582303],
        [0.5082658227848498, 1.733673417721512, 10.833610126582304],
        [7.42726582278485, 7.67067341772151, -5.757389873417694],
        [-12.59373417721515, -6.5263265822784895, 5.217610126582304],
        [-6.309734177215148, -4.75632658227849, 2.355610126582306],
        [-3.6647341772151485, -8.990326582278488, -8.014389873417695],
        [11.789265822784852, -10.62332658227849, 3.966610126582303],
        [6.895265822784854, -12.33832658227849, 1.2686101265823062],
        [5.7642658227848536, 16.110673417721507, 1.3446101265823032],
        [10.306265822784848, -12.226326582278489, 5.028610126582304],
        ]

    observed_U = protprot_caprimodule.kabsch(P, Q)

    expected_U = [
        [0.38845040189428726, 0.38160307568742435, -0.8387403518932808],
        [-0.20700756088693778, 0.9230933741949091, 0.32410876608493205],
        [0.8979165634978604, 0.047725414019726436, 0.43757071411697324],
        ]

    np.testing.assert_allclose(np.asarray(expected_U), observed_U)


def test_get_atoms():
    """Test the identification of atoms."""
    pdb_list = [
        Path(golden_data, "protein.pdb"),
        Path(golden_data, "dna.pdb"),
        Path(golden_data, "ligand.pdb"),
        ]
    observed_atom_dic = get_atoms(pdb_list)
    expected_atom_dic = {
        "ALA": ["C", "N", "CA", "O"],
        "ARG": ["C", "N", "CA", "O"],
        "ASN": ["C", "N", "CA", "O"],
        "ASP": ["C", "N", "CA", "O"],
        "CYS": ["C", "N", "CA", "O"],
        "GLN": ["C", "N", "CA", "O"],
        "GLU": ["C", "N", "CA", "O"],
        "GLY": ["C", "N", "CA", "O"],
        "HIS": ["C", "N", "CA", "O"],
        "ILE": ["C", "N", "CA", "O"],
        "LEU": ["C", "N", "CA", "O"],
        "LYS": ["C", "N", "CA", "O"],
        "MET": ["C", "N", "CA", "O"],
        "PHE": ["C", "N", "CA", "O"],
        "PRO": ["C", "N", "CA", "O"],
        "SER": ["C", "N", "CA", "O"],
        "THR": ["C", "N", "CA", "O"],
        "TRP": ["C", "N", "CA", "O"],
        "TYR": ["C", "N", "CA", "O"],
        "VAL": ["C", "N", "CA", "O"],
        "DA": [
            "C5",
            "N9",
            "N2",
            "C8",
            "O2",
            "N4",
            "N7",
            "C7",
            "N1",
            "N6",
            "C2",
            "O4",
            "C6",
            "N3",
            "C4",
            "O6",
            ],
        "DC": [
            "C5",
            "N9",
            "N2",
            "C8",
            "O2",
            "N4",
            "N7",
            "C7",
            "N1",
            "N6",
            "C2",
            "O4",
            "C6",
            "N3",
            "C4",
            "O6",
            ],
        "DT": [
            "C5",
            "N9",
            "N2",
            "C8",
            "O2",
            "N4",
            "N7",
            "C7",
            "N1",
            "N6",
            "C2",
            "O4",
            "C6",
            "N3",
            "C4",
            "O6",
            ],
        "DG": [
            "C5",
            "N9",
            "N2",
            "C8",
            "O2",
            "N4",
            "N7",
            "C7",
            "N1",
            "N6",
            "C2",
            "O4",
            "C6",
            "N3",
            "C4",
            "O6",
            ],
        "G39": [
            "C1",
            "O1A",
            "O1B",
            "C2",
            "C3",
            "C4",
            "N4",
            "C5",
            "N5",
            "C6",
            "C7",
            "O7",
            "C8",
            "C9",
            "C10",
            "O10",
            "C11",
            "C81",
            "C82",
            "C91",
            ],
        }

    assert observed_atom_dic == expected_atom_dic


def test_centroid(protprot_caprimodule):
    """Test the centroid calculation."""
    X = [
        [6.0952658227848495, -5.630326582278489, 8.033610126582303],
        [0.5082658227848498, 1.733673417721512, 10.833610126582304],
        [7.42726582278485, 7.67067341772151, -5.757389873417694],
        ]

    observed_centroid = protprot_caprimodule.centroid(X)
    observed_centroid = list(observed_centroid)
    expected_centroid = [
        4.6769324894515165,
        1.2580067510548443,
        4.369943459915638,
        ]

    assert observed_centroid == expected_centroid


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


def test_load_coords(protprot_caprimodule):
    """Test the loading of coordinates."""
    # pdb_f = protprot_input_list[0]
    pdb_f = Path(golden_data, "protein.pdb")
    (
        observed_coord_dic,
        observed_chain_ranges,
        ) = protprot_caprimodule.load_coords(pdb_f)
    observed_keys = list(observed_coord_dic.keys())
    expected_keys = [
        ("B", 1, "C"),
        ("B", 1, "O"),
        ("B", 1, "N"),
        ("B", 1, "CA"),
        ("B", 2, "N"),
        ("B", 2, "CA"),
        ("B", 2, "C"),
        ("B", 2, "O"),
        ("B", 3, "N"),
        ("B", 3, "CA"),
        ("B", 3, "C"),
        ("B", 3, "O"),
        ("B", 4, "N"),
        ("B", 4, "CA"),
        ("B", 4, "C"),
        ("B", 4, "O"),
        ("B", 5, "N"),
        ("B", 5, "CA"),
        ("B", 5, "C"),
        ("B", 5, "O"),
        ]

    assert observed_keys == expected_keys

    observed_coords = array_to_list(observed_coord_dic.values())
    expected_coords = [
        [2.76, 8.901, -10.955],
        [3.081, 10.085, -10.981],
        [3.315, 8.47, -13.254],
        [3.439, 7.91, -11.913],
        [1.846, 8.365, -10.156],
        [1.091, 9.158, -9.167],
        [1.355, 8.587, -7.772],
        [1.393, 7.371, -7.574],
        [1.625, 9.513, -6.851],
        [2.045, 9.187, -5.471],
        [1.203, 9.91, -4.411],
        [0.519, 10.884, -4.708],
        [1.149, 9.269, -3.244],
        [0.62, 9.867, -2.011],
        [1.447, 9.408, -0.806],
        [1.613, 8.206, -0.591],
        [2.009, 10.387, -0.089],
        [2.761, 10.16, 1.151],
        [1.832, 10.156, 2.373],
        [1.352, 11.174, 2.853],
        ]

    assert observed_coords == expected_coords

    expected_chain_ranges = {"B": (0, 19)}

    assert observed_chain_ranges == expected_chain_ranges


def test_pdb2fastadic():
    """Test the generation of the fastadic."""
    protein_f = Path(golden_data, "protein.pdb")
    dna_f = Path(golden_data, "dna.pdb")
    ligand_f = Path(golden_data, "ligand.pdb")

    observed_prot_fastadic = pdb2fastadic(protein_f)
    expected_prot_fastadic = {"B": {1: "M", 2: "F", 3: "Q", 4: "Q", 5: "E"}}

    assert observed_prot_fastadic == expected_prot_fastadic

    observed_dna_fastadic = pdb2fastadic(dna_f)
    expected_dna_fastadic = {
        "B": {
            1: "A",
            2: "G",
            3: "T",
            4: "A",
            5: "C",
            28: "A",
            29: "A",
            30: "G",
            31: "T",
            32: "T",
            }
        }

    assert observed_dna_fastadic == expected_dna_fastadic

    observed_ligand_fastadic = pdb2fastadic(ligand_f)
    expected_ligand_fastadic = {"B": {500: "X"}}

    assert observed_ligand_fastadic == expected_ligand_fastadic


def test_get_align():
    """Test the selection of the align function."""
    align_func = get_align(method="sequence")
    assert callable(align_func)

    align_func = get_align(method="structure", **{"lovoalign_exec": ""})
    assert callable(align_func)


# Need dependency to test this
# def test_align_strct():
#     pass


def test_align_seq():
    """Test the sequence alignment."""
    ref = Path(golden_data, "protein.pdb")
    mod = Path(golden_data, "protein_renumb.pdb")

    with tempfile.TemporaryDirectory() as tmpdirname:

        observed_numb_dic = align_seq(ref, mod, tmpdirname)
        expected_numb_dic = {"B": {1: 101, 2: 102, 3: 110, 5: 112}}

        assert observed_numb_dic == expected_numb_dic

        expected_aln_f = Path(tmpdirname, "blosum62_B.aln")

        assert expected_aln_f.exists()

        observed_aln = open(expected_aln_f).readlines()
        expected_aln = [
            f"MFQQE{os.linesep}",
            f"|||-|{os.linesep}",
            f"MFQ-E{os.linesep}",
            ]

        assert observed_aln == expected_aln


def test_make_range():
    """Test the expansion of a chain dic into ranges."""
    chain_range_dic = {"A": [1, 2, 4], "B": [100, 110, 200]}
    observed_range_dic = make_range(chain_range_dic)
    expected_range_dic = {"A": (1, 4), "B": (100, 200)}
    assert observed_range_dic == expected_range_dic


def test_dump_as_izone():
    """Test the generation of .izone file."""
    numb_dic = {"B": {1: 101, 2: 102, 3: 110, 5: 112}}
    with tempfile.NamedTemporaryFile() as fp:

        dump_as_izone(fp.name, numb_dic)

        assert Path(fp.name).stat().st_size != 0

        observed_izone = open(fp.name).readlines()
        expected_izone = [
            f"ZONE B1:B101{os.linesep}",
            f"ZONE B2:B102{os.linesep}",
            f"ZONE B3:B110{os.linesep}",
            f"ZONE B5:B112{os.linesep}",
            ]

        assert observed_izone == expected_izone


# def test_write_coord_dic():
#     pass


# def test_write_coords():
#     pass


# def test_write_pymol_viz():
#     pass
