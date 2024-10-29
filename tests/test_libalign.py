"""Test the libalign library."""

import os
import tempfile
from pathlib import Path

import numpy as np
import pytest

from haddock.libs.libalign import (
    ALIGNError,
    align_seq,
    calc_rmsd,
    centroid,
    check_chains,
    check_common_atoms,
    dump_as_izone,
    get_align,
    get_atoms,
    kabsch,
    load_coords,
    make_range,
    pdb2fastadic,
    rearrange_xyz_files,
    )

from . import golden_data


def array_to_list(np_array):
    """Transform a numpy array in a nested list."""
    return [list(e) for e in np_array]


def test_kabsch():
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

    observed_U = kabsch(P, Q)

    expected_U = [
        [0.38845040189428726, 0.38160307568742435, -0.8387403518932808],
        [-0.20700756088693778, 0.9230933741949091, 0.32410876608493205],
        [0.8979165634978604, 0.047725414019726436, 0.43757071411697324],
    ]

    np.testing.assert_allclose(np.asarray(expected_U), observed_U)


def test_calc_rmsd():
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

    rmsd = calc_rmsd(V, W)

    assert round(rmsd, 2) == 12.02


def test_centroid():
    """Test the centroid calculation."""
    X = [
        [6.0952658227848495, -5.630326582278489, 8.033610126582303],
        [0.5082658227848498, 1.733673417721512, 10.833610126582304],
        [7.42726582278485, 7.67067341772151, -5.757389873417694],
    ]

    observed_centroid = centroid(X)
    observed_centroid = list(observed_centroid)
    expected_centroid = [
        4.6769324894515165,
        1.2580067510548443,
        4.369943459915638,
    ]

    assert observed_centroid == expected_centroid


def test_load_coords():
    """Test the loading of coordinates."""
    # pdb_f = protprot_input_list[0]
    pdb_f = Path(golden_data, "protein.pdb")
    atoms = get_atoms(pdb_f)
    (
        observed_coord_dic,
        observed_chain_ranges,
    ) = load_coords(pdb_f, atoms)
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


def test_wrong_filtered_resid_error_load_coords():
    """Test the residue matching error with an uncompatible resdic."""
    filter_resdic_wrongres = {"B": [7, 8, 9]}  # protein has only residues 1-5
    pdb_f = Path(golden_data, "protein.pdb")
    atoms = get_atoms(pdb_f)
    with pytest.raises(ALIGNError):
        load_coords(pdb_f, atoms, filter_resdic=filter_resdic_wrongres)


def test_wrong_filtered_chain_error_load_coords():
    """Test the chain matching error with an uncompatible resdic."""
    filter_resdic_wrongchain = {"A": [1, 2, 3]}  # protein has only chain B
    pdb_f = Path(golden_data, "protein.pdb")
    atoms = get_atoms(pdb_f)
    with pytest.raises(ALIGNError):
        load_coords(pdb_f, atoms, filter_resdic=filter_resdic_wrongchain)


def test_get_atoms():
    """Test the identification of atoms."""
    pdb_list = [
        Path(golden_data, "protein.pdb"),
        Path(golden_data, "dna.pdb"),
        Path(golden_data, "ligand.pdb"),
    ]
    observed_atom_dic = {}
    for p in pdb_list:
        observed_atom_dic.update(get_atoms(p))
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
        "A": ["P", "O5'", "C5'", "C4'", "C3'", "O3'"],
        "G": ["P", "O5'", "C5'", "C4'", "C3'", "O3'"],
        "C": ["P", "O5'", "C5'", "C4'", "C3'", "O3'"],
        "U": ["P", "O5'", "C5'", "C4'", "C3'", "O3'"],
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
    align_func = get_align(method="sequence", lovoalign_exec="")
    assert callable(align_func)

    align_func = get_align(method="structure", lovoalign_exec="")
    assert callable(align_func)


# Need dependency to test this
# def test_align_strct():
#     pass


def test_align_seq():
    """Test the sequence alignment."""
    ref = Path(golden_data, "protein.pdb")
    mod = Path(golden_data, "protein_renumb.pdb")

    with tempfile.TemporaryDirectory() as tmpdirname:

        observed_numb_dic, observed_chm_dict = align_seq(ref, mod, tmpdirname)
        expected_numb_dic = {"B": {101: 1, 102: 2, 110: 3, 112: 5}}
        expected_chm_dict = {"B": "B"}

        assert observed_numb_dic == expected_numb_dic
        assert observed_chm_dict == expected_chm_dict

        expected_aln_f = Path(tmpdirname, "blosum62_B.aln")

        assert expected_aln_f.exists()

        observed_aln = open(expected_aln_f).readlines()
        expected_aln = [
            f"target            0 MFQQE 5{os.linesep}",
            f"                  0 |||-| 5{os.linesep}",
            f"query             0 MFQ-E 4{os.linesep}",
        ]

        assert observed_aln == expected_aln


def test_align_seq_chm():
    """Test the sequence alignment with chain matching."""
    ref = Path(golden_data, "protein.pdb")
    mod = Path(golden_data, "protein_segid.pdb")

    with tempfile.TemporaryDirectory() as tmpdirname:

        observed_numb_dic, observed_chm_dict = align_seq(ref, mod, tmpdirname)
        expected_numb_dic = {"B": {1: 1, 2: 2, 3: 3, 4: 4, 5: 5}}
        expected_chm_dict = {"X": "B"}

        assert observed_numb_dic == expected_numb_dic
        assert observed_chm_dict == expected_chm_dict


def test_align_seq_inverted():
    """Test the sequence alignment with inverted chain."""
    ref = Path(golden_data, "protprot_complex_1.pdb")
    mod = Path(golden_data, "protprot_complex_2_inverted.pdb")

    with tempfile.TemporaryDirectory() as tmpdirname:

        observed_numb_dic, observed_chm_dict = align_seq(ref, mod, tmpdirname)
        print(f"observed_numb_dic: {observed_numb_dic}")
        print(f"observed_chm_dict: {observed_chm_dict}")
        expected_numb_keys = ["A", "B"]
        expected_chm_dict = {"A": "A", "B": "B"}

        assert list(observed_numb_dic.keys()) == expected_numb_keys
        assert observed_chm_dict == expected_chm_dict


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

    chm_ref2model_dict = {"B": "X"}
    with tempfile.NamedTemporaryFile() as fp:

        dump_as_izone(fp.name, numb_dic, chm_ref2model_dict)

        assert Path(fp.name).stat().st_size != 0

        observed_izone = open(fp.name).readlines()
        expected_izone = [
            f"ZONE B1:X101{os.linesep}",
            f"ZONE B2:X102{os.linesep}",
            f"ZONE B3:X110{os.linesep}",
            f"ZONE B5:X112{os.linesep}",
        ]

        assert observed_izone == expected_izone


def test_check_common_atoms():
    """Test the identification of common atoms."""
    ref = Path(golden_data, "protprot_complex_1.pdb")
    mod = Path(golden_data, "protprot_complex_2.pdb")
    models = [ref, mod]

    n_atoms, obs_common_keys = check_common_atoms(models, None, False, 90.0)
    assert n_atoms == 950
    assert len(obs_common_keys) == 950
    assert ("B", 74, "N") in obs_common_keys

    models.append(Path(golden_data, "protein.pdb"))
    with pytest.raises(ALIGNError):
        n_atoms, obs_common_keys = check_common_atoms(models, None, False, 90.0)


def test_rearrange_xyz_files():
    """Test the rearrange_xyz_files function."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        ncores = 4
        # Create a temporary directory with some files
        for i in range(ncores):
            with open(Path(tmpdirname, f"file_{i}.xyz"), "w") as f:
                f.write(f"{i} 0 0 0\n")

        # Test the function
        rearrange_xyz_files("file.xyz", path=tmpdirname, ncores=ncores)

        # Check the files have been renamed
        assert not Path(tmpdirname, "file_0.xyz").exists()
        assert Path(tmpdirname, "file.xyz").exists()
        # Check the content of the file
        with open(Path(tmpdirname, "file.xyz"), "r") as f:
            obs_content = f.read()
        exp_content = os.linesep.join([f"{i} 0 0 0" for i in range(ncores)])
        exp_content += os.linesep
        assert obs_content == exp_content


def test_check_chains():
    """Test correct checking of chains."""
    obs_ch = [
        ["A", "C"],
        ["A", "B"],
        ["S", "E", "B", "A"],
        ["S", "E", "P", "A"],
        ["C", "D"],
    ]

    inp_receptor_chains = ["A", "A", "A", "A", "C"]
    inp_ligand_chains = [
        [],
        [],
        ["B", "E"],
        ["B"],
        ["B"],
    ]

    # assuming exp chains are A and B
    exp_ch = [
        ["A", ["C"]],  # C becomes the ligand
        ["A", ["B"]],  # C becomes the ligand
        ["A", ["B", "E"]],  # S is ignored (B,E are present)
        ["A", ["S", "E", "P"]],  # B is not there, S-E-P become the ligands
        ["C", ["D"]],
    ]  # B is not there, D becomes the ligand

    for n in range(len(obs_ch)):
        obs_r_chain, obs_l_chain = check_chains(
            obs_ch[n], inp_receptor_chains[n], inp_ligand_chains[n]
        )
        exp_r_chain, exp_l_chain = exp_ch[n][0], exp_ch[n][1]
        assert obs_r_chain == exp_r_chain
        assert obs_l_chain == exp_l_chain

