"""
Library of functions to perform sequence and structural alignments.

Main functions
--------------

* :py:func:`calc_rmsd`
* :py:func:`centroid`
* :py:func:`kabsch`
* :py:func:`load_coords`
* :py:func:`pdb2fastadic`
* :py:func:`get_atoms`
* :py:func:`get_align`
* :py:func:`align_struct`
* :py:func:`align_seq`
* :py:func:`make_range`
* :py:func:`dump_as_izone`
"""
import os
import shlex
import subprocess
from functools import partial
from pathlib import Path

import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq

from haddock import log
from haddock.core.typing import AtomsDict, FilePath, Literal, NDFloat, Optional
from haddock.libs.libio import pdb_path_exists
from haddock.libs.libontology import PDBFile, PDBPath
from haddock.libs.libpdb import split_by_chain


RES_TO_BE_IGNORED = ["SHA", "WAT"]

PROT_RES = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
]

DNA_RES = ["DA", "DC", "DT", "DG"]
# Backbone
PROT_ATOMS = ["C", "N", "CA", "O"]
# Side chains
PROT_SIDE_CHAINS_DICT = {
    "ALA": ["C", "N", "CA", "O", "CB"],
    "ARG": ["C", "N", "CA", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "ASN": ["C", "N", "CA", "O", "CB", "CG", "OD1", "ND2"],
    "ASP": ["C", "N", "CA", "O", "CB", "CG", "OD1", "OD2"],
    "CYS": ["C", "N", "CA", "O", "CB", "SG"],
    "GLN": ["C", "N", "CA", "O", "CB", "CG", "CD", "OE1", "NE2"],
    "GLU": ["C", "N", "CA", "O", "CB", "CG", "CD", "OE1", "OE2"],
    "GLY": ["C", "N", "CA", "O"],
    "HIS": ["C", "N", "CA", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"],
    "ILE": ["C", "N", "CA", "O", "CB", "CG1", "CG2", "CD1"],
    "LEU": ["C", "N", "CA", "O", "CB", "CG", "CD1", "CD2"],
    "LYS": ["C", "N", "CA", "O", "CB", "CG", "CD", "CE", "NZ"],
    "MET": ["C", "N", "CA", "O", "CB", "CG", "SD", "CE"],
    "PHE": ["C", "N", "CA", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PRO": ["C", "N", "CA", "O", "CB", "CG", "CD"],
    "SER": ["C", "N", "CA", "O", "CB", "OG"],
    "THR": ["C", "N", "CA", "O", "CB", "OG1", "CG2"],
    "TRP": [
        "C",
        "N",
        "CA",
        "O",
        "CB",
        "CG",
        "CD1",
        "CD2",
        "NE1",
        "CE2",
        "CE3",
        "CZ2",
        "CZ3",
        "CH2",
    ],  # noqa: E501
    "TYR": [
        "C",
        "N",
        "CA",
        "O",
        "CB",
        "CG",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CZ",
        "OH",
    ],  # noqa: E501
    "VAL": ["C", "N", "CA", "O", "CB", "CG1", "CG2"],
}

# Bases
DNA_ATOMS = [
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
]

RNA_RES = ["A", "G", "C", "U"]
RNA_ATOMS = ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]

DNA_FULL_DICT = {
    "DA": [
        "P",
        "O1P",
        "O2P",
        "O5'",
        "C5'",
        "C4'",
        "O4'",
        "C1'",
        "N9",
        "C4",
        "N3",
        "C2",
        "N1",
        "C6",
        "N6",
        "C5",
        "N7",
        "C8",
        "C2'",
        "C3'",
        "O3'",
    ],  # noqa: E501
    "DG": [
        "P",
        "O1P",
        "O2P",
        "O5'",
        "C5'",
        "C4'",
        "O4'",
        "C1'",
        "N9",
        "C4",
        "N3",
        "C2",
        "N2",
        "N1",
        "C6",
        "O6",
        "C5",
        "N7",
        "C8",
        "C2'",
        "C3'",
        "O3'",
    ],  # noqa: E501
    "DC": [
        "P",
        "O1P",
        "O2P",
        "O5'",
        "C5'",
        "C4'",
        "O4'",
        "C1'",
        "N1",
        "C6",
        "C2",
        "O2",
        "N3",
        "C4",
        "N4",
        "C5",
        "C2'",
        "C3'",
        "O3'",
    ],
    "DT": [
        "P",
        "O1P",
        "O2P",
        "O5'",
        "C5'",
        "C4'",
        "O4'",
        "C1'",
        "N1",
        "C6",
        "C2",
        "O2",
        "N3",
        "C4",
        "O4",
        "C5",
        "C7",
        "C2'",
        "C3'",
        "O3'",
    ],
}

RNA_FULL_DICT = {
    "A": [
        "P",
        "OP1",
        "OP2",
        "O5'",
        "C5'",
        "C4'",
        "O4'",
        "C3'",
        "O3'",
        "C2'",
        "O2'",
        "C1'",
        "N9",
        "C8",
        "N7",
        "C5",
        "C6",
        "N6",
        "N1",
        "C2",
        "N3",
        "C4",
    ],  # noqa: E501
    "G": [
        "P",
        "OP1",
        "OP2",
        "O5'",
        "C5'",
        "C4'",
        "O4'",
        "C3'",
        "O3'",
        "C2'",
        "O2'",
        "C1'",
        "N9",
        "C8",
        "N7",
        "C5",
        "C6",
        "O6",
        "N1",
        "C2",
        "N2",
        "N3",
        "C4",
    ],  # noqa: E501
    "C": [
        "P",
        "OP1",
        "OP2",
        "O5'",
        "C5'",
        "C4'",
        "O4'",
        "C3'",
        "O3'",
        "C2'",
        "O2'",
        "C1'",
        "N9",
        "C5",
        "C6",
        "O6",
        "N1",
        "C2",
        "N2",
        "N3",
        "C4",
        "N4",
    ],  # noqa: E501
    "U": [
        "P",
        "OP1",
        "OP2",
        "O5'",
        "C5'",
        "C4'",
        "O4'",
        "C3'",
        "O3'",
        "C2'",
        "O2'",
        "C1'",
        "N1",
        "C2",
        "O2",
        "N3",
        "C4",
        "O4",
        "C5",
        "C6",
    ],
}


class ALIGNError(Exception):
    """Raised when something goes wrong with the ALIGNMENT library."""

    def __init__(self, msg: object = "") -> None:
        self.msg = msg
        super().__init__(self.msg)


def calc_rmsd(V: NDFloat, W: NDFloat) -> float:
    """
    Calculate the RMSD from two vectors.

    Parameters
    ----------
    V : np.array dtype=float, shape=(n_atoms,3)
    W : np.array dtype=float, shape=(n_atoms,3)

    Returns
    -------
    rmsd : float
    """
    diff = np.array(V) - np.array(W)
    N = len(V)
    rmsd = np.sqrt((diff * diff).sum() / N)
    return rmsd


def kabsch(P: NDFloat, Q: NDFloat) -> NDFloat:
    """
    Find the rotation matrix using Kabsch algorithm.

    Parameters
    ----------
    P : np.array dtype=float, shape=(n_atoms,3)
    Q : np.array dtype=float, shape=(n_atoms,3)

    Returns
    -------
    U : np.array dtype=float, shape=(3,3)
    """
    # Covariance matrix
    P = np.array(P)
    Q = np.array(Q)
    C = np.dot(np.transpose(P), Q)
    # use SVD
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    # Create Rotation matrix U
    U = np.dot(V, W)
    return U


def centroid(X: NDFloat) -> NDFloat:
    """
    Get the centroid.

    Parameters
    ----------
    X : np.array dtype=float, shape=(n_atoms,3)

    Returns
    -------
    C : np.array dtype=float, shape=(3,)
    """
    X = np.array(X)
    C = X.mean(axis=0)
    return C


def load_coords(
    pdb_f, atoms, filter_resdic=None, numbering_dic=None, model2ref_chain_dict=None
):
    """
    Load coordinates from PDB.

    Parameters
    ----------
    pdb_f : PDBFile

    atoms : dict
        dictionary of atoms

    filter_resdic : dict
        dictionary of residues to be loaded (one list per chain)

    numbering_dic : dict
        dict of numbering dictionaries (one dictionary per chain)

    Returns
    -------
    coord_dic : dict
        dictionary of coordinates (one per chain)

    chain_ranges: dict
        dictionary of chain ranges
    """
    coord_dic: CoordsDict = {}
    chain_dic: ResDict = {}
    idx = 0
    if isinstance(pdb_f, PDBFile):
        pdb_f = pdb_f.rel_path
    with open(pdb_f, "r") as fh:
        for line in fh.readlines():
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                resname = line[17:20].strip()
                if resname in RES_TO_BE_IGNORED:
                    continue
                if model2ref_chain_dict:
                    chain = model2ref_chain_dict[line[21]]
                else:
                    chain = line[21]
                resnum = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords = np.asarray([x, y, z])
                if numbering_dic and model2ref_chain_dict:
                    try:
                        resnum = numbering_dic[chain][resnum]
                    except KeyError:
                        # this residue is not matched, and so it should
                        #  not be considered
                        # self.log(
                        #     f"WARNING: {chain}.{resnum}.{atom_name}"
                        #     " was not matched!"
                        #     )
                        continue
                # identifier = f"{chain}.{resnum}.{atom_name}"
                identifier = (chain, resnum, atom_name)
                if atom_name not in atoms[resname]:
                    continue
                if chain not in chain_dic:
                    chain_dic[chain] = []
                if filter_resdic:
                    # Only retrieve coordinates from the filter_resdic
                    if chain in filter_resdic and resnum in filter_resdic[chain]:
                        coord_dic[identifier] = coords
                        chain_dic[chain].append(idx)
                        idx += 1
                else:
                    # retrieve everything
                    coord_dic[identifier] = coords
                    chain_dic[chain].append(idx)
                    idx += 1
    chain_ranges: ChainsRange = {}
    for chain, indice in chain_dic.items():
        if not indice:
            # this may happen when filter_resdic is defined on a different set
            # of chains
            raise ALIGNError(f"Chain matching error on {pdb_f}, chain {chain}")
        else:
            min_idx = min(indice)
            max_idx = max(indice)
            chain_ranges[chain] = (min_idx, max_idx)
    return coord_dic, chain_ranges


def get_atoms(pdb: PDBPath, full: bool = False) -> AtomsDict:
    """
    Identify what is the molecule type of each PDB.

    Parameters
    ----------
    pdb : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`
        PDB file to have its atoms identified

    Returns
    -------
    atom_dic : dict
        dictionary of atoms
    """
    atom_dic: AtomsDict = {}
    atom_dic.update((r, PROT_ATOMS) for r in PROT_RES)
    atom_dic.update((r, DNA_ATOMS) for r in DNA_RES)
    atom_dic.update((r, RNA_ATOMS) for r in RNA_RES)
    if full:
        atom_dic.update(PROT_SIDE_CHAINS_DICT)
        atom_dic.update(DNA_FULL_DICT)
        atom_dic.update(RNA_FULL_DICT)

    if isinstance(pdb, PDBFile):
        pdb = pdb.rel_path

    exists, msg = pdb_path_exists(pdb)
    if not exists:
        raise Exception(msg)

    with open(pdb) as fh:
        for line in fh.readlines():
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()
                atom_name = line[12:16].strip()
                element = line[76:78].strip()
                if (
                    resname not in PROT_RES
                    and resname not in DNA_RES
                    and resname not in RNA_RES
                    and resname not in RES_TO_BE_IGNORED
                ):
                    # its neither DNA/RNA nor protein, use the heavy atoms
                    # WARNING: Atoms that belong to unknown residues must
                    #  be bound to a residue name;
                    #   For example: residue NEP, also contains
                    #  CB and CG atoms, if we do not bind it to the
                    #  residue name, the next functions will include
                    #  CG and CG atoms in the calculations for all
                    #  other residue names
                    if element != "H":
                        if resname not in atom_dic:
                            atom_dic[resname] = []
                        if atom_name not in atom_dic[resname]:
                            atom_dic[resname].append(atom_name)
    return atom_dic


ResCode = Literal[
    "C",
    "D",
    "S",
    "Q",
    "K",
    "I",
    "P",
    "T",
    "F",
    "N",
    "G",
    "H",
    "L",
    "R",
    "W",
    "A",
    "V",
    "E",
    "Y",
    "M",
    "A",
    "G",
    "C",
    "T",
    "X",
]
"""
The single letter code of a residue.

Unrecognized residues' code is `X`.
"""


def pdb2fastadic(pdb_f: PDBPath) -> dict[str, dict[int, str]]:
    """
    Write the sequence as a fasta.

    Parameters
    ----------
    pdb_f : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`

    Returns
    -------
    seq_dic : dict
        dict of fasta sequences (one per chain)
    """
    res_codes = dict(
        [
            ("CYS", "C"),
            ("ASP", "D"),
            ("SER", "S"),
            ("GLN", "Q"),
            ("LYS", "K"),
            ("ILE", "I"),
            ("PRO", "P"),
            ("THR", "T"),
            ("PHE", "F"),
            ("ASN", "N"),
            ("GLY", "G"),
            ("HIS", "H"),
            ("LEU", "L"),
            ("ARG", "R"),
            ("TRP", "W"),
            ("ALA", "A"),
            ("VAL", "V"),
            ("GLU", "E"),
            ("TYR", "Y"),
            ("MET", "M"),
            ("DA", "A"),
            ("DG", "G"),
            ("DC", "C"),
            (
                "DT",
                "T",
            ),  # 9/8/2023: adding non-standard amino-acids (src/haddock/cns/toppar/protein-allhdg5-4.top) # noqa: E501
            ("ALY", "K"),
            ("ASH", "D"),
            ("CFE", "C"),
            ("CSP", "C"),
            ("CYC", "C"),
            ("CYF", "C"),
            ("CYM", "C"),
            ("DDZ", "A"),
            ("GLH", "E"),
            ("HLY", "P"),
            ("HY3", "P"),
            ("HYP", "P"),
            ("M3L", "K"),
            ("MLY", "K"),
            ("MLZ", "K"),
            ("MSE", "M"),
            ("NEP", "H"),
            ("PNS", "S"),
            ("PTR", "Y"),
            ("SEP", "S"),
            ("TOP", "T"),
            ("TYP", "Y"),
            ("TYS", "Y"),
        ]
    )

    seq_dic: dict[str, dict[int, str]] = {}

    if isinstance(pdb_f, PDBFile):
        pdb_f = pdb_f.rel_path

    with open(pdb_f) as fh:
        for line in fh.readlines():
            if line.startswith("ATOM"):
                res_num = int(line[22:26])
                res_name = line[17:20].strip()
                chain = line[21]
                if res_name in RES_TO_BE_IGNORED:
                    continue
                try:
                    one_letter = res_codes[res_name]
                except KeyError:
                    one_letter = "X"
                if chain not in seq_dic:
                    seq_dic[chain] = {}
                seq_dic[chain][res_num] = one_letter
    return seq_dic


def get_align(
    method: str, lovoalign_exec: FilePath
) -> partial[dict[str, dict[int, int]]]:
    """
    Get the alignment function.

    Parameters
    ----------
    method : str
        Available options: ``sequence`` and ``structure``.

    lovoalign_exec : str
        Path to the lovoalign executable.

    Returns
    -------
    align_func : functools.partial
        desired alignment function
    """
    log.debug(f"Using {method} alignment")
    if method == "structure":
        align_func = partial(align_strct, lovoalign_exec=lovoalign_exec)
    elif method == "sequence":
        align_func = partial(align_seq)
    else:
        available_alns = ("sequence", "structure")
        raise ValueError(
            f"Alignment method {method!r} not recognized. "
            f"Available options are {', '.join(available_alns)}"
        )
    return align_func


def align_strct(
    reference: PDBFile,
    model: PDBFile,
    output_path: FilePath,
    lovoalign_exec: Optional[FilePath] = None,
) -> dict[str, dict[int, int]]:
    """
    Structuraly align and get numbering relationship.

    Parameters
    ----------
    reference : :py:class:`haddock.libs.libontology.PDBFile`

    model : :py:class:`haddock.libs.libontology.PDBFile`

    output_path : Path

    lovoalign_exec : Path
        lovoalign executable

    Returns
    -------
    numbering_dic : dict
        dict of numbering dictionaries (one dictionary per chain)
    """
    if lovoalign_exec is None:
        log.error(
            "Structural alignment needs LovoAlign " "get it at github.com/m3g/lovoalign"
        )
        raise ALIGNError("Path to LovoAlign executable required.")

    if not lovoalign_exec:
        raise ALIGNError("lovoalign_exec parameter not defined ")

    if not os.access(lovoalign_exec, os.X_OK):
        raise ALIGNError(f"{lovoalign_exec!r} for LovoAlign is not executable")

    numbering_dic: dict[str, dict[int, int]] = {}
    protein_a_dic = {
        e.stem.split("_")[-1]: e for e in split_by_chain(reference.rel_path)
    }
    protein_b_dic = {e.stem.split("_")[-1]: e for e in split_by_chain(model.rel_path)}

    # check if chain ids match
    if protein_a_dic.keys() != protein_b_dic.keys():
        # TODO: Make this a clearer raise
        return numbering_dic

    for chain in protein_a_dic.keys():
        pa_seqdic = pdb2fastadic(protein_a_dic[chain])
        pb_seqdic = pdb2fastadic(protein_b_dic[chain])
        # logging.debug(f"Structurally aligning chain {chain}")
        numbering_dic[chain] = {}
        cmd = (
            f"{lovoalign_exec} -p1 {protein_a_dic[chain]} "
            f"-p2 {protein_b_dic[chain]} "
            f"-c1 {chain} -c2 {chain}"
        )

        # logging.debug(f"Command is: {cmd}")
        p = subprocess.run(shlex.split(cmd), capture_output=True, text=True)
        lovoalign_out = p.stdout.split(os.linesep)

        # we don"t need the splitted proteins anymore
        protein_a_dic[chain].unlink()
        protein_b_dic[chain].unlink()

        # find out where the alignment starts and ends
        alignment_pass = True
        for i, line in enumerate(lovoalign_out):
            if "SEQUENCE ALIGNMENT" in line:
                # there are 2 extra white lines after this header
                alignment_start_index = i + 2
            elif "FINAL" in line:
                # there are 2 extra white lines after this header
                alignment_end_index = i - 2
            elif "ERROR" in line:
                failed_pdb = line.split()[-1]
                _msg = f"LovoAlign could not read {failed_pdb} " "is it a ligand?"
                log.warning(_msg)
                alignment_pass = False

                for elem in [k for k in pa_seqdic[chain]]:
                    numbering_dic[chain][elem] = elem

        if not alignment_pass:
            # This alignment failed, move on to the next
            log.warning(
                f"Skipping alignment of chain {chain}, " "used sequential matching"
            )
            continue

        aln_l = lovoalign_out[alignment_start_index:alignment_end_index]

        # dump this alignment to a file
        aln_fname = Path(output_path, f"lovoalign_{chain}.aln")
        log.debug(f"Writing alignment to {aln_fname.name}")
        with open(aln_fname, "w") as fh:
            fh.write(os.linesep.join(aln_l))

        # remove the line between the alignment segments
        alignment = [aln_l[i : i + 3][:2] for i in range(0, len(aln_l), 3)]
        # 100% (5 identical nucleotides / min(length(A),length(B))).
        len_seq_a = len(pa_seqdic[chain])
        len_seq_b = len(pb_seqdic[chain])
        identity = (
            (len_seq_a - sum([e[0].count("-") for e in alignment]))
            / min(len_seq_a, len_seq_b)
            * 100
        )

        if identity <= 40.0:
            log.warning(
                f'"Structural" identity of chain {chain} is {identity:.2f}%,'
                " please check the results carefully"
            )
        else:
            log.info(f'"Structural" identity of chain {chain} is {identity:.2f}%')

        # logging.debug("Reading alignment and matching numbering")
        for element in alignment:
            line_a, line_b = element

            resnum_a, seq_a, _ = line_a.split()
            resnum_b, seq_b, _ = line_b.split()

            resnum_a = int(resnum_a) - 1  # type: ignore
            resnum_b = int(resnum_b) - 1  # type: ignore

            for resname_a, resname_b in zip(seq_a, seq_b):
                if resname_a != "-":
                    resnum_a += 1  # type: ignore

                if resname_b != "-":
                    resnum_b += 1  # type: ignore

                if resname_a != "-" and resname_b != "-":
                    numbering_dic[chain][resnum_b] = resnum_a  # type: ignore

    izone_fname = Path(output_path, "lovoalign.izone")
    log.debug(f"Saving .izone to {izone_fname.name}")
    dump_as_izone(izone_fname, numbering_dic)

    return numbering_dic


def write_alignment(top_aln, output_path, ref_ch):
    """
    Write the alignment to a file.

    Parameters
    ----------
    top_aln : Bio.Align.PairwiseAlignments
        alignment object

    ref_ch : str
        reference chain
    """
    aln_fname = Path(output_path, f"blosum62_{ref_ch}.aln")
    log.debug(f"Writing alignment to {aln_fname.name}")
    with open(aln_fname, "w") as fh:
        fh.write(str(top_aln))
    return aln_fname


def sequence_alignment(seq_ref, seq_model):
    """
    Perform a sequence alignment.

    Parameters
    ----------
    seq_ref : str
        reference sequence

    seq_model : str
        model sequence

    Returns
    -------
    identity : float
        sequence identity

    top_aln : Bio.Align.PairwiseAlignments
        alignment object

    aln_ref_seg : tuple
        aligned reference segment

    aln_mod_seg : tuple
        aligned model segment
    """
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    alns = aligner.align(seq_ref, seq_model)
    top_aln = alns[0]

    aln_denom = min(len(seq_ref), len(seq_model))
    identity = (str(top_aln).count("|") / float(aln_denom)) * 100

    aln_ref_seg, aln_mod_seg = top_aln.aligned
    # check alignment length
    len_ref = len(aln_ref_seg)
    len_mod = len(aln_mod_seg)
    if len_ref != len_mod:
        raise ALIGNError(f"Alignment length mismatch ({len_ref} != {len_mod})")

    return identity, top_aln, aln_ref_seg, aln_mod_seg


class SeqAlign:
    """SeqAlign class."""

    def __init__(
        self,
    ):
        """Initialize the class."""
        self.align_dic = {}  # the alignment dictionary, important for CAPRI
        self.model2ref_chain_dict = {}  # model to reference chain dictionary
        self.ref2model_chain_dict = {}  # reference to model chain dictionary
        self.seqdic_ref = None  # reference sequence dictionary
        self.seqdic_model = None  # model sequence dictionary
        self.seqs_ref = {}  # reference sequences
        self.seqs_model = {}  # model sequences
        self.identities = []  # list of identity values
        self.aln_model_segs = []  # list of aligned model segments
        self.aln_ref_segs = []  # list of aligned reference segments
        self.top_alns = []  # list of alignment objects

    def postprocess_alignment(self, ref_ch, mod_ch, align_id):
        """
        Postprocess the alignment.

        Parameters
        ----------
        ref_ch : str
            reference chain

        mod_ch : str
            model chain

        align_id : int
            alignment id (index of the alignment)
        """
        self.align_dic[ref_ch] = {}
        if not any(e for e in self.top_alns[align_id].aligned):
            # No alignment!
            log.warning(
                f"No alignment for chain {ref_ch} is it protein/dna-rna? "
                "Matching sequentially"
            )
            if all("X" in s for s in self.seqs_ref[ref_ch]) and all(
                "X" in s for s in self.seqs_model[mod_ch]
            ):
                # this sequence contains only ligands, do it manually
                if len(self.seqs_ref[ref_ch]) != len(self.seqs_model[mod_ch]):
                    # we cannot handle this
                    # FIXME: This should raise a proper exception instead
                    raise f"Cannot align chain {mod_ch}"  # noqa: B016
                for ref_res, model_res in zip(
                    self.seqdic_ref[ref_ch], self.seqdic_model[mod_ch]
                ):
                    self.align_dic[ref_ch].update({model_res: ref_res})
        else:
            identity = self.identities[align_id]
            if identity <= 40.0:
                # Identity is very low
                log.warning(
                    f"Sequence identity of chain {ref_ch} is "
                    f"{identity:.2f}%, please check the results carefully"
                )
            else:
                log.debug(
                    f"Sequence identity between chain {ref_ch} "
                    f" of reference and {mod_ch} of model is "
                    f"{identity:.2f}%"
                )
            for ref_segment, model_segment in zip(
                self.aln_ref_segs[align_id], self.aln_model_segs[align_id]
            ):
                start_ref_segment, end_ref_segment = ref_segment
                start_model_segment, end_model_segment = model_segment
                reslist_ref = list(self.seqdic_ref[ref_ch].keys())[
                    start_ref_segment:end_ref_segment
                ]
                reslist_model = list(self.seqdic_model[mod_ch].keys())[
                    start_model_segment:end_model_segment
                ]
                for _ref_res, _model_res in zip(reslist_ref, reslist_model):
                    self.align_dic[ref_ch].update({_model_res: _ref_res})


def align_seq(reference, model, output_path):
    """
    Sequence align and get the numbering relationship.

    Parameters
    ----------
    reference : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`

    model : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`

    output_path : Path

    Returns
    -------
    align_dic : dict
        dictionary of sequence alignments (one per chain)
    """
    # print(f"running align_seq on {reference} and {model}")
    SeqAln = SeqAlign()
    SeqAln.seqdic_ref = pdb2fastadic(reference)
    SeqAln.seqdic_model = pdb2fastadic(model)

    # assign sequences
    for ref_ch in SeqAln.seqdic_ref.keys():
        ref_seq = Seq("".join(SeqAln.seqdic_ref[ref_ch].values()))
        SeqAln.seqs_ref[ref_ch] = ref_seq
    for mod_ch in SeqAln.seqdic_model.keys():
        mod_seq = Seq("".join(SeqAln.seqdic_model[mod_ch].values()))
        SeqAln.seqs_model[mod_ch] = mod_seq

    # check if chain ids match
    if SeqAln.seqdic_ref.keys() != SeqAln.seqdic_model.keys():
        # they do not match, we need to do chain matching
        n_partners = min(len(SeqAln.seqdic_ref.keys()), len(SeqAln.seqdic_model.keys()))

        # unique combinations of chains
        combs = []
        for ref_key in SeqAln.seqdic_ref.keys():
            for mod_key in SeqAln.seqdic_model.keys():
                combs.append((ref_key, mod_key))
        identities = []
        aln_model_segs = []
        aln_ref_segs = []
        top_alns = []
        # loop over combinations
        for ref_ch, mod_ch in combs:
            # align
            identity, top_aln, aln_ref_seg, aln_mod_seg = sequence_alignment(
                SeqAln.seqs_ref[ref_ch], SeqAln.seqs_model[mod_ch]
            )
            # append values to lists
            identities.append(identity)
            aln_model_segs.append(aln_mod_seg)
            aln_ref_segs.append(aln_ref_seg)
            top_alns.append(top_aln)
        # chain matching
        matches = 0  # counter for matched chains
        while matches < n_partners:
            # get the best alignment
            max_identity = max(identities)
            max_idx = identities.index(max_identity)
            # assigning chains
            ref_ch, mod_ch = combs[max_idx]
            SeqAln.model2ref_chain_dict[mod_ch] = ref_ch
            SeqAln.ref2model_chain_dict[ref_ch] = mod_ch
            SeqAln.aln_ref_segs.append(aln_ref_segs[max_idx])
            SeqAln.aln_model_segs.append(aln_model_segs[max_idx])
            SeqAln.identities.append(identities[max_idx])
            SeqAln.top_alns.append(top_alns[max_idx])

            # writing the alignment
            write_alignment(top_alns[max_idx], output_path, ref_ch)

            # postprocess alignment
            SeqAln.postprocess_alignment(ref_ch, mod_ch, matches)
            # update identities to avoid double matches
            identities = [
                identities[n] if combs[n][0] != ref_ch and combs[n][1] != mod_ch else -1
                for n in range(len(combs))
            ]
            matches += 1
        log.info(f"model2ref chain matching is {SeqAln.model2ref_chain_dict}")
    else:
        # chains do match. no need to do chain matching
        matches = 0
        for ref_ch, mod_ch in zip(SeqAln.seqdic_ref, SeqAln.seqdic_model):
            SeqAln.model2ref_chain_dict[mod_ch] = ref_ch
            SeqAln.ref2model_chain_dict[ref_ch] = mod_ch
            if ref_ch != mod_ch:
                raise ALIGNError(f"Chain mismatch: {ref_ch} != {mod_ch}")
            # extract sequences
            seq_ref = SeqAln.seqs_ref[ref_ch]
            seq_model = SeqAln.seqs_model[mod_ch]
            # sequence alignment
            identity, top_aln, aln_ref_seg, aln_mod_seg = sequence_alignment(
                seq_ref, seq_model
            )
            # update quantities
            SeqAln.identities.append(identity)
            SeqAln.aln_model_segs.append(aln_mod_seg)
            SeqAln.aln_ref_segs.append(aln_ref_seg)
            SeqAln.top_alns.append(top_aln)
            # write alignment
            write_alignment(top_aln, output_path, ref_ch)
            # postprocess alignment
            SeqAln.postprocess_alignment(ref_ch, mod_ch, matches)
            matches += 1
    # dump the .izone file
    izone_fname = Path(output_path, "blosum62.izone")
    log.debug(f"Saving .izone to {izone_fname.name}")
    dump_as_izone(izone_fname, SeqAln.align_dic, SeqAln.ref2model_chain_dict)

    return SeqAln.align_dic, SeqAln.model2ref_chain_dict


def make_range(chain_range_dic: dict[str, list[int]]) -> dict[str, tuple[int, int]]:
    """
    Expand a chain dictionary into ranges.

    Parameters
    ----------
    chain_range_dic : dict
        dictionary of chain indexes (one list per chain)

    Returns
    -------
    chain_ranges : dict
        dictionary of chain ranges (one tuple per chain)
    """
    chain_ranges: dict[str, tuple[int, int]] = {}
    for chain in chain_range_dic:
        min_idx = min(chain_range_dic[chain])
        max_idx = max(chain_range_dic[chain])
        chain_ranges[chain] = (min_idx, max_idx)
    return chain_ranges


def dump_as_izone(fname, numbering_dic, model2ref_chain_dict=None):
    """
    Dump the numbering dictionary as .izone.

    Parameters
    ----------
    fname : str
        output filename

    numbering_dic : dict
        dict of numbering dictionaries (one dictionary per chain)
    """
    # FIXME: Collapse the izones so its faster to load in profit
    with open(fname, "w") as fh:
        for chain in numbering_dic:
            for bound_res in numbering_dic[chain]:
                unbound_res = numbering_dic[chain][bound_res]
                unb_chain = chain
                if model2ref_chain_dict:
                    unb_chain = model2ref_chain_dict[chain]
                #
                izone_str = (
                    "ZONE "
                    f"{chain}{bound_res}:{unb_chain}{unbound_res}"
                    f"{os.linesep}"
                )
                fh.write(izone_str)
