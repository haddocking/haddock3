
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
from haddock.libs.libontology import PDBFile
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


class ALIGNError(Exception):
    """Raised when something goes wrong with the ALIGNMENT library."""

    def __init__(self, msg=""):
        self.msg = msg
        super().__init__(self.msg)


def calc_rmsd(V, W):
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


def kabsch(P, Q):
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


def centroid(X):
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


def load_coords(pdb_f, atoms, filter_resdic=None, numbering_dic=None):
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
    coord_dic = {}
    chain_dic = {}
    idx = 0
    if isinstance(pdb_f, PDBFile):
        pdb_f = pdb_f.rel_path
    with open(pdb_f, "r") as fh:
        for line in fh.readlines():
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                resname = line[17:20].strip()
                chain = line[21]
                resnum = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords = np.asarray([x, y, z])
                if numbering_dic:
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
                    if (
                            chain in filter_resdic
                            and resnum in filter_resdic[chain]
                            ):
                        coord_dic[identifier] = coords
                        chain_dic[chain].append(idx)
                        idx += 1
                else:
                    # retrieve everything
                    coord_dic[identifier] = coords
                    chain_dic[chain].append(idx)
                    idx += 1
    chain_ranges = {}
    for chain in chain_dic:
        if not chain_dic[chain]:
            # this may happen when filter_resdic is defined on a different set
            # of chains
            raise ALIGNError(f"Chain matching error on {pdb_f}, chain {chain}")
        else:
            min_idx = min(chain_dic[chain])
            max_idx = max(chain_dic[chain])
            chain_ranges[chain] = (min_idx, max_idx)
    return coord_dic, chain_ranges


def get_atoms(pdb):
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
    atom_dic = {}
    atom_dic.update(dict((r, PROT_ATOMS) for r in PROT_RES))
    atom_dic.update(dict((r, DNA_ATOMS) for r in DNA_RES))

    if isinstance(pdb, PDBFile):
        pdb = pdb.rel_path

    with open(pdb) as fh:
        for line in fh.readlines():
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()
                atom_name = line[12:16].strip()
                element = line[76:78].strip()
                if (
                        resname not in PROT_RES
                        and resname not in DNA_RES
                        and resname not in RES_TO_BE_IGNORED
                        ):
                    # its neither DNA nor protein, use the heavy atoms
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


def pdb2fastadic(pdb_f):
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
            ("DT", "T"),
            ]
        )
    seq_dic = {}

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


def get_align(method, lovoalign_exec):
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
        align_func = partial(
            align_strct,
            lovoalign_exec=lovoalign_exec
            )
    elif method == "sequence":
        align_func = partial(align_seq)
    else:
        available_alns = ("sequence", "structure")
        raise ValueError(
            f"Alignment method {method!r} not recognized. "
            f"Available options are {', '.join(available_alns)}"
            )
    return align_func


def align_strct(reference, model, output_path, lovoalign_exec=None):
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
            "Structural alignment needs LovoAlign "
            "get it at github.com/m3g/lovoalign"
            )
        raise ALIGNError("Path to LovoAlign executable required.")

    if not lovoalign_exec:
        raise ALIGNError("lovoalign_exec parameter not defined ")

    if not os.access(lovoalign_exec, os.X_OK):
        raise ALIGNError(f"{lovoalign_exec!r} for LovoAlign is not executable")

    numbering_dic = {}
    protein_a_dic = dict(
        (str(e.stem).split("_")[-1], e) for e in split_by_chain(reference)
        )
    protein_b_dic = dict(
        (str(e.stem).split("_")[-1], e) for e in split_by_chain(model)
        )

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
                _msg = (
                    f"LovoAlign could not read {failed_pdb} " "is it a ligand?"
                    )
                log.warning(_msg)
                alignment_pass = False

                for elem in [k for k in pa_seqdic[chain]]:
                    numbering_dic[chain][elem] = elem

        if not alignment_pass:
            # This alignment failed, move on to the next
            log.warning(
                f"Skipping alignment of chain {chain}, "
                "used sequential matching"
                )
            continue

        aln_l = lovoalign_out[alignment_start_index:alignment_end_index]

        # dump this alignment to a file
        aln_fname = Path(output_path, f"lovoalign_{chain}.aln")
        log.debug(f"Writing alignment to {aln_fname.name}")
        with open(aln_fname, "w") as fh:
            fh.write(os.linesep.join(aln_l))

        # remove the line between the alignment segments
        alignment = [aln_l[i: i + 3][:2] for i in range(0, len(aln_l), 3)]
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
                f"\"Structural\" identity of chain {chain} is {identity:.2f}%,"
                " please check the results carefully"
                )
        else:
            log.info(
                f"\"Structural\" identity of chain {chain} is {identity:.2f}%"
                )

        # logging.debug("Reading alignment and matching numbering")
        for element in alignment:
            line_a, line_b = element

            resnum_a, seq_a, _ = line_a.split()
            resnum_b, seq_b, _ = line_b.split()

            resnum_a = int(resnum_a) - 1
            resnum_b = int(resnum_b) - 1

            for resname_a, resname_b in zip(seq_a, seq_b):
                if resname_a != "-":
                    resnum_a += 1

                if resname_b != "-":
                    resnum_b += 1

                if resname_a != "-" and resname_b != "-":
                    numbering_dic[chain][resnum_b] = resnum_a

    izone_fname = Path(output_path, "lovoalign.izone")
    log.debug(f"Saving .izone to {izone_fname.name}")
    dump_as_izone(izone_fname, numbering_dic)

    return numbering_dic


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
    seqdic_ref = pdb2fastadic(reference)
    seqdic_model = pdb2fastadic(model)

    if seqdic_ref.keys() != seqdic_model.keys():
        # TODO: Implement chain-matching here
        return False

    align_dic = {}
    for ref_chain, model_chain in zip(seqdic_ref, seqdic_model):

        if ref_chain != model_chain:
            raise AlignError(
                f"Chain mismatch: {ref_chain} != {model_chain}"
                )

        align_dic[ref_chain] = {}

        seq_ref = Seq("".join(seqdic_ref[ref_chain].values()))
        seq_model = Seq("".join(seqdic_model[model_chain].values()))

        aligner = Align.PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        alns = aligner.align(seq_ref, seq_model)
        top_aln = alns[0]

        aln_fname = Path(output_path, f"blosum62_{ref_chain}.aln")
        log.debug(f"Writing alignment to {aln_fname.name}")
        with open(aln_fname, "w") as fh:
            fh.write(str(top_aln))
        aligned_ref_segment, aligned_model_segment = top_aln.aligned

        # this should always be true
        assert len(aligned_ref_segment) == len(aligned_model_segment)

        identity = (
            str(top_aln).count("|") / float(min(len(seq_ref), len(seq_model)))
            ) * 100

        if not any(e for e in top_aln.aligned):
            # No alignment!
            log.warning(
                f"No alignment for chain {ref_chain} is it protein/dna? "
                "Matching sequentially"
                )
            if all(
                    "X" in s for s in seq_ref) and all(
                    "X" in s for s in seq_model):
                # this sequence contains only ligands, do it manually
                if len(seq_ref) != len(seq_model):
                    # we cannot handle this
                    raise f"Cannot align chain {model_chain}"
                for ref_res, model_res in zip(
                        seqdic_ref[ref_chain],
                        seqdic_model[model_chain]):

                    align_dic[ref_chain].update({model_res: ref_res})
        else:
            if identity <= 40.0:
                # Identity is very low
                log.warning(
                    f"Sequence identity of chain {ref_chain} is "
                    f"{identity:.2f}%, please check the results carefully")
                log.warning(
                    "Please use alignment_method = \"structure\" instead")
            else:
                log.debug(
                    f"Sequence identity between chain {ref_chain} "
                    f" of {reference} and {model} is "
                    f"{identity:.2f}%")
            for ref_segment, model_segment in zip(
                    aligned_ref_segment, aligned_model_segment):

                start_ref_segment, end_ref_segment = ref_segment
                start_model_segment, end_model_segment = model_segment

                reslist_ref = list(seqdic_ref[ref_chain].keys())[
                    start_ref_segment:end_ref_segment]

                reslist_model = list(seqdic_model[model_chain].keys())[
                    start_model_segment:end_model_segment]

                for _ref_res, _model_res in zip(reslist_ref, reslist_model):
                    align_dic[ref_chain].update({_model_res: _ref_res})

    izone_fname = Path(output_path, "blosum62.izone")
    log.debug(f"Saving .izone to {izone_fname.name}")
    dump_as_izone(izone_fname, align_dic)

    return align_dic


def make_range(chain_range_dic):
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
    chain_ranges = {}
    for chain in chain_range_dic:
        min_idx = min(chain_range_dic[chain])
        max_idx = max(chain_range_dic[chain])
        chain_ranges[chain] = (min_idx, max_idx)
    return chain_ranges


def dump_as_izone(fname, numbering_dic):
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
                #
                izone_str = (
                    "ZONE "
                    f"{chain}{bound_res}:{chain}{unbound_res}"
                    f"{os.linesep}"
                    )
                fh.write(izone_str)


class AlignError(Exception):
    """Raised when something goes wrong with the Alignment library."""

    def __init__(self, msg=""):
        self.msg = msg
        super().__init__(self.msg)
