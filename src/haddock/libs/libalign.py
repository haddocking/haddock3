
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
from haddock.libs.libio import pdb_path_exists
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
    "TRP": ["C", "N", "CA", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],  # noqa: E501
    "TYR": ["C", "N", "CA", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],  # noqa: E501
    "VAL": ["C", "N", "CA", "O", "CB", "CG1", "CG2"]
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
    "DA": ["P", "O1P", "O2P", "O5'", "C5'", "C4'", "O4'", "C1'", "N9", "C4",
           "N3", "C2", "N1", "C6", "N6", "C5", "N7", "C8", "C2'", "C3'", "O3'"],  # noqa: E501
    "DG": ["P", "O1P", "O2P", "O5'", "C5'", "C4'", "O4'", "C1'", "N9", "C4",
           "N3", "C2", "N2", "N1", "C6", "O6", "C5", "N7", "C8", "C2'", "C3'", "O3'"],  # noqa: E501
    "DC": ["P", "O1P", "O2P", "O5'", "C5'", "C4'", "O4'", "C1'", "N1", "C6",
           "C2", "O2", "N3", "C4", "N4", "C5", "C2'", "C3'", "O3'"],
    "DT": ["P", "O1P", "O2P", "O5'", "C5'", "C4'", "O4'", "C1'", "N1", "C6",
           "C2", "O2", "N3", "C4", "O4", "C5", "C7", "C2'", "C3'", "O3'"]
    }

RNA_FULL_DICT = {
    "A": ["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'",
          "O2'", "C1'", "N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"],  # noqa: E501
    "G": ["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'",
          "O2'", "C1'", "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"],  # noqa: E501
    "C": ["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'",
          "O2'", "C1'", "N9", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4", "N4"],  # noqa: E501
    "U": ["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'",
          "O2'", "C1'", "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"]
    }


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


def load_coords(pdb_f, atoms, filter_resdic=None, numbering_dic=None, model2ref_chain_dict=None):
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
    print(f"load_coords on pdb_f {pdb_f} with chain dict {model2ref_chain_dict}")
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
                #print(f"identifier {identifier}")
                if atom_name not in atoms[resname]:
                    #print(f"not in atoms[resname] {atoms[resname]}")
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
    #print(f"coord_dic {coord_dic} chain_ranges {chain_ranges}")
    return coord_dic, chain_ranges


def get_atoms(pdb, full=False):
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
    atom_dic.update(dict((r, RNA_ATOMS) for r in RNA_RES))
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
    log.info(f"Using {method} alignment")
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

def write_alignment(top_aln, output_path, ref_chain):
    """
    Write the alignment to a file.

    Parameters
    ----------
    top_aln : Bio.Align.PairwiseAlignments
        alignment object
    
    ref_chain : str
        reference chain
    """
    aln_fname = Path(output_path, f"blosum62_{ref_chain}.aln")
    log.debug(f"Writing alignment to {aln_fname.name}")
    with open(aln_fname, "w") as fh:
        fh.write(str(top_aln))
    return aln_fname


def sequence_alignment(seq_ref, seq_model):
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    alns = aligner.align(seq_ref, seq_model)
    top_aln = alns[0]

    identity = (
                str(top_aln).count("|") / float(min(len(seq_ref), len(seq_model)))
                ) * 100
    
    aligned_ref_segment, aligned_model_segment = top_aln.aligned
    #print(f"aligned_ref_segment {aligned_ref_segment} and aligned_model_segment {aligned_model_segment}")
    # this should always be true
    assert len(aligned_ref_segment) == len(aligned_model_segment)

    return identity, top_aln, aligned_ref_segment, aligned_model_segment


#def postprocess_alignment(top_aln, identity, seq_ref, seq_model, ref_chain, model_chain, aligned_ref_segment, aligned_model_segment, align_dic, seqdic_ref, seqdic_model):
#    print(f"running postprocess_alignment for {ref_chain} - {model_chain}")
#    align_dic[ref_chain] = {}
#    if not any(e for e in top_aln.aligned):
#        # No alignment!
#        log.warning(
#            f"No alignment for chain {ref_chain} is it protein/dna-rna? "
#            "Matching sequentially"
#            )
#        if all(
#                "X" in s for s in seq_ref) and all(
#                "X" in s for s in seq_model):
#            # this sequence contains only ligands, do it manually
#            if len(seq_ref) != len(seq_model):
#                # we cannot handle this
#                # FIXME: This should raise a proper exception instead
#                raise f"Cannot align chain {model_chain}"  # noqa: B016
#            for ref_res, model_res in zip(
#                    seqdic_ref[ref_chain],
#                    seqdic_model[model_chain]):
#                align_dic[ref_chain].update({model_res: ref_res})
#    else:
#        if identity <= 40.0:
#            # Identity is very low
#            log.warning(
#                f"Sequence identity of chain {ref_chain} is "
#                f"{identity:.2f}%, please check the results carefully")
#            # log.warning(
#            #    "Please use alignment_method = \"structure\" instead")
#        else:
#            log.debug(
#                f"Sequence identity between chain {ref_chain} "
#                f" of reference and {model_chain} of model is "
#                f"{identity:.2f}%")
#        for ref_segment, model_segment in zip(
#                aligned_ref_segment, aligned_model_segment):
#            print(f"ref_segment {ref_segment} and model_segment {model_segment}")
#            start_ref_segment, end_ref_segment = ref_segment
#            start_model_segment, end_model_segment = model_segment
#            reslist_ref = list(seqdic_ref[ref_chain].keys())[
#                start_ref_segment:end_ref_segment]
#            reslist_model = list(seqdic_model[model_chain].keys())[
#                start_model_segment:end_model_segment]
#            for _ref_res, _model_res in zip(reslist_ref, reslist_model):
#                align_dic[ref_chain].update({_model_res: _ref_res})

class SeqAlign:
    """SeqAlign class."""

    def __init__(
            self,
            ):
        """
        Initialize the class.

        Parameters
        ----------
        """
        self.align_dic = {}
        self.model2ref_chain_dict = {}
        self.seqdic_ref = None
        self.seqdic_model = None
        self.seqs_ref = {}
        self.seqs_model = {}
        self.identities = []
        self.aligned_model_segments = []
        self.aligned_ref_segments = []
        self.top_alns = []

    def postprocess_alignment(
        self,
        ref_chain,
        model_chain,
        align_id
    ):
        """
        Postprocess the alignment.

        Parameters
        ----------
        ref_chain : str
            reference chain
        
        model_chain : str
            model chain
        
        align_id : int
            alignment id (index of the alignment)
        """
        print(f"running postprocess_alignment for {ref_chain} - {model_chain}")
        self.align_dic[ref_chain] = {}
        if not any(e for e in self.top_alns[align_id].aligned):
            # No alignment!
            log.warning(
                f"No alignment for chain {ref_chain} is it protein/dna-rna? "
                "Matching sequentially"
                )
            if all(
                    "X" in s for s in self.seq_refs[align_id]) and all(
                    "X" in s for s in self.seq_model[align_id]):
                # this sequence contains only ligands, do it manually
                if len(self.seqs_ref[ref_chain]) != len(self.seqs_model[model_chain]):
                    # we cannot handle this
                    # FIXME: This should raise a proper exception instead
                    raise f"Cannot align chain {model_chain}"  # noqa: B016
                for ref_res, model_res in zip(
                        self.seqdic_ref[ref_chain],
                        self.seqdic_model[model_chain]):
                    self.align_dic[ref_chain].update({model_res: ref_res})
        else:
            identity = self.identities[align_id] 
            if identity <= 40.0:
                # Identity is very low
                log.warning(
                    f"Sequence identity of chain {ref_chain} is "
                    f"{identity:.2f}%, please check the results carefully")
                # log.warning(
                #    "Please use alignment_method = \"structure\" instead")
            else:
                log.debug(
                    f"Sequence identity between chain {ref_chain} "
                    f" of reference and {model_chain} of model is "
                    f"{identity:.2f}%")
            for ref_segment, model_segment in zip(
                    self.aligned_ref_segments[align_id], self.aligned_model_segments[align_id]):
                print(f"ref_segment {ref_segment} and model_segment {model_segment}")
                start_ref_segment, end_ref_segment = ref_segment
                start_model_segment, end_model_segment = model_segment
                reslist_ref = list(self.seqdic_ref[ref_chain].keys())[
                    start_ref_segment:end_ref_segment]
                reslist_model = list(self.seqdic_model[model_chain].keys())[
                    start_model_segment:end_model_segment]
                for _ref_res, _model_res in zip(reslist_ref, reslist_model):
                    self.align_dic[ref_chain].update({_model_res: _ref_res})

    

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
    print(f"running align_seq on {reference} and {model}")
    SeqAlign_obj = SeqAlign()
    SeqAlign_obj.seqdic_ref = pdb2fastadic(reference)
    SeqAlign_obj.seqdic_model = pdb2fastadic(model)
    model2ref_chain_dict = {}
    
    # assign sequences
    seqs_ref = {}
    seqs_model = {}
    for ref_chain in SeqAlign_obj.seqdic_ref.keys():
        SeqAlign_obj.seqs_ref[ref_chain] = Seq("".join(SeqAlign_obj.seqdic_ref[ref_chain].values()))
    for model_chain in SeqAlign_obj.seqdic_model.keys():
        SeqAlign_obj.seqs_model[model_chain] = Seq("".join(SeqAlign_obj.seqdic_model[model_chain].values()))
    
    # check if chain ids match
    if SeqAlign_obj.seqdic_ref.keys() != SeqAlign_obj.seqdic_model.keys():
        # they do not match, we need to do chain matching
        n_partners = min(len(SeqAlign_obj.seqdic_ref.keys()), len(SeqAlign_obj.seqdic_model.keys()))

        # unique combinations of chains
        combs = []
        for ref_key in SeqAlign_obj.seqdic_ref.keys():
            for mod_key in SeqAlign_obj.seqdic_model.keys():
                combs.append((ref_key, mod_key))
        print(f"combs {combs}")
        identities = []
        aligned_model_segments = []
        aligned_ref_segments = []
        top_alns = []
        # align
        for ref_chain, model_chain in combs:
            print(f"ref_chain {ref_chain} and model_chain {model_chain}")

            identity, top_aln, aligned_ref_segment, aligned_model_segment = sequence_alignment(SeqAlign_obj.seqs_ref[ref_chain], SeqAlign_obj.seqs_model[model_chain])
            #print(f"identity {identity} aligned_ref_segment {aligned_ref_segment} aligned_model_segment {aligned_model_segment}")
            identities.append(identity)
            aligned_model_segments.append(aligned_model_segment)
            aligned_ref_segments.append(aligned_ref_segment)
            top_alns.append(top_aln)
        matches = 0
        align_dic = {}
        
        while matches < n_partners:
            # get the best alignment
            max_identity = max(identities)
            max_identity_index = identities.index(max_identity)
            # assignin chains
            ref_chain, model_chain = combs[max_identity_index]
            SeqAlign_obj.model2ref_chain_dict[model_chain] = ref_chain
            SeqAlign_obj.aligned_ref_segments.append(aligned_ref_segments[max_identity_index])
            SeqAlign_obj.aligned_model_segments.append(aligned_model_segments[max_identity_index])
            SeqAlign_obj.identities.append(identities[max_identity_index])
            SeqAlign_obj.top_alns.append(top_alns[max_identity_index])
            
            # writing the alignment
            write_alignment(top_alns[max_identity_index], output_path, ref_chain)

            print(f"ref_chain {ref_chain} and model_chain {model_chain} have highest identity")
            # postprocess alignment
            SeqAlign_obj.postprocess_alignment(ref_chain, model_chain, matches)
            # update identities to avoid double matches
            identities = [identities[n] if combs[n][0] != ref_chain and combs[n][1] != model_chain else -1 for n in range(len(combs))]
            matches += 1
        return SeqAlign_obj.align_dic, SeqAlign_obj.model2ref_chain_dict
    else:
        # chains do match. no need to do chain matching
        matches = 0
        for ref_chain, model_chain in zip(SeqAlign_obj.seqdic_ref, SeqAlign_obj.seqdic_model):
            SeqAlign_obj.model2ref_chain_dict[model_chain] = ref_chain
            if ref_chain != model_chain:
                raise AlignError(
                    f"Chain mismatch: {ref_chain} != {model_chain}"
                    )

            #align_dic[ref_chain] = {}

            seq_ref = SeqAlign_obj.seqs_ref[ref_chain]
            seq_model = SeqAlign_obj.seqs_model[model_chain]
            identity, top_aln, aligned_ref_segment, aligned_model_segment = sequence_alignment(seq_ref, seq_model)
            SeqAlign_obj.identities.append(identity)
            SeqAlign_obj.aligned_model_segments.append(aligned_model_segment)
            SeqAlign_obj.aligned_ref_segments.append(aligned_ref_segment)
            SeqAlign_obj.top_alns.append(top_aln)
            # write alignment
            write_alignment(top_aln, output_path, ref_chain)
            # postprocess alignment
            SeqAlign_obj.postprocess_alignment(ref_chain, model_chain, matches)
            matches += 1
            #if not any(e for e in top_aln.aligned):
            #    # No alignment!
            #    log.warning(
            #        f"No alignment for chain {ref_chain} is it protein/dna-rna? "
            #        "Matching sequentially"
            #        )
            #    if all(
            #            "X" in s for s in seq_ref) and all(
            #            "X" in s for s in seq_model):
            #        # this sequence contains only ligands, do it manually
            #        if len(seq_ref) != len(seq_model):
            #            # we cannot handle this
            #            # FIXME: This should raise a proper exception instead
            #            raise f"Cannot align chain {model_chain}"  # noqa: B016
            #        for ref_res, model_res in zip(
            #                seqdic_ref[ref_chain],
            #                seqdic_model[model_chain]):
#
            #            align_dic[ref_chain].update({model_res: ref_res})
            #else:
            #    if identity <= 40.0:
            #        # Identity is very low
            #        log.warning(
            #            f"Sequence identity of chain {ref_chain} is "
            #            f"{identity:.2f}%, please check the results carefully")
            #        log.warning(
            #            "Please use alignment_method = \"structure\" instead")
            #    else:
            #        log.debug(
            #            f"Sequence identity between chain {ref_chain} "
            #            f" of reference and {model} is "
            #            f"{identity:.2f}%")
            #    for ref_segment, model_segment in zip(
            #            aligned_ref_segment, aligned_model_segment):
#
            #        start_ref_segment, end_ref_segment = ref_segment
            #        start_model_segment, end_model_segment = model_segment
#
            #        reslist_ref = list(seqdic_ref[ref_chain].keys())[
            #            start_ref_segment:end_ref_segment]
#
            #        reslist_model = list(seqdic_model[model_chain].keys())[
            #            start_model_segment:end_model_segment]
#
            #        for _ref_res, _model_res in zip(reslist_ref, reslist_model):
            #            align_dic[ref_chain].update({_model_res: _ref_res})

        izone_fname = Path(output_path, "blosum62.izone")
        log.debug(f"Saving .izone to {izone_fname.name}")
        dump_as_izone(izone_fname, SeqAlign_obj.align_dic)

        return SeqAlign_obj.align_dic, SeqAlign_obj.model2ref_chain_dict


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
