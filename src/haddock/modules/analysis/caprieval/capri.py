"""CAPRI module."""
import os
import shlex
import shutil
import subprocess
import tempfile
from functools import partial
from pathlib import Path

import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from fccpy import read_pdb
from fccpy.contacts import get_intermolecular_contacts
from pdbtools import pdb_segxchain

from haddock import log
from haddock.libs.libontology import PDBFile
from haddock.libs.libpdb import split_by_chain


IGNORE_RES = ["SHA"]

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


class CAPRI:
    """CAPRI class."""

    def __init__(self,
                 reference,
                 model_list,
                 receptor_chain,
                 ligand_chain,
                 aln_method,
                 path,
                 **params):
        self.reference = reference
        self.model_list = []
        self.irmsd_dic = {}
        self.lrmsd_dic = {}
        self.ilrmsd_dic = {}
        self.fnat_dic = {}
        self.atoms = get_atoms(model_list)
        self.r_chain = receptor_chain
        self.l_chain = ligand_chain
        self.score_dic = {}
        self.path = path

        # TODO: For scoring we might need to get one alignment per model
        reference = str(reference)
        model = model_list[0].rel_path
        align_func = get_align(aln_method, **params)
        self.numbering_dic = align_func(reference, model, path)
        if not self.numbering_dic:
            raise CAPRIError("Could not align reference and model")

        # Load the models in the class
        for struct in model_list:
            pdb_f = struct.rel_path
            pdb_w_chain = self.add_chain_from_segid(pdb_f)
            self.model_list.append(pdb_w_chain)
            self.score_dic[pdb_f] = struct.score

    def irmsd(self, cutoff=5.0):
        """Calculate the I-RMSD."""
        # Identify reference interface
        ref_interface_resdic = self.identify_interface(self.reference, cutoff)

        # Load interface coordinates
        ref_coord_dic, _ = self.load_coords(
            self.reference, ref_interface_resdic, match=False)

        for model in self.model_list:

            mod_coord_dic, _ = self.load_coords(
                model, ref_interface_resdic, match=True)

            # Here _coord_dic keys are matched
            #  and formatted as (chain, resnum, atom)
            #  we will use atoms that are present in both
            P = []
            Q = []
            for k in ref_coord_dic.keys() & mod_coord_dic.keys():
                ref_xyz = ref_coord_dic[k]
                mod_xyz = mod_coord_dic[k]

                Q.append(ref_xyz)
                P.append(mod_xyz)

            Q = np.asarray(Q)
            P = np.asarray(P)
            # write_coords('model.pdb', P)
            # write_coords('ref.pdb', Q)

            Q = Q - self.centroid(Q)
            P = P - self.centroid(P)
            U = self.kabsch(P, Q)
            P = np.dot(P, U)
            i_rmsd = self.calc_rmsd(P, Q)
            # write_coords('model_aln.pdb', P)
            # write_coords('ref_aln.pdb', Q)

            self.irmsd_dic[model] = i_rmsd

        return self.irmsd_dic

    def lrmsd(self):
        """Calculate the L-RMSD."""
        ref_coord_dic, _ = self.load_coords(self.reference)

        for model in self.model_list:

            mod_coord_dic, _ = self.load_coords(model, match=True)

            Q = []
            P = []
            # Note: this MUST be sorted since we will use the indexes to
            #  separate between receptor and ligand coordinates
            instersection = sorted(ref_coord_dic.keys() & mod_coord_dic.keys())

            chain_ranges = {}
            for i, segment in enumerate(instersection):
                chain, _, _ = segment
                if chain not in chain_ranges:
                    chain_ranges[chain] = []
                chain_ranges[chain].append(i)

            chain_ranges = make_range(chain_ranges)
            r_start, r_end = chain_ranges[self.r_chain]
            l_start, l_end = chain_ranges[self.l_chain]

            for k in instersection:
                ref_xyz = ref_coord_dic[k]
                mod_xyz = mod_coord_dic[k]

                Q.append(ref_xyz)
                P.append(mod_xyz)

            Q = np.asarray(Q)
            P = np.asarray(P)

            # write_coord_dic('ref.pdb', ref_coord_dic)
            # write_coord_dic('model.pdb', mod_coord_dic)

            # write_coords('ref.pdb', Q)
            # write_coords('model.pdb', P)

            # # move to the origin
            Q = Q - self.centroid(Q)
            P = P - self.centroid(P)

            # get receptor coordinates
            Q_r = Q[r_start:r_end - 1]
            P_r = P[r_start:r_end - 1]

            # Center receptors and get rotation matrix
            # Q_r = Q_r - self.centroid(Q_r)
            # P_r = P_r - self.centroid(P_r)

            U_r = self.kabsch(P_r, Q_r)

            # Center complexes at receptor centroids
            Q = Q - self.centroid(Q_r)
            P = P - self.centroid(P_r)

            # Apply rotation to complex
            #  - complex are now aligned by the receptor
            P = np.dot(P, U_r)

            # write_coords('ref.pdb', Q)
            # write_coords('model.pdb', P)

            # Identify the ligand coordinates
            Q_l = Q[l_start:l_end - 1]
            P_l = P[l_start:l_end - 1]

            # write_coords('ref_l.pdb', Q_l)
            # write_coords('model_l.pdb', P_l)

            # Calculate the RMSD of the ligands
            l_rmsd = self.calc_rmsd(P_l, Q_l)

            # write_coords('ref.pdb', Q)
            # write_coords('model.pdb', P)

            self.lrmsd_dic[model] = l_rmsd

        return self.lrmsd_dic

    def ilrmsd(self, cutoff=10.0):
        """Calculate the Interface Ligand RMSD."""
        # Identify interface
        ref_interface_resdic = self.identify_interface(self.reference, cutoff)

        # Load interface coordinates
        ref_coord_dic, _ = self.load_coords(self.reference, match=False)

        ref_int_coord_dic, _ = self.load_coords(
            self.reference, ref_interface_resdic, match=False
            )

        for model in self.model_list:

            mod_coord_dic, _ = self.load_coords(model, match=True)

            mod_int_coord_dic, _ = self.load_coords(
                model, ref_interface_resdic, match=True
                )

            # write_coord_dic('ref.pdb', ref_int_coord_dic)
            # write_coord_dic('model.pdb', mod_int_coord_dic)

            # find atoms present in both interfaces
            Q_int = []
            P_int = []
            for k in sorted(ref_int_coord_dic.keys()
                            & mod_int_coord_dic.keys()):
                ref_xyz = ref_int_coord_dic[k]
                mod_xyz = mod_int_coord_dic[k]

                Q_int.append(ref_xyz)
                P_int.append(mod_xyz)

            Q_int = np.asarray(Q_int)
            P_int = np.asarray(P_int)

            # write_coords('ref.pdb', Q_int)
            # write_coords('model.pdb', P_int)

            # find atoms present in both molecules
            Q = []
            P = []
            intersection = sorted(ref_coord_dic.keys() & mod_coord_dic.keys())
            chain_ranges = {}
            for i, segment in enumerate(intersection):
                chain, _, _ = segment
                if chain not in chain_ranges:
                    chain_ranges[chain] = []
                chain_ranges[chain].append(i)

            chain_ranges = make_range(chain_ranges)
            l_start, l_end = chain_ranges[self.l_chain]

            for k in sorted(ref_coord_dic.keys() & mod_coord_dic.keys()):
                ref_xyz = ref_coord_dic[k]
                mod_xyz = mod_coord_dic[k]

                Q.append(ref_xyz)
                P.append(mod_xyz)

            Q = np.asarray(Q)
            P = np.asarray(P)

            # write_coords('ref.pdb', Q)
            # write_coords('model.pdb', P)

            # put system at origin
            Q_int = Q_int - self.centroid(Q_int)
            P_int = P_int - self.centroid(P_int)

            # # put interfaces at the origin
            # Q_int = Q_int - self.centroid(Q_int)
            # P_int = P_int - self.centroid(P_int)

            # find the rotation that minimizes the interface rmsd
            U_int = self.kabsch(P_int, Q_int)
            P_int = np.dot(P_int, U_int)

            # write_coords('ref.pdb', Q_int)
            # write_coords('model.pdb', P_int)

            # # move the system to the centroid of the interfaces
            Q = Q - self.centroid(Q)
            P = P - self.centroid(P)

            # write_coords('ref_1.pdb', Q)
            # write_coords('model_1.pdb', P)

            Q = Q - self.centroid(Q_int)
            P = P - self.centroid(P_int)

            # write_coords('ref_2.pdb', Q)
            # write_coords('model_2.pdb', P)

            # apply this rotation to the model
            #  - complexes are now aligned by the interfaces
            P = np.dot(P, U_int)

            # write_coords('ref_i.pdb', Q)
            # write_coords('model_i.pdb', P)

            # Calculate the rmsd of the ligand
            Q_l = Q[l_start:l_end - 1]
            P_l = P[l_start:l_end - 1]

            # write_coords('ref_l.pdb', Q_l)
            # write_coords('model_l.pdb', P_l)

            # this will be the interface-ligand-rmsd
            i_l_rmsd = self.calc_rmsd(P_l, Q_l)
            self.ilrmsd_dic[model] = i_l_rmsd

        return self.ilrmsd_dic

    def fnat(self, cutoff=5.0):
        """Calculate the frequency of native contacts."""
        ref_contacts = self.load_contacts(self.reference, cutoff)
        for model in self.model_list:
            model_contacts = self.load_contacts(model, cutoff)
            intersection = ref_contacts & model_contacts
            fnat = len(intersection) / float(len(ref_contacts))
            self.fnat_dic[model] = fnat
        return self.fnat_dic

    def output(self, output_f, sortby_key, sort_ascending, rankby_key,
               rank_ascending):
        """Output the CAPRI results to a .tsv file."""
        output_l = []
        for model in self.model_list:
            data = {}
            # keep always 'model' the first key
            data["model"] = model
            # create the empty rank here so that it will appear
            #  as the second column
            data["rank"] = None
            data["score"] = self.score_dic[model]
            if model in self.irmsd_dic:
                data["irmsd"] = self.irmsd_dic[model]
            if model in self.fnat_dic:
                data["fnat"] = self.fnat_dic[model]
            if model in self.lrmsd_dic:
                data["lrmsd"] = self.lrmsd_dic[model]
            if model in self.ilrmsd_dic:
                data["ilrmsd"] = self.ilrmsd_dic[model]
            # list of dictionaries
            output_l.append(data)

        # Get the ranking of each model
        rankkey_values = [(i, k[rankby_key]) for i, k in enumerate(output_l)]
        rankkey_values.sort(key=lambda x: x[1], reverse=not rank_ascending)
        for i, k in enumerate(rankkey_values, start=1):
            idx, _ = k
            output_l[idx]["rank"] = i

        # Sort the column
        key_values = [(i, k[sortby_key]) for i, k in enumerate(output_l)]
        key_values.sort(key=lambda x: x[1], reverse=not sort_ascending)

        max_model_space = max(len(str(_d["model"])) for _d in output_l) + 2
        hmodel = "model".center(max_model_space, " ")
        header = hmodel + "".join(
            _.rjust(10, " ") for _ in list(output_l[0].keys())[1:]
            )

        with open(output_f, "w") as out_fh:
            out_fh.write(header + os.linesep)
            for idx, _ in key_values:
                row_l = []
                for value in output_l[idx].values():
                    if isinstance(value, Path):
                        row_l.append(str(value).ljust(max_model_space, " "))
                    elif isinstance(value, int):
                        row_l.append(f"{value}".rjust(10, " "))
                    else:
                        row_l.append(f"{value:.3f}".rjust(10, " "))
                out_fh.write("".join(row_l) + os.linesep)

    @staticmethod
    def identify_interface(pdb_f, cutoff=5.0):
        """Identify the interface."""
        pdb = read_pdb(pdb_f)
        interface_resdic = {}
        for atom_i, atom_j in get_intermolecular_contacts(pdb, cutoff):

            if atom_i.chain not in interface_resdic:
                interface_resdic[atom_i.chain] = []
            if atom_j.chain not in interface_resdic:
                interface_resdic[atom_j.chain] = []

            if atom_i.resid not in interface_resdic[atom_i.chain]:
                interface_resdic[atom_i.chain].append(atom_i.resid)
            if atom_j.resid not in interface_resdic[atom_j.chain]:
                interface_resdic[atom_j.chain].append(atom_j.resid)

        return interface_resdic

    @staticmethod
    def load_contacts(pdb_f, cutoff=5.0):
        """Load residue-based contacts."""
        con_list = []
        structure = read_pdb(pdb_f)
        for atom_i, atom_j in get_intermolecular_contacts(structure, cutoff):
            con = (atom_i.chain, atom_i.resid, atom_j.chain, atom_j.resid)
            con_list.append(con)
        return set(con_list)

    @staticmethod
    def calc_rmsd(V, W):
        """Calculate the RMSD from two vectors."""
        diff = np.array(V) - np.array(W)
        N = len(V)
        return np.sqrt((diff * diff).sum() / N)

    @staticmethod
    def kabsch(P, Q):
        """Find the rotation matrix using Kabsch algorithm."""
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

    @staticmethod
    def centroid(X):
        """Get the centroid."""
        X = np.array(X)
        return X.mean(axis=0)

    @staticmethod
    def add_chain_from_segid(pdb_path):
        """Replace the chainID with the segID."""
        temp_f = tempfile.NamedTemporaryFile(delete=False, mode="w+t")
        with open(pdb_path) as fh:
            for line in list(pdb_segxchain.run(fh)):
                temp_f.writelines(line)
        temp_f.close()
        # REPLACE!
        new_pdb_path = shutil.move(temp_f.name, pdb_path)
        return new_pdb_path

    def load_coords(self, pdb_f, filter_resdic=None, match=False):
        """Load coordinates from PDB."""
        coord_dic = {}
        chain_dic = {}
        idx = 0
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

                    if match:
                        try:
                            resnum = self.numbering_dic[chain][resnum]
                        except KeyError:
                            # this residue is not matched, and so it should
                            #  not be considered
                            # self.log(
                            #     f'WARNING: {chain}.{resnum}.{atom_name}'
                            #     ' was not matched!'
                            #     )
                            continue

                    # identifier = f'{chain}.{resnum}.{atom_name}'
                    identifier = (chain, resnum, atom_name)

                    if atom_name not in self.atoms[resname]:
                        continue

                    if chain not in chain_dic:
                        chain_dic[chain] = []

                    if filter_resdic:
                        # Only retrieve coordinates from the filter_resdic
                        if (chain in filter_resdic
                                and resnum in filter_resdic[chain]):
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
            min_idx = min(chain_dic[chain])
            max_idx = max(chain_dic[chain])
            chain_ranges[chain] = (min_idx, max_idx)

        return coord_dic, chain_ranges


class CAPRIError(Exception):
    """Raised when something goes wrong with the CAPRI class."""

    def __init__(self, msg=""):
        self.msg = msg
        super().__init__(self.msg)


def get_atoms(pdb_list):
    """Identify what is the molecule type of each PDB."""
    atom_dic = {}
    atom_dic.update(dict((r, PROT_ATOMS) for r in PROT_RES))
    atom_dic.update(dict((r, DNA_ATOMS) for r in DNA_RES))
    for pdb in pdb_list:
        if isinstance(pdb, PDBFile):
            pdb = pdb.rel_path
        with open(pdb) as fh:
            for line in fh.readlines():
                if line.startswith("ATOM"):
                    resname = line[17:20].strip()
                    atom_name = line[12:16].strip()
                    element = line[76:78].strip()
                    if (resname not in PROT_RES
                            and resname not in DNA_RES
                            and resname not in IGNORE_RES):
                        # its neither DNA nor protein, use the heavy atoms
                        # WARNING: Atoms that belong to unknown residues mutt
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
    """Write the sequence as a fasta."""
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
    with open(pdb_f) as fh:
        for line in fh.readlines():
            if line.startswith("ATOM"):
                res_num = int(line[22:26])
                res_name = line[17:20].strip()
                chain = line[21]
                if res_name in IGNORE_RES:
                    continue
                try:
                    one_letter = res_codes[res_name]
                except KeyError:
                    one_letter = "X"
                if chain not in seq_dic:
                    seq_dic[chain] = {}
                seq_dic[chain][res_num] = one_letter
    return seq_dic


def get_align(method, **kwargs):
    """Get the alignment function."""
    log.info(f"Using {method} alignment")
    if method == "structure":
        return partial(align_strct, lovoalign_exec=kwargs['lovoalign_exec'])
    elif method == "sequence":
        return partial(align_seq)
    else:
        available_alns = ("sequence", "structure")
        raise ValueError(
            f"Alignment method {method!r} not recognized. "
            f"Available options are {', '.join(available_alns)}"
            )


def align_strct(reference, model, output_path, lovoalign_exec=None):
    """Structuraly align and get numbering relationship."""
    if lovoalign_exec is None:
        log.error(
            "Structural alignment needs LovoAlign "
            "get it at github.com/m3g/lovoalign"
            )
        raise CAPRIError("Path to LovoAlign executable required.")

    if not lovoalign_exec:
        raise CAPRIError("lovoalign_exec parameter not defined ")

    if not os.access(lovoalign_exec, os.X_OK):
        raise CAPRIError(f"{lovoalign_exec!r} for LovoAlign is not executable")

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
        # logging.debug(f'Structurally aligning chain {chain}')
        numbering_dic[chain] = {}
        cmd = (
            f"{lovoalign_exec} -p1 {protein_a_dic[chain]} "
            f"-p2 {protein_b_dic[chain]} "
            f"-c1 {chain} -c2 {chain}"
            )

        # logging.debug(f'Command is: {cmd}')
        p = subprocess.run(shlex.split(cmd), capture_output=True, text=True)
        lovoalign_out = p.stdout.split(os.linesep)

        # we don't need the splitted proteins anymore
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
        alignment = [aln_l[i:i + 3][:2] for i in range(0, len(aln_l), 3)]
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
            log.info(
                f'"Structural" identity of chain {chain} is {identity:.2f}%'
                )

        # logging.debug('Reading alignment and matching numbering')
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
    """Sequence align and get the numbering relationship."""
    seqdic_a = pdb2fastadic(reference)
    seqdic_b = pdb2fastadic(model)

    if seqdic_a.keys() != seqdic_b.keys():
        # TODO: Implement chain-matching here
        return False

    align_dic = {}
    for a, b in zip(seqdic_a, seqdic_b):

        align_dic[a] = {}

        seq_a = Seq("".join(seqdic_a[a].values()))
        seq_b = Seq("".join(seqdic_b[b].values()))

        aligner = Align.PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        alns = aligner.align(seq_a, seq_b)
        top_aln = alns[0]

        aln_fname = Path(output_path, f"blosum62_{a}.aln")
        log.debug(f"Writing alignment to {aln_fname.name}")
        with open(aln_fname, "w") as fh:
            fh.write(str(top_aln))
        aligned_seg_a, aligned_seg_b = top_aln.aligned

        # this should always be true
        assert len(aligned_seg_a) == len(aligned_seg_b)

        identity = (
            str(top_aln).count("|") / float(min(len(seq_a), len(seq_b)))
            ) * 100

        if not any(e for e in top_aln.aligned):
            # No alignment!
            log.warning(
                f"No alignment for chain {a} is it protein/dna? "
                "Matching sequentially"
                )
            if all("X" in s for s in seq_a) and all("X" in s for s in seq_b):
                # this sequence contains only ligands, do it manually
                if len(seq_a) != len(seq_b):
                    # we cannot handle this
                    raise f"Cannot align chain {b}"
                for res_a, res_b in zip(seqdic_a[a], seqdic_b[b]):
                    align_dic[a].update({res_a: res_b})
        else:
            if identity <= 40.0:
                # Identity is very low
                log.warning(
                    f"Sequence identity of chain {a} is {identity:.2f}%,"
                    " please check the results carefully"
                    )
                log.warning('Please use alignment_method = "structure" instead')
            else:
                log.info(f"Sequence identity of chain {a} is {identity:.2f}%")
            for seg_a, seg_b in zip(aligned_seg_a, aligned_seg_b):
                start_a, end_a = seg_a
                start_b, end_b = seg_b

                reslist_a = list(seqdic_a[a].keys())[start_a:end_a]
                reslist_b = list(seqdic_b[b].keys())[start_b:end_b]

                align_dic[a].update(
                    dict((i, j) for i, j in zip(reslist_a, reslist_b))
                    )
    izone_fname = Path(output_path, "blosum62.izone")
    log.debug(f"Saving .izone to {izone_fname.name}")
    dump_as_izone(izone_fname, align_dic)

    return align_dic


def make_range(chain_range_dic):
    """Expand a chain dictionary into ranges."""
    chain_ranges = {}
    for chain in chain_range_dic:
        min_idx = min(chain_range_dic[chain])
        max_idx = max(chain_range_dic[chain])
        chain_ranges[chain] = (min_idx, max_idx)
    return chain_ranges


def dump_as_izone(fname, numbering_dic):
    """Dump the numbering dictionary as .izone."""
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


# # debug only
# def write_coord_dic(output_name, coord_dic):
#     """Add a dummy atom to a PDB file according to a list of coordinates."""
#     with open(output_name, "w") as fh:
#         for i, k in enumerate(coord_dic):
#             atom_num = f"{i+1}".rjust(4, " ")
#             chain, resnum, atom = k
#             resnum = int(resnum)
#             resnum = f"{resnum}".rjust(3, " ")
#             atom_name = f"{atom}".rjust(3, " ")
#             x, y, z = coord_dic[k]
#             dum_x = f"{x:.3f}".rjust(7, " ")
#             dum_y = f"{y:.3f}".rjust(7, " ")
#             dum_z = f"{z:.3f}".rjust(7, " ")
#             dummy_line = (
#                 f"ATOM   {atom_num} {atom_name}  DUM {chain} {resnum}   "
#                 f"  {dum_x} {dum_y} {dum_z}  1.00  1.00   "
#                 "        H  " + os.linesep
#                 )
#             fh.write(dummy_line)


# # debug only
# def write_coords(output_name, coor_list):
#     """Add a dummy atom to a PDB file according to a list of coordinates."""
#     with open(output_name, "w") as fh:
#         for i, dummy_coord in enumerate(coor_list):
#             atom_num = f"{i}".rjust(4, " ")
#             resnum = f"{i}".rjust(3, " ")
#             dum_x = f"{dummy_coord[0]:.3f}".rjust(7, " ")
#             dum_y = f"{dummy_coord[1]:.3f}".rjust(7, " ")
#             dum_z = f"{dummy_coord[2]:.3f}".rjust(7, " ")
#             dummy_line = (
#                 f"ATOM   {atom_num}  H   DUM X {resnum}   "
#                 f"  {dum_x} {dum_y} {dum_z}  1.00  1.00   "
#                 "        H  " + os.linesep
#                 )
#             fh.write(dummy_line)


# # debug only
# def write_pymol_viz(resdic):
#     """Write PyMol vizualitation."""
#     for k in resdic:
#         reslist = "+".join(map(str, resdic[k]))
#         cmd = f"sele {k}, chain {k} and resid {reslist}"
#         print(cmd)
