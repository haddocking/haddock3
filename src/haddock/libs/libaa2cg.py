#  CODE TAKEN FROM MARTINIZE 1.1 ##
# Please refer to it for more information
"""
Reduces complexity of protein residue to the MARTINI coarse grained model:
CA, O, Bead(s) in specific atom location.

Reference:
Monticelli et al. The MARTINI coarse-grained force field: extension to proteins. 
J. Chem. Theory Comput. (2008) vol. 4 (5) pp. 819-834

Martinize Script from Tserk Wassenaar

Uses Biopython to parse the structure and DSSP output.
Uses pieces of the martinize-1.1.py script to convert the SS types

Outputs a coarse grained pdb file (*_ss.pdb) with assigned bfactors.
Outputs a tbl file to map the beads to the atoms they represent.

Updates
 - Updated python version to 2.6 to support isdisjoint() set method in DSSP.py (JR Apr 2012)
 - Residues that DSSP can't handle (incomplete backbone f ex) treated as coil  (JR Apr 2012)
 - Update to_one_letter_code library to protein_letters_3to1 (Jorge Roel 2017)
 - Inclusion of fake beads for corresponding amino-acids <SCd> (Jorge Roel 2017)
 - Changed the mapping routine to include DNA bead types (Rodrigo Honorato 2018)
 - Implemented feature to check if nucleic acid is a candidate for hbond (Rodrigo Honorato 2018)
"""

import itertools
import os
import random
import subprocess
import warnings
from pathlib import Path
from haddock import log
import tempfile
import collections
import math

from Bio.PDB import Entity
from Bio.PDB import PDBIO
from Bio.PDB import PDBParser
from Bio.PDB.StructureBuilder import StructureBuilder

from haddock.core.exceptions import ModuleError

warnings.filterwarnings("ignore")

CRYST_LINE = "CRYST1 " + os.linesep

def norm(a):
    """

    Args:
        a:

    Returns:

    """
    return math.sqrt(norm2(a))


def norm2(a):
    """

    Args:
        a:

    Returns:

    """
    return sum([i * i for i in a])


def pat(x, c="."):
    """
    Reformats pattern strings.

    Args:
        x:
        c:

    Returns:

    """
    return x.replace(c, "\x00").split()


def hash(x, y):
    """
    Makes a dictionary from two lists.

    Args:
        x:
        y:

    Returns:

    """
    return dict(zip(x, y))


def spl(x):
    """
    Splits a string.

    Args:
        x:

    Returns:

    """
    return x.split()


def tt(program):
    """

    Args:
        program:

    Returns:

    """
    return "".join([ssd[program].get(chr(i), "C") for i in range(256)])


def typesub(seq, patterns, types):
    """
    Pattern substitutions.

    Args:
        seq:
        patterns:
        types:

    Returns:

    """
    for i, j in zip(patterns, types):
        seq = seq.replace(i, j)
    return seq


def ss_classification(ss, program="dssp"):
    """
   Translates a string encoding the secondary structure to a string of corresponding Martini types, taking the
   origin of the secondary structure into account, and replacing termini if requested.

    Args:
        ss:
        program:

    Returns:

    """
    # Translate dssp/pymol/gmx ss to Martini ss
    ss = ss.translate(sstt[program])
    # Separate the different secondary structure types
    sep = dict([(i, ss.translate(sstd[i])) for i in sstd.keys()])
    # Do type substitutions based on patterns
    # If the ss type is not in the patterns lists, do not substitute
    # (use empty lists for substitutions)

    typ = [typesub(sep[i], patterns.get(i, []), pattypes.get(i, []))
           for i in sstd.keys()]
    # Translate all types to numerical values
    typ = [[ord(j) for j in list(i)] for i in typ]
    # Sum characters back to get a full typed sequence
    typ = "".join([chr(sum(i)) for i in zip(*typ)])
    # Return both the actual as well as the fully typed sequence
    return ss, typ

# ----+--------------------------------------+
#   A | SECONDARY STRUCTURE TYPE DEFINITIONS |
# ----+--------------------------------------+


ssdefs = {
    "dssp": list(".HGIBETSC~"),  # DSSP one letter secondary structure code     #@#
    "pymol": list(".H...S...L"),  # Pymol one letter secondary structure code    #@#
    "gmx": list(".H...ETS.C"),  # Gromacs secondary structure dump code        #@#
    "self": list("FHHHEETSCC")  # Internal CG secondary structure codes        #@#
}
cgss = list("FHHHEETSCC")  # Corresponding CG secondary structure types   #@#

patterns = {
    "H": pat(".H. .HH. .HHH. .HHHH. .HHHHH. .HHHHHH. .HHHHHHH. .HHHH HHHH.")  # @#
}
pattypes = {
    "H": pat(".3. .33. .333. .3333. .13332. .113322. .1113222. .1111 2222.")  # @#
}

ss_to_code = {"C": 1,  # Free,
              " ": 1,
              "S": 2,
              "H": 3,
              "1": 4,
              "2": 5,
              "3": 6,
              "E": 7,  # Extended
              "T": 8,  # Turn
              "F": 9  # Fibril
              }

ss_eq = list("CBHHHHBTF")

# List of programs for which secondary structure definitions can be processed
programs = ssdefs.keys()

# Dictionaries mapping ss types to the CG ss types
ssd = dict([(i, hash(ssdefs[i], cgss)) for i in programs])

# The translation table depends on the program used to obtain the
# secondary structure definitions
sstt = dict([(i, tt(i)) for i in programs])

# The following translation tables are used to identify stretches of
# a certain type of secondary structure.
null = "\x00"
sstd = dict([(i, ord(i) * null + i + (255 - ord(i)) * null) for i in cgss])

# ==========================================================================================#
# ==========================================================================================#
# ==========================================================================================#

# CG MAPPING INFORMATION

bb = "CA C N O "
prot_atoms = {"ALA": [bb + "CB"],
              "CYS": [bb, "CB SG"],
              "ASP": [bb, "CB CG OD1 OD2"],
              "GLU": [bb, "CB CG CD OE1 OE2"],
              "PHE": [bb, "CB CG CD1", "CD2 CE2", "CE1 CZ"],
              "GLY": [bb],
              "HIS": [bb, "CB CG", "CD2 NE2", "ND1 CE1"],
              "ILE": [bb, "CB CG1 CG2 CD1"],
              "LYS": [bb, "CB CG CD", "CE NZ"],
              "LEU": [bb, "CB CG CD1 CD2"],
              "MET": [bb, "CB CG SD CE"],
              "ASN": [bb, "CB CG ND1 ND2 OD1 OD2"],  # ND1?
              "PRO": [bb, "CB CG CD"],
              "GLN": [bb, "CB CG CD OE1 OE2 NE1 NE2"],
              "ARG": [bb, "CB CG CD", "NE CZ NH1 NH2"],
              "SER": [bb, "CB OG"],
              "THR": [bb, "CB OG1 CG2"],
              "VAL": [bb, "CB CG1 CG2"],
              "TRP": [bb, "CB CG CD2", "CD1 NE1 CE2", "CE3 CZ3", "CZ2 CH2"],
              "TYR": [bb, "CB CG CD1", "CD2 CE2", "CE1 CZ OH"]}

bead_names = ["BB", "SC1", "SC2", "SC3", "SC4"]

# insert beads into the data structure
cg_mapping = {}
for res in prot_atoms:
    cg_mapping[res] = collections.OrderedDict()
    for i, atom_l in enumerate(prot_atoms[res]):
        bead = bead_names[i]
        cg_mapping[res][atom_l] = bead

######################################
# Nucleotide mapping,
# This is a custom naming convention
# but the atom mapping is defined in
# 10.1021/acs.jctc.5b00286 -  S1
######################################

DA_beads = collections.OrderedDict()
DA_beads["O3'* P O1P O2P O5' OP1 OP2"] = "BB1"
DA_beads["C5' O4' C4'"] = "BB2"
DA_beads["C3' C2' C1'"] = "BB3"
DA_beads["N9 C4"] = "SC1"
DA_beads["C2 N3"] = "SC2"
DA_beads["C6 N6 N1"] = "SC3"
DA_beads["C8 N7 C5"] = "SC4"

DC_beads = collections.OrderedDict()
DC_beads["O3'* P O1P O2P O5' OP1 OP2"] = "BB1"
DC_beads["C5' O4' C4'"] = "BB2"
DC_beads["C3' C2' C1'"] = "BB3"
DC_beads["N1 C6"] = "SC1"
DC_beads["N3 C2 O2"] = "SC2"
DC_beads["C5 C4 N4"] = "SC3"

DG_beads = collections.OrderedDict()
DG_beads["O3'* P O1P O2P O5' OP1 OP2"] = "BB1"
DG_beads["C5' O4' C4'"] = "BB2"
DG_beads["C3' C2' C1'"] = "BB3"
DG_beads["N9 C4"] = "SC1"
DG_beads["C2 N2 N3"] = "SC2"
DG_beads["C6 O6 N1"] = "SC3"
DG_beads["C8 N7 C5"] = "SC4"

DT_beads = collections.OrderedDict()
DT_beads["O3'* P O1P O2P O5' OP1 OP2"] = "BB1"
DT_beads["C5' O4' C4'"] = "BB2"
DT_beads["C3' C2' C1'"] = "BB3"
DT_beads["N1 C6"] = "SC1"
DT_beads["N3 C2 O2"] = "SC2"
DT_beads["C5 C4 O4 C7"] = "SC3"

A_beads = collections.OrderedDict()
A_beads["O3'* P O1P O2P O5' OP1 OP2"] = "BB1"
A_beads["C5' O4' C4'"] = "BB2"
A_beads["C3' C2' O2' C1'"] = "BB3"
A_beads["N9 C4"] = "SC1"
A_beads["C2 N3"] = "SC2"
A_beads["C6 N6 N1"] = "SC3"
A_beads["C8 N7 C5"] = "SC4"

C_beads = collections.OrderedDict()
C_beads["O3'* P O1P O2P O5' OP1 OP2"] = "BB1"
C_beads["C5' O4' C4'"] = "BB2"
C_beads["C3' C2' O2' C1'"] = "BB3"
C_beads["N1 C6"] = "SC1"
C_beads["N3 C2 O2"] = "SC2"
C_beads["C5 C4 N4"] = "SC3"

G_beads = collections.OrderedDict()
G_beads["O3'* P O1P O2P O5' OP1 OP2"] = "BB1"
G_beads["C5' O4' C4'"] = "BB2"
G_beads["C3' C2' O2' C1'"] = "BB3"
G_beads["N9 C4"] = "SC1"
G_beads["C2 N2 N3"] = "SC2"
G_beads["C6 O6 N1"] = "SC3"
G_beads["C8 N7 C5"] = "SC4"

U_beads = collections.OrderedDict()
U_beads["O3'* P O1P O2P O5' OP1 OP2"] = "BB1"
U_beads["C5' O4' C4'"] = "BB2"
U_beads["C3' C2' O2' C1'"] = "BB3"
U_beads["N1 C6"] = "SC1"
U_beads["N3 C2 O2"] = "SC2"
U_beads["C5 C4 O4"] = "SC3"

cg_mapping["DA"] = DA_beads
cg_mapping["DC"] = DC_beads
cg_mapping["DT"] = DT_beads
cg_mapping["DG"] = DG_beads

cg_mapping["A"] = A_beads
cg_mapping["C"] = C_beads
cg_mapping["U"] = U_beads
cg_mapping["G"] = G_beads

pairing = {
    ("DG", "DC"): [("N2", "O2"), ("N1", "N3"), ("O6", "N4")],
    ("DC", "DG"): [("O2", "N2"), ("N3", "N1"), ("N4", "O6")],
    ("DA", "DT"): [("N6", "O4"), ("N1", "N3")],
    ("DT", "DA"): [("O4", "N6"), ("N3", "N1")],
    #
    ("G", "C"): [("N2", "O2"), ("N1", "N3"), ("O6", "N4")],
    ("C", "G"): [("O2", "N2"), ("N3", "N1"), ("N4", "O6")],
    ("A", "U"): [("N6", "O4"), ("N1", "N3")],
    ("U", "A"): [("O4", "N6"), ("N3", "N1")],

}

polar = ["GLN", "ASN", "SER", "THR"]
charged = ["ARG", "LYS", "ASP", "GLU"]


def add_dummy(bead_list, dist=0.11, n=2):
    """

    Args:
        bead_list:
        dist:
        n:

    Returns:

    """
    new_bead_dic = {}

    # Generate a random vector in a sphere of -1 to +1, to add to the bead position
    v = [random.random() * 2. - 1, random.random() * 2. - 1, random.random() * 2. - 1]

    # Calculated the length of the vector and divide by the final distance of the dummy bead
    norm_v = norm(v) / dist

    # Resize the vector
    vn = [i / norm_v for i in v]

    # m sets the direction of the added vector, currently only works when adding one or two beads.
    m = 1
    for j in range(n):  # create two new beads
        bead_s = str(j + 1)
        new_name = f"SCD{bead_s}"  # set the name of the new bead
        new_bead_dic[new_name] = [i + (m * j) for i, j in zip(bead_list[-1][1], vn)]
        m *= -2
    return new_bead_dic


def map_cg(chain):
    """

    Args:
        chain:

    Returns:

    """
    m_dic = collections.OrderedDict()

    for aares in chain:
        m_dic[aares] = collections.OrderedDict()

        resn = aares.resname.split()[0]  # resname
        segid = aares.segid.strip()
        resi = aares.id[1]

        if resn not in cg_mapping:
            log.warning(f"Residue {resn} not in cg_mapping — skipping.")
            continue  # skip to next residue

        # for each atom segment, calculate its center of mass and map the correct bead
        for atom_segment in cg_mapping[resn]:
            atoms = [aares[a] for a in atom_segment.split() if a in aares.child_dict]

            if atoms:
                if "*" in atom_segment:  # this is important to correctly place CG DNA beads

                    # this * means it belongs to the previous residue... find it!
                    target_previous_atom_list = [a for a in atom_segment.split() if "*" in a]

                    for target_atom in target_previous_atom_list:
                        # does it exist?
                        target_atom_name = target_atom.split("*")[0]
                        try:
                            previous_atom = chain[resi - 1][target_atom_name]
                            # how far away the previous atom is from this atom segment?
                            #  if it is too far away this could be the next chain...!
                            minimum_dist = min([(a - previous_atom) for a in atoms])
                            if minimum_dist < 2.0:  # 2.0 A is very permissive
                                atoms.append(previous_atom)
                        except KeyError:
                            # previous atom not found, move on
                            pass

            if not atoms:
                log.warning('Residue {} {:d} of chain {} cannot be processed: missing atoms {} '.
                                 format(resn, resi, aares.parent.id, atom_segment))
                continue

            bead_name = cg_mapping[resn][atom_segment]

            # get center of mass
            code = list(set([a.bfactor for a in aares if a.bfactor != 0]))

            if len(code) > 1:
                emsg = "Something is wrong with HADDOCK codes"
                raise ModuleError(emsg)

            if not code:
                code = 0.0
            else:
                code = code[0]

            bead_coord = center_of_mass(atoms)

            atom_segment = " ".join(list(atom_segment.split())).replace("*", "")

            # restrain for backmapping
            restrain = "assign (segid {}CG and resid {:d} and name {})"\
                " (segid {} and resid {:d} and (name {})) 0 0 0".\
                format(segid, resi, bead_name, segid, resi, " or name ".join(atom_segment.split()))

            m_dic[aares][bead_name] = bead_coord, code, restrain

    # add dummy beads whenever its needed
    for r in m_dic:

        if r.resname in polar:
            d = 0.14  # distance
            n = 2  # number of dummy beads to be placed

        elif r.resname in charged:
            d = 0.11  # distance
            n = 1  # number of dummy beads to be placed

        else:
            continue

        # add to data structure
        # this special beads have no HADDOCK code
        bead_list = [(b, m_dic[r][b][0]) for b in m_dic[r]]
        dummy_bead_dic = add_dummy(bead_list, dist=d, n=n)
        for db in dummy_bead_dic:
            db_coords = dummy_bead_dic[db]
            # code should be the same as the residue
            code = m_dic[r][list(m_dic[r])[0]][1]

            m_dic[r][db] = (db_coords, code, None)

    return m_dic


def center_of_mass(entity, geometric=False):
    """
    Returns gravitic [default] or geometric center of mass of an Entity.
    Geometric assumes all masses are equal (geometric=True)

    Args:
        entity:
        geometric:

    Returns:

    """
    # Structure, Model, Chain, Residue
    if isinstance(entity, Entity.Entity):
        atom_list = entity.get_atoms()
    # List of Atoms
    elif hasattr(entity, "__iter__") and [x for x in entity if x.level == "A"]:
        atom_list = entity
    else:  # Some other weirdo object
        raise ValueError("Center of Mass can only be calculated from the following objects:\n"
                         "Structure, Model, Chain, Residue, list of Atoms.")

    masses = []
    positions = [[], [], []]  # [ [X1, X2, ..] , [Y1, Y2, ...] , [Z1, Z2, ...] ]

    for atom in atom_list:
        masses.append(atom.mass)

        for i, coord in enumerate(atom.coord.tolist()):
            positions[i].append(coord)

    # If there is a single atom with undefined mass complain loudly.
    if "ukn" in set(masses) and not geometric:
        raise ValueError(f"Some Atoms don't have an element assigned.{os.linesep}"
                         "Try adding them manually or calculate the geometrical "
                         "center of mass instead.")

    if geometric:
        return [sum(coord_list) / len(masses) for coord_list in positions]

    w_pos = [[], [], []]
    for atom_index, atom_mass in enumerate(masses):
        w_pos[0].append(positions[0][atom_index] * atom_mass)
        w_pos[1].append(positions[1][atom_index] * atom_mass)
        w_pos[2].append(positions[2][atom_index] * atom_mass)

    return [sum(coord_list) / sum(masses) for coord_list in w_pos]


def determine_hbonds(structure):
    """

    Args:
        structure:

    Returns:

    """
    nuc = ["DA", "DC", "DG", "DT", "A", "C", "G", "U"]
    aa = ["ALA", "CYS", "ASP", "GLU", "PHE",
          "GLY", "HIS", "ILE", "LYS", "LEU",
          "MET", "ASN", "PRO", "GLN", "ARG",
          "SER", "THR", "VAL", "TRP", "TYR"]

    pair_list = []
    for model in structure:

        dna_chain_l = []

        for chain in model:

            prot_comp = len([r for r in chain.get_residues() if r.resname.split()[0] in aa])
            dna_comp = len([r for r in chain.get_residues() if r.resname.split()[0] in nuc])

            if prot_comp:
                # protein
                pass

            if dna_comp:
                # nucleic
                dna_chain_l.append(chain)

        if len(dna_chain_l) == 1:
            log.warning('Only one DNA/RNA chain detected, is this correct?')

            chain_a = dna_chain_l[0]
            reslist_a = [r for r in chain_a.get_residues()]

            for ra, rb in itertools.combinations(reslist_a, 2):
                pair = identify_pairing(ra, rb)
                if pair:
                    pair_list.append(pair)

        if len(dna_chain_l) > 1:  # list sizes could be different, this might be improbable
            for chain_a, chain_b in itertools.combinations(dna_chain_l, 2):
                reslist_a = [r for r in chain_a.get_residues()]
                reslist_b = [r for r in chain_b.get_residues()]
                for ra in reslist_a:
                    for rb in reslist_b:
                        pair = identify_pairing(ra, rb)
                        pair_list.append(pair)
    return pair_list


def identify_pairing(ra, rb):
    """

    Args:
        ra:
        rb:

    Returns:

    """
    pair = []

    # check if the pairing is correct
    ra_name = ra.resname.split()[0]
    rb_name = rb.resname.split()[0]

    try:
        atom_pair_list = pairing[ra_name, rb_name]
    except KeyError:
        # pairing not possible
        return

    # check if distances are ok
    distance_l = []

    for atom_list in atom_pair_list:

        try:
            a = ra[atom_list[0]]
            b = rb[atom_list[1]]
            distance_l.append(a - b)
        except KeyError:
            # residue does not have the necessary sidechain atoms
            #  assume its not a pair
            return

    # check P-P distances to make sure its the opposite base
    # distances for perfect DNA:
    #  opposite = 18.8A
    #  sequential = 6.6A
    #
    # 10.0A should be sufficient
    p_cutoff = 10.0

    # Basedist_cutoff = 3.5
    basedist_cutoff = 3.5
    try:
        pa = ra.child_dict["P"]
        pb = rb.child_dict["P"]
        p_distance = pa - pb
    except KeyError:
        # some base is missing its P, use the geometric center instead
        cen_a = center_of_mass(ra.child_dict.values())
        cen_b = center_of_mass(rb.child_dict.values())
        p_distance = math.sqrt(sum([(a - b) ** 2 for a, b in zip(cen_a, cen_b)]))

    if p_distance > p_cutoff:
        resnum_a = ra.id[1]
        resnum_b = rb.id[1]

        if all(e < basedist_cutoff for e in distance_l):  # if ALL bonds are within range
            # KEEP IN MIND THAT this will not account for badly paired DNA
            # Implement a way to search for the closest possible pair? ##
            segid_a = ra.get_segid().split()[0]
            segid_b = rb.get_segid().split()[0]
            pair = (resnum_a, segid_a), (resnum_b, segid_b)

            # special atoms, mark them!
            for atom_pair in atom_pair_list:
                atom_a, atom_b = atom_pair

                ra[atom_a].bfactor = 1
                rb[atom_b].bfactor = 1
    return pair


def output_cg_restraints(pair_list):
    """

    Args:
        pair_list:

    Returns:

    """
    out = open("dna_restraints.def", "w")
    for i, e in enumerate(pair_list):
        idx = i + 1
        res_a = e[0][0]
        segid_a = e[0][1]
        res_b = e[1][0]
        segid_b = e[1][1]
        out.write(f"{{===>}} base_a_{idx}=(resid {res_a} and segid {segid_a});\n"
                  f"{{===>}} base_b_{idx}=(resid {res_b} and segid {segid_b});\n\n")
    out.close()


def extract_groups(pair_list):
    """

    Args:
        pair_list:

    Returns:

    """
    # this will be used to define AA restraints
    out = open("dna-aa_groups.dat", "w")
    # extract groups
    group_a = [a[0][0] for a in pair_list]
    segid_a = list(set([a[0][1] for a in pair_list]))

    group_b = [a[1][0] for a in pair_list]
    segid_b = list(set([a[0][1] for a in pair_list]))

    if len(segid_a) != 1:
        emsg = "Something is wrong with SEGID A"
        raise ModuleError(emsg)

    if len(segid_b) != 1:
        emsg = "Something is wrong with SEGID B"
        raise ModuleError(emsg)

    segid_a = segid_a[0]
    segid_b = segid_b[0]

    group_a.sort()
    group_b.sort()
    out.write(f"{group_a[0]}:{group_a[-1]}\n{segid_a}\n{group_b[0]}:{group_b[-1]}\n{segid_b}")
    out.close()


def create_file_with_cryst(pdb_file: str) -> None:
    """
    This function creates a new pdb because the CRYST line is missing from the pdf file.
    This line is necessary for DSSP.

    Args:
        output: str
        pdb_file: str

    Returns:
        pdb_file_copy: str

    """
    with open(pdb_file, "r") as file_in,\
    tempfile.NamedTemporaryFile(mode = "w", delete=False) as file_out:
        content = file_in.read()
        file_out.write(CRYST_LINE + content)

    return file_out.name


def determine_ss(structure, skipss, pdbf_path):
    """

    Args:
        structure:
        skipss:
        pdbf_path:

    Returns:
        structure:

    """
    # calculate SS
    for model in structure:

        if skipss:
            continue
        try:
            tmp_file_name = create_file_with_cryst(pdbf_path)
            p = subprocess.Popen(["dssp", tmp_file_name, "--output-format", "dssp"],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 text=True)
            dssp_raw, _ = p.communicate()
            dssp_raw = dssp_raw.split('#')[1].split('\n')[1:-1]
        except:  # TODO: think about making this exception more specific
            # no secondary structure detected for this model
            log.warning('SS could not be assigned, assigning code 1 to all residues')
            continue
        finally:
            if Path(tmp_file_name).exists(): Path(tmp_file_name).unlink()

        dssp = {}

        for line in dssp_raw:
            var_a, var_b = line[11], line[16]
            dssp.setdefault(var_a, []).append(var_b)

        # Get SS information and translate it:
        # DSSP > MARTINI > HADDOCK
        # this could still be improved
        for chain in model:
            dssp_ss = "".join(dssp[chain.id])
            _, martini_types = ss_classification(dssp_ss)  # ancestral function, keep it there

            # transform MARTINI > HADDOCK
            # and add it to the bfactor col
            for residue, ss in zip(chain, martini_types): # chain and martini_types order must match
                code = ss_to_code[ss]
                # for atom in residue.get_atoms():
                for atom in residue:
                    atom.bfactor = code
    return structure


def rename_nucbases(structure):
    """

    Args:
        structure:

    Returns:

    """
    chainresdic = dict([(c.get_id(),
                       [r.get_resname() for r in c.get_residues()]) for m in structure for c in m])

    nucleotide_list = ["CYT", "C", "DC", "THY", "T", "DT", "ADE",
                       "A", "DA", "G", "GUA", "DG", "U", "URI"]

    if [True for c in chainresdic for e in chainresdic[c] if e in nucleotide_list]:

        if [True for c in chainresdic for e in chainresdic[c] if e in ["U", "URI"]]:
            # CG needs 1 letter for RNA
            ref_dic = {"CYT": "C", "URI": "U", "ADE": "A", "GUA": "G"}
        else:
            # CG needs 2 letters for DNA
            ref_dic = {"CYT": "DC", "THY": "DT", "ADE": "DA", "GUA": "DG"}

        for model in structure:
            for chain in model:
                for r in chain.get_residues():
                    if r.resname in ref_dic.keys():
                        # rename
                        r.resname = ref_dic[r.resname]


def martinize(input_pdb, output_path, skipss):
    """
    
    Args:
        input_pdb: str
        output_path: str
        skipss: boolean
    
    Returns:
        cg_pdb_name: str
    """

    if not input_pdb:
        emsg = "No input file detected"
        raise ModuleError(emsg)

    p = PDBParser()
    io = PDBIO()

    # Parse PDB and run DSSP
    pdbf_path = os.path.realpath(input_pdb)
    aa_model = p.get_structure("aa_model", pdbf_path)

    # set ALL bfactors to 1
    for model in aa_model:
        for chain in model:
            if chain.id == " ":
                emsg = "Empty chain id detected"
                raise ModuleError(emsg)

            for residue in chain:
                for atom in residue:
                    atom.bfactor = 1.0

    # Assign HADDOCK code according to SS (1-9)
    determine_ss(structure=aa_model, skipss=skipss, pdbf_path=pdbf_path)

    # Strandardize naming
    # WARNING, THIS ASSUMES THAT INPUT DNA/RNA IS 3-LETTER CODE
    rename_nucbases(aa_model)

    # Assign HADDOCK code for hydrogen bonding capable nucleotides (0-1)
    pair_list = determine_hbonds(aa_model)
    if pair_list:
        output_cg_restraints(pair_list)

    # Map CG beads to AA structure
    structure_builder = StructureBuilder()
    structure_builder.init_structure("cg_model")
    structure_builder.init_seg(" ")  # Empty SEGID

    tbl_cg_to_aa = []
    restrain_counter = 0
    for model in aa_model:

        structure_builder.init_model(model.id)

        for chain in model:

            structure_builder.init_chain(chain.id)
            structure_builder.init_seg(chain.id)

            mapping_dic = map_cg(chain)

            for residue in mapping_dic:
                if residue.id[0] != " ":  # filter HETATMS
                    continue

                structure_builder.init_residue(residue.resname, residue.id[0],
                                               residue.id[1], residue.id[2])

                for i, bead in enumerate(mapping_dic[residue]):

                    bead_name = bead
                    bead_coord = mapping_dic[residue][bead_name][0]
                    haddock_code = mapping_dic[residue][bead_name][1]
                    restrain = mapping_dic[residue][bead_name][2]

                    structure_builder.init_atom(
                        bead_name,
                        bead_coord,
                        haddock_code,
                        1.00,
                        " ",
                        bead_name,
                        i)

                    tbl_cg_to_aa.append(restrain)
                    restrain_counter += 1

    cg_model = structure_builder.get_structure()

    # Write CG structure
    cg_pdb_name = f"../{output_path}/{pdbf_path.split('/')[-1][:-4]}_cg.pdb"
    io.set_structure(cg_model)
    io.save("temp.pdb", write_end=1)

    # make sure atom names are in the correct place
    # .BB. .BB1. .BB2. and not BB.. BB1.. BB2..
    out = open(cg_pdb_name, "w")
    for line in open("temp.pdb", "r"):
        if "ATOM" in line[:4]:
            atom_name = line[12:16].split()[0]
            # mind the spacing
            if len(atom_name) == 3:
                n_l = f"{line[:12]} {atom_name}{line[16:]}"
            elif len(atom_name) == 2:
                n_l = f"{line[:12]} {atom_name} {line[16:]}"
            elif len(atom_name) == 1:
                n_l = f"{line[:12]} {atom_name}  {line[16:]}"
            else:
                n_l = line
        else:
            n_l = line
        out.write(n_l)
    out.close()
    Path("temp.pdb").unlink(missing_ok=True)

    # Write Restraints
    tbl_file_name = f"../{output_path}/{pdbf_path.split('/')[-1][:-4]}_cg_to_aa.tbl"
    tbl_file = open(tbl_file_name, "w")
    tbl_str = "\n".join([tbl for tbl in tbl_cg_to_aa if tbl])
    tbl_file.write(f"\n{tbl_str}")
    tbl_file.close()

    return cg_pdb_name

