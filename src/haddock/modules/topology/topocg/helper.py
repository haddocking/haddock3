import collections
import math

#  CODE TAKEN FROM MARTINIZE 1.1 ##
"""
Reduces complexity of protein residue to the MARTINI coarse grained model:
CA, O, Bead(s) in specific atom location.

Reference:
Monticelli et al. The MARTINI coarse-grained force field: extension to proteins. 
J. Chem. Theory Comput. (2008) vol. 4 (5) pp. 819-834

Martinize Script from Tserk Wassenaar
"""

# Please refer to it for more information
# SECONDARY STRUCTURE DEFINITION


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
