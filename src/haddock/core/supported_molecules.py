"""
HADDOCK 3 supported residues.

https://bianca.science.uu.nl/haddock2.4/library
"""
import re
import itertools as it
from pathlib import Path
from collections import namedtuple

from haddock import toppar_path


# must be defined as HETATM
supported_carbohydrates = (
    'A2G',   # 2-N-acetyl-alpha-D-glucopyranose, different stereochemistry at C4
    'BGC',  # beta-D-glucopyranose
    'BMA',  # beta-D-mannopyronose
    'FCA',  # alpha-D-fucopyranose
    'FCB',  # beta-D-fucopyranose
    'FUC',  # alpha-L-fucopyranose
    'FUL',  # beta-L-fucopyranose
    'GAL',  # (GLB): beta-D-galactopyranose
    'GLA',  # alpha-D-galactopyranose
    'GLC',  # alpha-D-glucopyranose
    'MAN',  # alpha-D-mannopyranose
    'NAG',  # 2-N-acetyl-beta-D-glucopyranose
    'NDG',  # 2-N-acetyl-alpha-D-glucopyranose
    'NGA',  # 2-N-acetyl-beta-D-galactopyranose
    'SIA',  # alpha-N-acetyl neuraminic acid
    'SIB',  # beta-N-acetyl neuraminic acid
    'XYP',  # beta-D-xylopyranose
    )

# must be defined as HETATM
#ion_charges = {
#    'AG': "+1",  # Silver
#    'AL': "+3",  # Aluminium
#    'AU': "+1",  # Gold
#    'BR': "-1",  # Bromine
#    'CA': "+2",  # Calcium
#    'CD': "+2",  # Cadmium
#    'CL': "-1",  # Chlore
#    'CO': "+2",  # Cobalt
#    'CR': "+3",  # Chromium
#    'CS': "+1",  # Cesium
#    'CU': "+1",  # Copper
#    'F': "-1",   # Fluor
#    'FE': "+2",  # Iron
#    'HG': "+2",  # Mercury
#    'HO': "+3",  # Holmium
#    'I': "-1",   # Iodine
#    'IR': "+3",  # Iridium
#    'K': "+1",   # Potassium
#    'KR': "+2",  # Krypton
#    'LI': "+1",  # Lithium
#    'MG': "+2",  # Magnesium
#    'MN': "+2",  # Manganese
#    'MO': "+2",  # Molybdenum
#    'NA': "+1",  # Sodium
#    'NI': "+2",  # Nickel
#    'OS': "+6",  # Osmium
#    'PB': "+2",  # Lead
#    'PT': "+2",  # Platinum
#    'SR': "+2",  # Strontium
#    'U': "+2",  # Uranium
#    'V': "+2",  # Vanadium
#    'YB': "+3",  # Ytterbium
#    'ZN': "+2",  # Zinc
#    }


def read_ion_topology(topfile):
    """
    Read ion topology from topfile.

    Expects a file with the format of `toppar/ion.top`.

    See regex:
        https://regex101.com/r/T9ak8y/1
    """
    Ion = namedtuple('Ion', ['resname', 'name', 'full', 'atom', 'charge'])
    regex = r'\nRESIdue (\w{1,4}) {(\w+ \d[\+\-])}\n.*\n *ATOM (\w{1,2}[\+|\-]?\d{1,2}) .*\nEND {\w*}\n'  # noqa: E501
    text = Path(topfile).read_text()
    groups = re.findall(regex, text)

    ions = tuple(
        Ion(
            resname=group[0],
            name=group[2][:-2],
            full=group[1],
            atom=group[2],
            charge=group[2][-2:],
            )
        for group in groups
        )

    return ions


supported_ions = read_ion_topology(Path(toppar_path, 'ion.top'))
supported_ion_names = tuple(i.name for i in supported_ions)
supported_ion_resnames = tuple(i.resname for i in supported_ions)


def read_multiatom_topology(topfile):
    """
    Read multi-atom ions topology from topfile.

    Expects a file with the format of `toppar/ion.top`.

    See regex:
        https://regex101.com/r/nkPMUp/1
    """
    Ion = namedtuple('Ion', ['resname', 'name'])
    regex = r'RESIdue (\w{1,4}) {(\w+)}\n  GROUP\n(?:    ATOM.*\n)+\n(?:  BOND.*\n)+'  # noqa: E501
    text = Path(topfile).read_text()
    groups = re.findall(regex, text)
    ions = tuple(Ion(resname=group[0], name=group[1]) for group in groups)
    return ions



# must be defined as HETATM
supported_multiatom_ions = read_multiatom_topology(Path(toppar_path, 'ion.top'))
supported_multiatom_ions_resnames = tuple(i.resname for i in supported_multiatom_ions)

# must be defined as HETATM
supported_cofactors = (
    'HEB',  # Heme-B
    'HEC',  # Heme-C
    )

# must be defined as ATOM
supported_nucleic_acid_bases = (
    # DNA
    'DA',  # Adenine
    'DC',  # Cytosine
    'DG',  # Guanidine
    'DT',  # Thymidine

    # RNA
    'A',  # Adenine
    'C',  # Cytosine
    'G',  # Guanidine
    'U',  # Uridine
    )

# must be defined as ATOM
supported_modified_amino_acids = (
    'ACE',  # N-terminal acetyl group
    'ALY',  # Acetylated LYS
    'ASH',  # protonated ASP
    'CFE',  # CYS with an iron sulfur cluster
    'CSP',  # phosphorylated CYS
    'CTN',  # C-terminal amide group
    'CYC',  # special for covalent docking - vdw of Sulphur reduced
    'CYF',  # CYS without the sulfur H (for coordinating metals)
    'CYM',  # CYS with MTSL grouped
    'DDZ',  # 3,3,-dihydrozy ALA
    'GLH',  # protonated GLU
    'HYP',  # 4R-hydroxyproline
    'M3L',  # trimethyl LYS
    'MLY',  # dimethyl LYS
    'MLZ',  # monomethyl LYS
    'MSE',  # Selenomethionine
    'NEP',  # NE phosphorylated HIS
    'NME',  # C-terminal N-Methyl
    'PNS',  # Phosphopanthenite Serine
    'PTR',  # O-Phosphotyrosine
    'SEP',  # phosphorylated SER
    'TOP',  # phosphorylated THR
    'TYP',  # phosphorylated TYR
    'PTR',  # phosphorylated TYR (also)
    'TYS',  # sulfonated TYR.
    )

supported_natural_amino_acids = (
    'ALA',
    'ARG',
    'ASN',
    'ASP',
    'CYS',
    'GLU',
    'GLN',
    'GLY',
    'HIS',
    'HIP',
    'HIE',
    'HID',
    'ILE',
    'LEU',
    'LYS',
    'MET',
    'PHE',
    'PRO',
    'SER',
    'THR',
    'TRP',
    'TYR',
    'VAL',
    )

supported_atom = set(it.chain(
    supported_natural_amino_acids,
    supported_modified_amino_acids,
    supported_nucleic_acid_bases,
    ))

supported_hetatm = set(it.chain(
    supported_cofactors,
    supported_multiatom_ions,
    supported_ions,
    supported_ion_resnames,
    supported_carbohydrates,
    ))
