"""
HADDOCK 3 supported residues.

https://bianca.science.uu.nl/haddock2.4/library
"""

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
supported_ions = (
    'AG',  # Silver
    'AL',  # Aluminium
    'AU',  # Gold
    'BR',  # Bromine
    'CA',  # Calcium
    'CD',  # Cadmium
    'CL',  # Chlore
    'CO',  # Cobalt
    'CR',  # Chromium
    'CS',  # Cesium
    'CU',  # Copper
    'F',   # Fluor
    'FE',  # Iron
    'HG',  # Mercury
    'HO',  # Holmium
    'I',   # Iodine
    'IR',  # Iridium
    'K',   # Potassium
    'KR',  # Krypton
    'LI',  # Adenine
    'MG',  # Magnesium
    'MN',  # Manganese
    'MO',  # Molybdenum
    'NA',  # Sodium
    'NI',  # Nickel
    'OS',  # Osmium
    'PB',  # Lead
    'PT',  # Platinum
    'SR',  # Strontium
    'U',  # Uranium
    'V',  # Vanadium
    'YB',  # Ytterbium
    'ZN',  # Zinc
    )

# must be defined as HETATM
supported_multiatom_ions = (
    'PO4',  # Phosphate
    'SO4',  # Sulphate
    'WO4',  # Tungstate
    )

supported_cofactors = (
    'HEB',  # Heme-B
    'HEC',  # Heme-C
    )

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
