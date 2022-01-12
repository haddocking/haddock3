"""
HADDOCK 3 supported residues.

https://bianca.science.uu.nl/haddock2.4/library


Write decent documentation .-)
"""
import itertools as it
import re
import string
from copy import copy
from pathlib import Path

from haddock import toppar_path


class CNSTopologyResidue:

    def __init__(self, resname, charge, atoms):
        """
        CNS topology residue.

        This is not an actual topology. Is a namespace storing the
        basic parameters defined in CNS residue topologies.

        So far, this is used by the `preprocessing` step (see
        `gear.preprocessing`; but this class can be extended with
        additional parameters if needed.

        Parameters
        ----------
        resname : str
            The name of the residue.

        charge : float
            The charge of the residue. The charge will be round to two
            decimal places.

        atoms : list-like
            The list of atom names. Will be converted to a tuple.
        """
        self._resname = resname
        self._charge = round(charge, 2)
        self._atoms = tuple(atoms)
        return

    def __repr__(self):
        return str(self)

    def __str__(self):
        return f"{self.__class__.__name__}({self.resname!r}, {self.charge!r}, {self.atoms!r})"

    @property
    def resname(self):
        return self._resname

    @property
    def charge(self):
        return self._charge

    @property
    def atoms(self):
        return self._atoms

    @property
    def elements(self):
        try:
            return self._elements
        except AttributeError as err:
            emsg = "`elements` is not defined. Use `make_elements` first."
            raise AttributeError(emsg) from err

    def make_elements(self, n=1):
        """Generate elements from atoms."""
        ele = [a.rstrip(string.digits + "+-")[:n] for a in self.atoms]
        self._elements = tuple(ele)
        return


def read_residues_from_top(topfile):
    """
    Read the residues defined in CNS topology file.

    Parameters
    ----------
    topfile : str or pathlib.Path
        Path to the topology file.

    Returns
    -------
    list
        List of :class:`CNSTopologyResidue` objects.
        Each object has information on the residue name, the atoms,
        and the total charge of the residue.
    """
    # regular expression to parse information from atom lines
    atom_regex = (
        r"^ATOM +([A-Z0-9\'\+\-]{1,4}) +.* +(?:charge|CHARGE|CHARge) *= *"
        r"([\-|+]?\d+.?\d*).*(?:END|end).*$"
        )

    # helper variables for the for-loop
    atoms, charges, residues = [], [], []
    in_residue = False

    # just to avoid identention with with
    fin = open(topfile, 'r')

    # all lines are striped before looping
    lines = map(str.strip, fin)

    for line in lines:
        if line.startswith(('RESI', 'residue')):
            in_residue = True  # defines we start a new residue

            # extracts the resname from the residue line
            resname = line.split()[1]

            # ensures there is no bug when changing residues
            # before starting a new residue `atoms` and `charges` should be
            # empty
            if atoms or charges:
                emsg = \
                    "Starting a new residue without emptying the previous one."
                raise ValueError(emsg)

        elif line.startswith('ATOM'):
            aname, acharge = re.findall(atom_regex, line)[0]
            atoms.append(aname)
            charges.append(acharge)

        elif line.startswith(("end", "END")) and in_residue:

            # the total charge of the residue
            total_charge = sum(map(float, charges))
            # need to be a copy to avoid unexpected referencing from the list
            _atoms = copy(atoms)

            resi = CNSTopologyResidue(resname, total_charge, copy(atoms))
            residues.append(resi)

            # clear help variables
            charges.clear()
            atoms.clear()
            in_residue = False

    fin.close()
    return residues


# what we will need to inspect if the residues are allowed or not are
# the residue names. So it is nice to have a function that retrieve them
# from the respective named tuples.
def get_resnames(group):
    """
    Make a set of residue names.

    Elements of the group must have a `resname` attribute.
    """
    return {i.resname: i for i in group}


# this function is necessary because of the `self_contained` option
# where the topology files are not defined in the haddock3 code but in
# the run folder
def read_supported_residues(source_path):
    # paths to the `.top` files
    carbo_top = Path(source_path, "carbohydrate.top")
    dna_rna_all_top = Path(source_path, "dna-rna-allatom-hj-opls-1.3.top")
    dna_rna_martini_top = Path(source_path, "dna-rna-CG-MARTINI-2-1p.top")
    fragment_top = Path(source_path, "fragment_probes.top")
    hemes_top = Path(source_path, "hemes-allhdg.top")
    ions_top = Path(source_path, "ion.top")
    protein_5_4_top = Path(source_path, "protein-allhdg5-4.top")
    protein_martini_2_top = Path(source_path, "protein-CG-Martini-2-2.top")
    protein_martini_top = Path(source_path, "protein-CG-Martini.top")
    solvent_top = Path(source_path, "solvent-allhdg5-4.top")


    # supported Residues (tuple of namedtuples)
    supported_carbohydrates = read_residues_from_top(carbo_top)

    _1 = read_residues_from_top(dna_rna_all_top)
    _2 = read_residues_from_top(dna_rna_martini_top)
    supported_nucleic = set(it.chain(_1, _2))

    supported_fragments = read_residues_from_top(fragment_top)

    supported_hemes = read_residues_from_top(hemes_top)

    # all ions are defined here
    supported_ions = read_residues_from_top(ions_top)

    # separates both kinds of ions
    supported_single_ions = \
        tuple(ion for ion in supported_ions if len(ion.atoms) == 1)
    supported_multiatom_ions = \
        tuple(ion for ion in supported_ions if len(ion.atoms) > 1)

    _l0 = len(supported_ions)
    _l1 = len(supported_single_ions)
    _l2 = len(supported_multiatom_ions)
    if _l0 == _l1 + _l2:
        emsg = (
            "Ions were not parsed correctly. Some atoms were left out "
            "when parsing single atom and multiatom ions."
            )
        raise ValueError(emsg)

    _1 = read_residues_from_top(protein_5_4_top)
    _2 = read_residues_from_top(protein_martini_2_top)
    _3 = read_residues_from_top(protein_martini_top)
    supported_aminoacids = set(it.chain(_1, _2, _3))

    supported_solvents = read_residues_from_top(solvent_top)

    # supported resnames
    supported_carbo_resnames = get_resnames(supported_carbohydrates)
    supported_nucleic_resnames = get_resnames(supported_nucleic)
    supported_fragments_resnames = get_resnames(supported_fragments)
    supported_hemes_resnames = get_resnames(supported_hemes)
    supported_single_ions_resnames = get_resnames(supported_single_ions)
    supported_multiatom_ions_resnames = get_resnames(supported_multiatom_ions)
    supported_aminoacids_resnames = get_resnames(supported_aminoacids)
    supported_solvents_resnames = get_resnames(supported_solvents)

    # other attributes
    for ion in supported_single_ions:
        ion.make_elements(n=2)

    supported_ions_elements = \
        {ion.elements[0]: ion for ion in supported_single_ions}

    supported_ions_atoms = {ion.atoms: ion for ion in supported_ions}

    return (
        supported_carbo_resnames,
        supported_nucleic_resnames,
        supported_fragments_resnames,
        supported_hemes_resnames,
        supported_ions_resnames,
        supported_multiatom_ions_resnames,
        supported_aminoacids_resnames,
        supported_solvents_resnames,
        )


supported_carbo_resnames, \
supported_nucleic_resnames, \
supported_fragments_resnames, \
supported_hemes_resnames, \
supported_ions_resnames, \
supported_multiatom_ions_resnames, \
supported_aminoacids_resnames, \
supported_solvents_resnames = read_supported_residues(toppar_path)

#
# Residues that must be set as ATOM
supported_ATOM = set(it.chain(
    supported_nucleic_resnames,
    supported_aminoacids_resnames,
    ))

# Residues that must be set as HETATM
supported_HETATM = set(it.chain(
    supported_carbo_resnames,
    supported_fragments_resnames,
    supported_hemes_resnames,
    supported_ions_resnames,
    supported_multiatom_ions_resnames,
    supported_solvents_resnames,
    ))
