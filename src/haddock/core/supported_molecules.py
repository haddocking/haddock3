"""
HADDOCK 3 supported residues.

The HADDOCK3 supported residues are defined by the different ``.top``
files in the ``cns/toppar`` source folder. This module will automatically
retrieve any new residue added to the ``.top`` defined in the folder and
add it to the variables here.

If a new ``.top`` file is added in the folder, that must be added here
inside the :py:func:`read_supported_residues` function, because file retrieve is
not done automatically on purpose.

``.top`` files considered here:

* ``carbohydrate.top``
* ``dna-rna-CG-MARTINI-2-1p.top``
* ``dna-rna-allatom-hj-opls-1.3.top``
* ``fragment_probes.top``
* ``hemes-allhdg.top``
* ``ion.top``
* ``protein-CG-Martini-2-2.top``
* ``protein-CG-Martini.top``
* ``protein-allhdg5-4.top``
* ``solvent-allhdg5-4.top``
* ``shape.top``
* ``cofactors.top``

DNA-related files are all concatenated into a single DNA-supported
residues datastructure. The same occurs for protein residues (natural or
modified).

You can import from this module the different sets containing the
supported residue names:

* :py:data:`supported_aminoacids_resnames`
* :py:data:`supported_carbo_resnames`
* :py:data:`supported_fragments_resnames`
* :py:data:`supported_hemes_resnames`
* :py:data:`supported_single_ions_resnames`
* :py:data:`supported_single_ions_elements`
* :py:data:`supported_single_ions_atoms`
* :py:data:`supported_multiatom_ions_resnames`
* :py:data:`supported_nucleic_resnames`
* :py:data:`supported_solvents_resnames`
* :py:data:`supported_shape`
* :py:data:`supported_cofactors_resnames`

* :py:data:`supported_ATOM` (contains residues and nucleic acids)
* :py:data:`supported_HETATM` (everything not contained in ATOM)
"""
import itertools as it
import re
import string
from copy import copy
from pathlib import Path

from haddock import toppar_path
from haddock.core.typing import FilePath, Iterable


class CNSTopologyResidue:
    """CNS topology residue."""

    def __init__(self, resname: str, charge: float,
                 atoms: Iterable[str]) -> None:
        """
        CNS topology residue.

        This is not an actual topology. It is a namespace storing the
        basic parameters defined in CNS residue topologies.

        So far, this is used by the
        :py:module:`haddock.gear.preprocessing` step; but this class can
        be extended with additional parameters if needed.

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

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return (
            f"{self.__class__.__name__}({self.resname!r}, "
            f"{self.charge!r}, {self.atoms!r})"
            )

    @property
    def resname(self) -> str:
        """Residue name."""
        return self._resname

    @property
    def charge(self) -> float:
        """Residue total charge."""
        return self._charge

    @property
    def atoms(self) -> tuple[str, ...]:
        """Tuple of residue atoms."""
        return self._atoms

    @property
    def elements(self) -> tuple[str, ...]:
        """
        List of elements for all atoms composing the molecule.

        This list is ordered the same as ``atoms``.

        By default ``elements`` is not defined and must be injected manually.

        See Also
        --------
        * :py:func:`get_ion_element_from_atoms`
        """
        try:
            return self._elements
        except AttributeError as err:
            emsg = "`elements` is not defined."
            raise AttributeError(emsg) from err

    @elements.setter
    def elements(self, elements: tuple[str, ...]) -> None:
        self._elements = elements


def get_ion_element_from_atoms(atoms: Iterable[str]) -> tuple[str, ...]:
    """
    Generate the element name from the atom name for single atom ions.

    Examples
    --------
    >>> get_ion_element_from_atoms(["ZN+2", "CU+2", "NI"])
    ("ZN", "CU", "NI")

    Parameters
    ----------
    atom : list
        A list of atom names.
    """
    ele = (a.rstrip(string.digits + "+-") for a in atoms)
    return tuple(ele)


def read_residues_from_top_file(topfile: FilePath) -> list[CNSTopologyResidue]:
    """
    Read the residues defined in CNS topology file.

    This is a general implementation to read the basic information from
    the CNS topology files. This basic information is the information
    needed to help the the :py:mod:`haddock.gear.preprocessing` step
    in identifying the supported residues and correctly format their PDB
    lines.

    Parameters
    ----------
    topfile : str or :external:py:class:`pathlib.Path`
        Path to the topology file.

    Returns
    -------
    list
        List of :class:`CNSTopologyResidue` objects.
        Each object has information on the residue name, the atoms,
        and the total charge of the residue.

    See Also
    --------
    * :py:func:`haddock.gear.preprocessing.read_additional_residues`
    """
    # just to avoid identention with with
    fin = open(topfile, 'r')
    lines = fin.readlines()
    fin.close()
    return _read_residues_from_top_file(lines)


def _read_residues_from_top_file(
        lines: Iterable[str]) -> list[CNSTopologyResidue]:
    # regular expression to parse information from atom lines
    atom_regex = (
        r"^ATOM +([A-Z0-9\'\+\-]{1,4}) +.* +(?:charge|CHARGE|CHARge) *= *"
        r"([\-|+]?\d+.?\d*).*(?:END|end).*$"
        )

    # helper variables for the for-loop
    atoms: list[str] = []
    charges: list[str] = []
    residues: list[CNSTopologyResidue] = []
    in_residue: bool = False

    # all lines are striped before looping
    for line in map(str.strip, lines):  # it is important to strip spaces
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

            resi = CNSTopologyResidue(resname, total_charge, _atoms)
            residues.append(resi)

            # clear help variables
            charges.clear()
            atoms.clear()
            in_residue = False

    return residues


# what we will need to inspect if the residues are allowed or not are
# the residue names. So it is nice to have a function that retrieve them
# from the respective named tuples.
def get_resnames(
        group: Iterable[CNSTopologyResidue]) -> dict[str, CNSTopologyResidue]:
    """
    Make a set of residue names.

    Elements of the group must have a `resname` attribute.

    Parameters
    ----------
    group : list of :py:class:`CNSTopologyResidue`

    Returns
    -------
    dict
        Dict of residue name and residue pairs.
    """
    return {i.resname: i for i in group}


# this function is necessary because of the `self_contained` option
# where the topology files are not defined in the haddock3 code but in
# the run folder
def read_supported_residues(
        source_path: FilePath) -> tuple[dict[str, CNSTopologyResidue], ...]:
    """
    Read supported residues from a folder containing ``.top`` files.

    The function expects the ``source_path`` to have the required files.

    Expected files are:

    * ``carbohydrate.top``
    * ``dna-rna-CG-MARTINI-2-1p.top``
    * ``dna-rna-allatom-hj-opls-1.3.top``
    * ``fragment_probes.top``
    * ``hemes-allhdg.top``
    * ``ion.top``
    * ``protein-CG-Martini-2-2.top``
    * ``protein-CG-Martini.top``
    * ``protein-allhdg5-4.top``
    * ``solvent-allhdg5-4.top``
    * ``shape.top``
    * ``cofactors.top``

    You need to edit this function to account for any additional
    ``.top`` file that is added to the ``source_path`` on top of the
    above ones. Otherwise, read the residues manually using
    :py:func:`read_residues_from_top_file` and :py:func:`get_resnames`
    and add them to the relevant set, or have them as a new set.

    Returns
    -------
    tuple
        A tuple of the residue names supported by HADDOCK3.
    """
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
    shape_top = Path(source_path, "shape.top")
    cofactors_top = Path(source_path, "cofactors.top")

    # supported Residues (tuple of namedtuples)
    supported_carbohydrates = read_residues_from_top_file(carbo_top)

    _1 = read_residues_from_top_file(dna_rna_all_top)
    _2 = read_residues_from_top_file(dna_rna_martini_top)
    supported_nucleic = set(it.chain(_1, _2))

    supported_fragments = read_residues_from_top_file(fragment_top)

    supported_hemes = read_residues_from_top_file(hemes_top)

    # all ions are defined here
    supported_ions = read_residues_from_top_file(ions_top)

    # separates both kinds of ions
    supported_single_ions = \
        tuple(ion for ion in supported_ions if len(ion.atoms) == 1)
    supported_multiatom_ions = \
        tuple(ion for ion in supported_ions if len(ion.atoms) > 1)

    _l0 = len(supported_ions)
    _l1 = len(supported_single_ions)
    _l2 = len(supported_multiatom_ions)
    if _l0 != _l1 + _l2:
        emsg = (
            "Ions were not parsed correctly. Some atoms were left out "
            "when parsing single atom and multiatom ions."
            )
        raise ValueError(emsg)

    _1 = read_residues_from_top_file(protein_5_4_top)
    _2 = read_residues_from_top_file(protein_martini_2_top)
    _3 = read_residues_from_top_file(protein_martini_top)
    supported_aminoacids = set(it.chain(_1, _2, _3))

    supported_solvents = read_residues_from_top_file(solvent_top)
    supported_shape = read_residues_from_top_file(shape_top)
    supported_cofactors = read_residues_from_top_file(cofactors_top)

    # supported resnames
    supported_carbo_resnames = get_resnames(supported_carbohydrates)
    supported_nucleic_resnames = get_resnames(supported_nucleic)
    supported_fragments_resnames = get_resnames(supported_fragments)
    supported_hemes_resnames = get_resnames(supported_hemes)
    supported_single_ions_resnames = get_resnames(supported_single_ions)
    supported_multiatom_ions_resnames = get_resnames(supported_multiatom_ions)
    supported_aminoacids_resnames = get_resnames(supported_aminoacids)
    supported_solvents_resnames = get_resnames(supported_solvents)
    supported_shape_resnames = get_resnames(supported_shape)
    supported_cofactors_resnames = get_resnames(supported_cofactors)

    # other attributes
    for ion in supported_single_ions:
        ion.elements = get_ion_element_from_atoms(ion.atoms)

    supported_ions_elements = \
        {ion.elements[0]: ion for ion in supported_single_ions}

    supported_ions_atoms = {ion.atoms[0]: ion for ion in supported_ions}

    return (
        supported_carbo_resnames,
        supported_nucleic_resnames,
        supported_fragments_resnames,
        supported_hemes_resnames,
        supported_single_ions_resnames,
        supported_ions_elements,
        supported_ions_atoms,
        supported_multiatom_ions_resnames,
        supported_aminoacids_resnames,
        supported_solvents_resnames,
        supported_shape_resnames,
        supported_cofactors_resnames,
        )


# reads and assigns the supported molecules variables to this module
_supported_carbo_resnames, \
    _supported_nucleic_resnames, \
    _supported_fragments_resnames, \
    _supported_hemes_resnames, \
    _supported_single_ions_resnames, \
    _supported_single_ions_elements, \
    _supported_single_ions_atoms, \
    _supported_multiatom_ions_resnames, \
    _supported_aminoacids_resnames, \
    _supported_solvents_resnames, \
    _supported_shape_resnames, \
    _supported_cofactors_resnames = read_supported_residues(toppar_path)

# render docstrings
supported_carbo_resnames = set(_supported_carbo_resnames)
"""Supported carbohydrate residues names."""

supported_nucleic_resnames = set(_supported_nucleic_resnames)
"""Supported nucleic acid residue names."""

supported_fragments_resnames = set(_supported_fragments_resnames)
"""Supported fragments residue names."""

supported_hemes_resnames = set(_supported_hemes_resnames)
"""Supported Hemes."""

supported_single_ions_resnames_map = _supported_single_ions_resnames
"""Supported single ion resname mapping."""

supported_single_ions_elements_map = _supported_single_ions_elements
"""Supported single ions elements mapping."""

supported_single_ions_atoms_map = _supported_single_ions_atoms
"""Supported single ion atom names mapping."""

supported_single_ions_resnames = set(_supported_single_ions_resnames)
"""Supported single residue names."""

supported_single_ions_elements = set(_supported_single_ions_elements)
"""Supported single ions elements."""

supported_single_ions_atoms = set(_supported_single_ions_atoms)
"""Supported single ion atom names."""

supported_multiatom_ions_resnames = set(_supported_multiatom_ions_resnames)
"""Supported multiatom ions residue names."""

supported_aminoacids_resnames = set(_supported_aminoacids_resnames)
"""Supported amino acids residue names."""

supported_solvents_resnames = set(_supported_solvents_resnames)
"""Supported solvents."""

supported_shape_resnames = set(_supported_shape_resnames)
"""Supported shape."""

supported_cofactors_resnames = set(_supported_cofactors_resnames)
"""Supported cofactors."""

#
# Residues that must be set as ATOM
supported_ATOM = set(it.chain(
    supported_nucleic_resnames,
    supported_aminoacids_resnames,
    supported_shape_resnames,
    ))
"""
Supported ``ATOM`` residues.

These residues must be defined as ``ATOM`` records
in PDB files for HADDOCK3.
"""

# non-ion residue names
supported_non_ions = set(it.chain(
    supported_ATOM,
    supported_carbo_resnames,
    supported_fragments_resnames,
    supported_hemes_resnames,
    supported_solvents_resnames,
    supported_shape_resnames,
    supported_cofactors_resnames,
    ))

# Residues that must be set as HETATM
supported_HETATM = set(it.chain(
    supported_carbo_resnames,
    supported_fragments_resnames,
    supported_hemes_resnames,
    supported_single_ions_resnames,
    supported_multiatom_ions_resnames,
    supported_solvents_resnames,
    supported_cofactors_resnames,
    ))
"""
Supported ``HETATM`` residues.

These residues must be defined as ``HETATM`` records
in PDB files for HADDOCK3.
"""

supported_residues = supported_ATOM.union(supported_HETATM)
"""All HADDOCK3 supported residues."""

