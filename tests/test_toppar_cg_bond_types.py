"""Test the bead types of the coarse-grained (Martini) nucleotides.

CNS looks bond parameters up by the pair of atom types, so two bonds of a
nucleotide that join the same pair of types are forced to the same length.
Adenine SC1 and SC2 used to share a type, which stretched the 1.62 A
SC1-SC4 bond to the 3.3 A of SC2-SC4.

The bonds of a base ring also fix its angles. Some of the original Martini
ring angles disagree with the ring bonds (adenine SC2-SC1-SC4 is 125 degrees
where the bonds give 111.8); these are replaced by the angles implied by the
bonds, so that the harmonic angle and bond terms of CNS do not compete.
"""

import math
import re
from pathlib import Path

import pytest

from haddock import toppar_path


TOPOLOGY = Path(toppar_path, "dna-rna-CG-MARTINI-2-1p.top")
PARAMETERS = Path(toppar_path, "dna-rna-CG-MARTINI-2-1p.param")
NUCLEOTIDES = ("DA", "DC", "DG", "DT", "A", "C", "G", "U")


def read_residue(text: str, keyword: str, name: str) -> str:
    """Return the body of a RESIdue or PRESidue entry."""
    match = re.search(rf"^{keyword} {name}\s*\n(.*?)\nEND", text, re.S | re.M)
    assert match, name
    return match.group(1)


def nucleotide_bonds(resname: str, patched: bool) -> tuple[list, dict]:
    """Return the intra-residue bonds and atom types of a CG nucleotide."""
    text = TOPOLOGY.read_text()
    body = read_residue(text, "RESIdue", resname)
    types = dict(re.findall(r"ATOM (\S+)\s+type=(\S+)", body, re.I))
    if patched:
        patch = read_residue(text, "PRESidue", f"{resname}1")
        types.update(re.findall(r"MODIFY ATOM (\S+)\s+TYPE=(\S+)", patch, re.I))
    bonds = re.findall(r"^\s*BOND (\S+) (\S+)", body, re.M)
    return bonds, types


def bond_length(type1: str, type2: str) -> float:
    """Return the equilibrium length of the bond between two atom types."""
    for line in PARAMETERS.read_text().splitlines():
        fields = line.split("!")[0].split()
        if fields[:1] == ["BOND"] and {*fields[1:3]} == {type1, type2}:
            return float(fields[4])
    raise KeyError((type1, type2))


def angle_value(type1: str, type2: str, type3: str) -> float:
    """Return the equilibrium value of the angle between three atom types."""
    for line in PARAMETERS.read_text().splitlines():
        fields = line.split("!")[0].split()
        if fields[:1] == ["ANGLE"] and fields[2] == type2:
            if (fields[1], fields[3]) in ((type1, type3), (type3, type1)):
                return float(fields[5])
    raise KeyError((type1, type2, type3))


def triangle_angle(side1: float, side2: float, opposite: float) -> float:
    """Return the angle between two sides of a triangle (law of cosines)."""
    cosine = (side1**2 + side2**2 - opposite**2) / (2 * side1 * side2)
    return math.degrees(math.acos(cosine))


def ring_angle_from_bonds(atoms: tuple, bonds: list, types: dict) -> float:
    """Return the angle of a base ring implied by its (planar) ring bonds.

    The angle is taken from the triangle it belongs to or, for the purine
    rings, from the two triangles that share the SC2-SC4 diagonal.
    """
    bonded = {frozenset(bond) for bond in bonds}

    def length(atom1, atom2):
        return bond_length(types[atom1], types[atom2])

    end1, apex, end2 = atoms
    if frozenset((end1, end2)) in bonded:
        return triangle_angle(
            length(apex, end1), length(apex, end2), length(end1, end2)
        )
    for middle in types:
        pairs = ((end1, middle), (apex, middle), (end2, middle))
        if middle not in atoms and all(frozenset(p) in bonded for p in pairs):
            return triangle_angle(
                length(apex, end1), length(apex, middle), length(end1, middle)
            ) + triangle_angle(
                length(apex, middle), length(apex, end2), length(middle, end2)
            )
    raise ValueError(atoms)


@pytest.mark.parametrize("resname", NUCLEOTIDES)
@pytest.mark.parametrize("patched", [False, True])
def test_cg_nucleotide_bonds_distinct(resname, patched):
    """Each bond of a CG nucleotide is listed once and has its own types."""
    bonds, types = nucleotide_bonds(resname, patched)
    atom_pairs = [frozenset(bond) for bond in bonds]
    assert len(set(atom_pairs)) == len(atom_pairs), bonds
    type_pairs = [frozenset(types[atom] for atom in bond) for bond in bonds]
    assert len(set(type_pairs)) == len(type_pairs), bonds


@pytest.mark.parametrize("resname", ["DA", "A"])
@pytest.mark.parametrize("patched", [False, True])
def test_cg_adenine_sc1_sc4_length(resname, patched):
    """Adenine SC1-SC4 follows the Martini 0.162 nm constraint."""
    _, types = nucleotide_bonds(resname, patched)
    assert bond_length(types["SC1"], types["SC4"]) == pytest.approx(1.62)


@pytest.mark.parametrize("resname", NUCLEOTIDES)
@pytest.mark.parametrize("patched", [False, True])
def test_cg_base_ring_angles_match_bonds(resname, patched):
    """Base ring angles agree with the ring bond lengths within 5 degrees."""
    bonds, types = nucleotide_bonds(resname, patched)
    body = read_residue(TOPOLOGY.read_text(), "RESIdue", resname)
    angles = re.findall(r"^\s*ANGLE (SC\d) (SC\d) (SC\d)", body, re.M)
    assert angles
    for atoms in angles:
        expected = ring_angle_from_bonds(atoms, bonds, types)
        observed = angle_value(*(types[atom] for atom in atoms))
        assert observed == pytest.approx(expected, abs=5.0), atoms
