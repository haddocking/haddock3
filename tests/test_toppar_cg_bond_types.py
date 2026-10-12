"""Test the bead types of the coarse-grained (Martini) nucleotides.

CNS looks bond parameters up by the pair of atom types, so two bonds of a
nucleotide that join the same pair of types are forced to the same length.
Adenine SC1 and SC2 used to share a type, which stretched the 1.62 A
SC1-SC4 bond to the 3.3 A of SC2-SC4.
"""

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
