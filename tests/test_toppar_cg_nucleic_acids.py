"""Test the coarse-grained (Martini) nucleic acid topology and parameters.

Every bond, angle and dihedral that the CG topology defines within a
nucleotide must have a parameter, also after the hydrogen-bond patches
(DA1, A1, ...) have retyped the base beads. A missing one makes CNS abort
with %CODBON-ERR when topocg builds the topology.
"""

import re
from pathlib import Path

import pytest

from haddock import toppar_path


TOPOLOGY = Path(toppar_path, "dna-rna-CG-MARTINI-2-1p.top")
PARAMETERS = Path(toppar_path, "dna-rna-CG-MARTINI-2-1p.param")
NUCLEOTIDES = ("DA", "DC", "DG", "DT", "A", "C", "G", "U")
N_ATOMS = {"BOND": 2, "ANGL": 3, "DIHE": 4, "IMPR": 4}


def read_parameters() -> dict[str, set[tuple[str, ...]]]:
    """Return the atom-type tuples parameterised for each term."""
    params: dict[str, set[tuple[str, ...]]] = {term: set() for term in N_ATOMS}
    for line in PARAMETERS.read_text().splitlines():
        fields = line.split("!")[0].split()
        term = fields[0][:4].upper() if fields else ""
        if term in N_ATOMS:
            types = tuple(t.upper() for t in fields[1 : 1 + N_ATOMS[term]])
            params[term].add(types)
    return params


def read_residues() -> dict[str, tuple[dict[str, str], list]]:
    """Return atom types and intra-residue bonded terms of each residue."""
    residues = {}
    text = TOPOLOGY.read_text()
    for match in re.finditer(r"RESIdue (\S+)(.*?)\nEND", text, re.S):
        name, body = match.groups()
        types = dict(re.findall(r"ATOM (\S+)\s+type=(\S+)", body, re.I))
        terms = []
        for term, n_atoms in N_ATOMS.items():
            pattern = rf"^\s*{term}\w*((?:\s+\S+){{{n_atoms}}})"
            for atoms in re.findall(pattern, body, re.M | re.I):
                terms.append((term, atoms.split()))
        residues[name] = (types, terms)
    return residues


def read_patches() -> dict[str, dict[str, str]]:
    """Return the atom retyping done by each patch."""
    text = TOPOLOGY.read_text()
    return {
        match.group(1): dict(
            re.findall(r"MODIFY ATOM (\S+)\s+TYPE=(\S+)", match.group(2), re.I)
        )
        for match in re.finditer(r"PRESidue (\S+)(.*?)\nEND", text, re.S)
    }


@pytest.mark.parametrize("resname", NUCLEOTIDES)
@pytest.mark.parametrize("patched", [False, True])
def test_cg_nucleotide_terms_parameterised(resname, patched):
    """All intra-residue bonded terms of a CG nucleotide have parameters."""
    params = read_parameters()
    types, terms = read_residues()[resname]
    if patched:
        types = {**types, **read_patches()[f"{resname}1"]}

    missing = []
    for term, atoms in terms:
        if any(atom.startswith(("+", "-")) for atom in atoms):
            continue  # inter-residue terms belong to the linkage
        atom_types = tuple(types[atom].upper() for atom in atoms)
        if atom_types not in params[term] and atom_types[::-1] not in params[term]:
            missing.append(f"{term} {' '.join(atoms)} ({' '.join(atom_types)})")
    assert not missing, missing
