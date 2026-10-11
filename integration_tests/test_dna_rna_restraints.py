"""Integration tests for the base pairing in dna-rna_restraints.cns.

The restraint script is run standalone in CNS on all-atom (topoaa) and
coarse-grained (topocg) topologies of:

- one duplex of 1SDR, an A-form RNA dodecamer (12 Watson-Crick pairs);
- 3GO3, B-DNA from an echinomycin complex, with 5 Watson-Crick and
  3 Hoogsteen pairs;
- 1D3X, a DNA triplex in which every purine of the central strand has a
  Watson-Crick and a Hoogsteen partner;
- 9TCG, an RNA triplex with 10 Watson-Crick and 2 Hoogsteen pairs.

The expected pairs were assigned from the hydrogen bonds of the all-atom
models, independently of the CNS script.
"""

import re
from collections import defaultdict
from pathlib import Path

import pytest

from haddock.libs.libsubprocess import CNSJob
from haddock.modules.refinement.emref import RECIPE_PATH as EMREF_PATH

from integration_tests import CNS_EXEC, GOLDEN_DATA


# Watson-Crick pairs as (A/DA or C/DC base, partner) and Hoogsteen pairs as
# (A/DA or G/DG base, partner), each as segid + resid, following the side
# the CNS script restrains them from.
WATSON_CRICK = {
    "rna_1SDR": {
        ("A2", "B23"),
        ("A3", "B22"),
        ("A6", "B19"),
        ("A11", "B14"),
        ("B13", "A12"),
        ("B15", "A10"),
        ("B16", "A9"),
        ("B17", "A8"),
        ("B18", "A7"),
        ("B20", "A5"),
        ("B21", "A4"),
        ("B24", "A1"),
    },
    "dna_3GO3": {
        ("A1", "B108"),
        ("A2", "B107"),
        ("A6", "B103"),
        ("B102", "A7"),
        ("B106", "A3"),
    },
    "dna_triplex_1D3X": {
        ("A1", "B16"),
        ("A3", "B14"),
        ("A5", "B12"),
        ("A7", "B10"),
        ("A8", "B9"),
        ("B11", "A6"),
        ("B13", "A4"),
        ("B15", "A2"),
    },
    "rna_triplex_9TCG": {
        ("A1", "B5"),
        ("A3", "B3"),
        ("A5", "B1"),
        ("B2", "A4"),
        ("B4", "A2"),
        ("C1", "D5"),
        ("C3", "D3"),
        ("C5", "D1"),
        ("D2", "C4"),
        ("D4", "C2"),
    },
}
HOOGSTEEN = {
    "rna_1SDR": set(),
    "dna_3GO3": {("A5", "B104"), ("B101", "A8"), ("B105", "A4")},
    "dna_triplex_1D3X": {
        ("A1", "C17"),
        ("A2", "C18"),
        ("A3", "C19"),
        ("A4", "C20"),
        ("A5", "C21"),
        ("A6", "C22"),
        ("A7", "C23"),
        ("A8", "C24"),
    },
    "rna_triplex_9TCG": {("A2", "D8"), ("C2", "B8")},
}

# restraints per pair, by representation, pair type and restrained base
RESTRAINTS_PER_PAIR = {
    "aa": {("WC", "ADE"): 4, ("WC", "CYT"): 5, ("HOOG", "ADE"): 4, ("HOOG", "GUA"): 4},
    "cg": {("WC", "ADE"): 2, ("WC", "CYT"): 3, ("HOOG", "ADE"): 3, ("HOOG", "GUA"): 3},
}
BASE = {"A": "ADE", "DA": "ADE", "C": "CYT", "DC": "CYT", "G": "GUA", "DG": "GUA"}

CNS_INPUT = """
set message=normal echo=off end
structure @GOLDEN:{name}.psf end
coor @GOLDEN:{name}.pdb
noe reset nrestraints=20000 ceiling=1000 end
restraints dihedral nassign=20000 end
@MODULE:dna-rna_restraints.cns
set echo=on end
noe print threshold=-1.0 end
stop
"""

NOE_RE = re.compile(
    r"set-i-atoms\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+"
    r"set-j-atoms\s+(\S+)\s+(\S+)\s+\S+\s+\S+"
)


def run_dna_rna_restraints(name: str) -> str:
    """Run dna-rna_restraints.cns on a golden topology, return the output."""
    job = CNSJob(
        input_file=CNS_INPUT.format(name=name),
        envvars={
            "MODULE": str(Path(EMREF_PATH, "cns")),
            "GOLDEN": str(GOLDEN_DATA),
        },
        cns_exec=CNS_EXEC,
    )
    return job.run().decode()


def restrained_pairs(output: str) -> dict:
    """
    Group the NOE restraints per (base, partner) pair.

    Returns {(pair type, base type): {(base, partner): n restraints}}. A
    Hoogsteen pair is recognised by its N7-edge atom on the purine side
    (N7 all-atom, SC4 CG), which no Watson-Crick restraint uses.
    """
    noe_table = output.split("noe print threshold")[-1]
    atoms = defaultdict(list)
    resnames = {}
    for seg1, res1, resn1, atom1, seg2, res2 in NOE_RE.findall(noe_table):
        pair = (f"{seg1}{res1}", f"{seg2}{res2}")
        atoms[pair].append(atom1)
        resnames[pair] = resn1
    pairs: dict = defaultdict(dict)
    for pair, pair_atoms in atoms.items():
        pair_type = "HOOG" if {"N7", "SC4"} & set(pair_atoms) else "WC"
        pairs[(pair_type, BASE[resnames[pair]])][pair] = len(pair_atoms)
    return pairs


@pytest.mark.parametrize("rep", ["aa", "cg"])
@pytest.mark.parametrize("structure", sorted(WATSON_CRICK))
def test_base_pairs(structure, rep):
    """Exactly the Watson-Crick and Hoogsteen pairs are restrained."""
    output = run_dna_rna_restraints(f"{structure}_{rep}")

    assert "ABORT mode" not in output
    cg_detected = "coarse-grained model detected" in output
    assert cg_detected == (rep == "cg")

    pairs = restrained_pairs(output)
    found = {
        pair_type: {
            pair
            for (ptype, _), group in pairs.items()
            if ptype == pair_type
            for pair in group
        }
        for pair_type in ("WC", "HOOG")
    }
    assert found["WC"] == WATSON_CRICK[structure]
    assert found["HOOG"] == HOOGSTEEN[structure]

    for key, group in pairs.items():
        for pair, n_restraints in group.items():
            assert n_restraints == RESTRAINTS_PER_PAIR[rep][key], (key, pair)


def test_dna_rna_restraints_copies_identical():
    """All modules must ship the same dna-rna_restraints.cns."""
    modules_path = EMREF_PATH.parents[1]
    copies = sorted(modules_path.glob("*/*/cns/dna-rna_restraints.cns"))
    assert len(copies) == 5
    reference = copies[0].read_text()
    for copy in copies[1:]:
        assert copy.read_text() == reference, copy
