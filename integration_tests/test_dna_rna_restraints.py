"""Integration tests for the Watson-Crick base pairing in dna-rna_restraints.cns.

The restraint script is run standalone in CNS on all-atom (topoaa) and
coarse-grained (topocg) topologies of:

- one duplex of 1SDR, an A-form RNA dodecamer (12 Watson-Crick pairs);
- 3GO3, B-DNA from an echinomycin complex, with 5 Watson-Crick and
  3 Hoogsteen pairs.

The expected pairs were assigned from the hydrogen bonds of the all-atom
models, independently of the CNS script.
"""

import re
from collections import Counter
from pathlib import Path

import pytest

from haddock.libs.libsubprocess import CNSJob
from haddock.modules.refinement.emref import RECIPE_PATH as EMREF_PATH

from integration_tests import CNS_EXEC, GOLDEN_DATA


# (purine/pyrimidine, partner) as (segid + resid) of the A/DA or C/DC base
RNA_1SDR_WC = {
    ("A2", "B23"), ("A3", "B22"), ("A6", "B19"), ("A11", "B14"),
    ("B13", "A12"), ("B15", "A10"), ("B16", "A9"), ("B17", "A8"),
    ("B18", "A7"), ("B20", "A5"), ("B21", "A4"), ("B24", "A1"),
    }
DNA_3GO3_WC = {
    ("A1", "B108"), ("A2", "B107"), ("A6", "B103"),
    ("B102", "A7"), ("B106", "A3"),
    }
DNA_3GO3_HOOGSTEEN = {("A5", "B104"), ("B101", "A8"), ("B105", "A4")}

# restraints defined per pair: all-atom A-T 4, C-G 5; CG A-T 2, C-G 3
RESTRAINTS_PER_PAIR = {
    "aa": {"ADE": 4, "CYT": 5},
    "cg": {"ADE": 2, "CYT": 3},
    }

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
    r"set-i-atoms\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+"
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


def restrained_pairs(output: str) -> Counter:
    """Count the NOE restraints per (base, partner) pair."""
    noe_table = output.split("noe print threshold")[-1]
    pairs: Counter = Counter()
    for seg1, res1, resn1, seg2, res2 in NOE_RE.findall(noe_table):
        base_type = "ADE" if resn1 in ("A", "DA") else "CYT"
        pairs[(f"{seg1}{res1}", f"{seg2}{res2}", base_type)] += 1
    return pairs


@pytest.mark.parametrize(
    "name,expected,excluded",
    [
        ("rna_1SDR_aa", RNA_1SDR_WC, set()),
        ("rna_1SDR_cg", RNA_1SDR_WC, set()),
        ("dna_3GO3_aa", DNA_3GO3_WC, DNA_3GO3_HOOGSTEEN),
        ("dna_3GO3_cg", DNA_3GO3_WC, DNA_3GO3_HOOGSTEEN),
        ],
    )
def test_watson_crick_pairs(name, expected, excluded):
    """Every Watson-Crick pair is restrained, nothing else is."""
    output = run_dna_rna_restraints(name)

    assert "ABORT mode" not in output
    is_cg = name.endswith("_cg")
    cg_detected = "coarse-grained model detected" in output
    assert cg_detected == is_cg

    pairs = restrained_pairs(output)
    found = {(base, partner) for base, partner, _ in pairs}
    assert found == expected
    assert not found & excluded

    per_pair = RESTRAINTS_PER_PAIR["cg" if is_cg else "aa"]
    for (_, _, base_type), n_restraints in pairs.items():
        assert n_restraints == per_pair[base_type]


def test_dna_rna_restraints_copies_identical():
    """All modules must ship the same dna-rna_restraints.cns."""
    modules_path = EMREF_PATH.parents[1]
    copies = sorted(modules_path.glob("*/*/cns/dna-rna_restraints.cns"))
    assert len(copies) == 5
    reference = copies[0].read_text()
    for copy in copies[1:]:
        assert copy.read_text() == reference, copy
