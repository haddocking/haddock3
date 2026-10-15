"""Integration tests for the interface restraints in contactairs.cns.

The restraint script is run standalone in CNS on the all-atom topology of an
antibody-antigen complex (segids B, H and L). Every residue with any atom
within 5 A of another molecule must be restrained to all residues of the
other molecules that are within 5 A of its own molecule (see issue #1727).

The expected restraints are computed from the coordinates with numpy,
independently of the CNS script.
"""

import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

from haddock.libs.libsubprocess import CNSJob
from haddock.modules.refinement.flexref import RECIPE_PATH as FLEXREF_PATH

from integration_tests import CNS_EXEC, GOLDEN_DATA, has_cns


STRUCTURE = "ab_ag_BHL"
SEGIDS = ("B", "H", "L")
CUTOFF = 5.0

CNS_INPUT = """
set message=normal echo=off end
structure @GOLDEN:{name}.psf end
coor @GOLDEN:{name}.pdb
evaluate ($data.ncomponents = {ncomponents})
{segids}
noe reset nrestraints=1200000 ceiling=1000 end
@MODULE:contactairs.cns
set echo=on end
noe print threshold=-1.0 end
stop
"""

ATOM_RE = re.compile(r"^\s+(\S+)\s+(\d+)\s+\S+\s+\S+\s*$", re.MULTILINE)


def run_contactairs(name: str, segids: tuple[str, ...]) -> str:
    """Run contactairs.cns on a golden topology, return the output."""
    segid_defs = "\n".join(
        f'evaluate ($prot_segid_{i} = "{segid}")'
        for i, segid in enumerate(segids, start=1)
    )
    job = CNSJob(
        input_file=CNS_INPUT.format(
            name=name, ncomponents=len(segids), segids=segid_defs
        ),
        envvars={
            "MODULE": str(Path(FLEXREF_PATH, "cns")),
            "GOLDEN": str(GOLDEN_DATA),
        },
        cns_exec=CNS_EXEC,
    )
    return job.run().decode()


def parse_restraints(output: str) -> dict:
    """
    Parse the NOE restraints printed by CNS.

    Returns {(segid, resid): set of (segid, resid) partner residues}.
    """
    noe_table = output.split("noe print threshold")[-1]
    restraints: dict = {}
    for block in noe_table.split("set-i-atoms")[1:]:
        i_atoms, j_atoms = block.split("R<average>")[0].split("set-j-atoms")
        (residue,) = {(seg, int(res)) for seg, res in ATOM_RE.findall(i_atoms)}
        assert residue not in restraints, f"duplicate restraint {residue}"
        restraints[residue] = {(seg, int(res)) for seg, res in ATOM_RE.findall(j_atoms)}
    return restraints


def expected_restraints(pdb: Path, segids: tuple[str, ...]) -> dict:
    """Compute the expected interface restraints from the coordinates."""
    residues, coords = [], []
    for line in pdb.read_text().splitlines():
        if line.startswith("ATOM") and line[72:76].strip() in segids:
            residues.append((line[72:76].strip(), int(line[22:26])))
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    xyz = np.array(coords)
    segid_of = np.array([segid for segid, _ in residues])

    # atoms of each segid within the cutoff of any atom of another segid
    near_other = {}
    for segid in segids:
        own, other = segid_of == segid, segid_of != segid
        dist = np.linalg.norm(xyz[own][:, None, :] - xyz[other][None, :, :], axis=-1)
        near = np.zeros(len(residues), dtype=bool)
        near[np.flatnonzero(own)[(dist <= CUTOFF).any(axis=1)]] = True
        near[np.flatnonzero(other)[(dist <= CUTOFF).any(axis=0)]] = True
        near_other[segid] = near

    expected = defaultdict(set)
    for segid in segids:
        interface = {residues[i] for i in np.flatnonzero(near_other[segid])}
        own_interface = {res for res in interface if res[0] == segid}
        partners = interface - own_interface
        for residue in own_interface:
            expected[residue] = partners
    return dict(expected)


@pytest.fixture(scope="module")
def restraints():
    """Restraints defined by contactairs.cns on the golden complex."""
    output = run_contactairs(STRUCTURE, SEGIDS)
    assert "ABORT mode" not in output
    return parse_restraints(output)


@has_cns
def test_contactairs_restraints_defined(restraints):
    """Interface restraints are defined for every molecule (issue #1727)."""
    assert restraints
    assert {segid for segid, _ in restraints} == set(SEGIDS)


@has_cns
def test_contactairs_interface(restraints):
    """Exactly the interface residues are restrained to the partner interface."""
    expected = expected_restraints(Path(GOLDEN_DATA, f"{STRUCTURE}.pdb"), SEGIDS)
    assert set(restraints) == set(expected)
    for residue, partners in restraints.items():
        assert partners == expected[residue], residue
        assert residue[0] not in {segid for segid, _ in partners}


def test_contactairs_copies_identical():
    """All modules must ship the same contactairs.cns."""
    modules_path = FLEXREF_PATH.parents[1]
    copies = sorted(modules_path.glob("*/*/cns/contactairs.cns"))
    assert len(copies) == 4
    reference = copies[0].read_text()
    for copy in copies[1:]:
        assert copy.read_text() == reference, copy
