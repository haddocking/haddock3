"""Integration tests for the coarse-grained nucleic acid chain break detection."""

import math
import re
import tempfile
import warnings
from pathlib import Path

import pytest
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from haddock.libs.libsubprocess import CNSJob
from haddock.modules.topology.topocg import DEFAULT_CONFIG as DEFAULT_TOPOCG_CONFIG
from haddock.modules.topology.topocg import HaddockModule as TopocgModule
from haddock.modules.topology.topocg import generate_topology

from . import CNS_EXEC, GOLDEN_DATA


ATOM_RE = re.compile(r"^(\d+) '([^']*)' '([^']*)' '[^']*' '([^']*)'")


@pytest.fixture
def topocg_module(monkeypatch):
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        yield TopocgModule(
            order=1, path=Path("."), initial_params=DEFAULT_TOPOCG_CONFIG
        )


def remove_residue(source: Path, target: Path, chain: str, resid: int) -> None:
    """Write `source` to `target` without residue `resid` of `chain`."""
    with open(source) as fin, open(target, "w") as fout:
        for line in fin:
            if (
                line.startswith(("ATOM", "HETATM"))
                and line[21] == chain
                and int(line[22:26]) == resid
            ):
                continue
            fout.write(line)


def run_topocg(module: TopocgModule, input_pdb: Path) -> tuple[Path, Path]:
    """Run the topocg CNS script on `input_pdb`, return the CG PSF and PDB."""
    defaults = {k: v for k, v in module.params.items() if k != "mol1"}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        inp = generate_topology(
            input_pdb=input_pdb,
            output_path=module.path,
            recipe_str=module.recipe_str,
            defaults=defaults,
            mol_params=module.params["mol1"],
            default_params_path=module.toppar_path,
            write_to_disk=False,
            force_field=module.params["cgffversion"],
            seed=module.params["iniseed"],
        )
    CNSJob(inp, envvars=module.default_envvars(), cns_exec=CNS_EXEC).run()
    cg_stem = f"{input_pdb.stem}_cg_{module.params['cgffversion']}"
    return Path(f"{cg_stem}.psf"), Path(f"{cg_stem}.pdb")


def backbone_links(psf: Path) -> set[tuple]:
    """BB3(i)-BB1(j) bonds of a CNS mtf file, as ((chain, i), (chain, j))."""
    atoms, bonds, section = {}, [], None
    with open(psf) as fin:
        for line in fin:
            if not line.strip():
                section = None
            elif line.startswith("_cns_mtf_atom."):
                section = "atom"
            elif line.startswith("_cns_mtf_bond."):
                section = "bond"
            elif section == "atom" and (match := ATOM_RE.match(line)):
                idx, segid, resid, name = match.groups()
                atoms[int(idx)] = (segid.strip(), int(resid), name)
            elif section == "bond" and not line.startswith("-1"):
                bonds.append(tuple(int(idx) for idx in line.split()))
    links = set()
    for idx1, idx2 in bonds:
        for at1, at2 in ((atoms[idx1], atoms[idx2]), (atoms[idx2], atoms[idx1])):
            if at1[2] == "BB3" and at2[2] == "BB1" and at1[:2] != at2[:2]:
                links.add((at1[:2], at2[:2]))
    return links


def bead_distance(pdb: Path, bead1: tuple, bead2: tuple) -> float:
    """Distance between two (chain, resid, name) beads of `pdb`."""
    coords = {}
    with open(pdb) as fin:
        for line in fin:
            if line.startswith("ATOM"):
                key = (line[72:76].strip(), int(line[22:26]), line[12:16].strip())
                coords[key] = [float(line[i : i + 8]) for i in (30, 38, 46)]
    return math.dist(coords[bead1], coords[bead2])


@pytest.mark.parametrize(
    "structure,removed",
    [("dna_3GO3_aa.pdb", 4), ("rna_1SDR_aa.pdb", 6)],
)
def test_topocg_detects_single_nucleotide_gap(topocg_module, structure, removed):
    """A single missing nucleotide is detected as a CG chain break."""
    gapped = Path(Path(structure).stem.replace("_aa", "_gap") + ".pdb")
    remove_residue(Path(GOLDEN_DATA, structure), gapped, "A", removed)

    psf, cg_pdb = run_topocg(topocg_module, gapped)

    links = backbone_links(psf)
    assert (("A", removed - 1), ("A", removed + 1)) not in links
    # The neighbouring linkages are kept
    assert (("A", removed - 2), ("A", removed - 1)) in links
    assert (("A", removed + 1), ("A", removed + 2)) in links
    # The gap must not be closed by the minimization in generate-topology.cns
    gap = bead_distance(cg_pdb, ("A", removed - 1, "BB3"), ("A", removed + 1, "BB1"))
    assert gap > 7.0


@pytest.mark.parametrize("structure", ["dna_3GO3_aa.pdb", "rna_1SDR_aa.pdb"])
def test_topocg_no_break_in_intact_strands(topocg_module, structure):
    """Intact nucleic acid strands are not broken."""
    psf, _ = run_topocg(topocg_module, Path(GOLDEN_DATA, structure))
    links = backbone_links(psf)
    residues = sorted({res for link in links for res in link})
    consecutive = {
        (res1, res2) for res1, res2 in zip(residues, residues[1:]) if res1[0] == res2[0]
    }
    assert links == consecutive
