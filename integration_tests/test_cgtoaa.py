import gzip
import shutil
from pathlib import Path

from haddock.libs.libaa2cg import martinize
from haddock.libs.libontology import PDBFile, TopologyFile
from haddock.modules.refinement.cgtoaa import DEFAULT_CONFIG as DEFAULT_CGTOAA_CONFIG
from haddock.modules.refinement.cgtoaa import HaddockModule as CgtoaaModule
from haddock.modules.topology.topoaa import DEFAULT_CONFIG as DEFAULT_TOPOAA_CONFIG
from haddock.modules.topology.topoaa import HaddockModule as TopoaaModule

from integration_tests import GOLDEN_DATA


class MockPreviousIO:
    """Returns a single CG model wired with AA topology and backmapping TBL."""

    def __init__(self, path, cg_pdb, cg_psf, aa_psf, cg_to_aa_tbl):
        self.path = path
        self._cg_pdb = cg_pdb
        self._cg_psf = cg_psf
        self._aa_psf = aa_psf
        self._cg_to_aa_tbl = cg_to_aa_tbl

    def retrieve_models(self, individualize=True):
        cg_pdb_dest = self.path / self._cg_pdb.name
        shutil.copy(self._cg_pdb, cg_pdb_dest)
        cg_psf_dest = self.path / self._cg_psf.name
        shutil.copy(self._cg_psf, cg_psf_dest)

        model = PDBFile(
            file_name=cg_pdb_dest.name,
            path=self.path,
            topology=TopologyFile(
                file_name=cg_psf_dest.name,
                path=self.path,
            ),
        )
        model.aa_topology = TopologyFile(
            file_name=self._aa_psf.name,
            path=self._aa_psf.parent,
        )
        model.cgtoaa_tbl = self._cg_to_aa_tbl
        model.seed = 42
        return [model]

    def output(self):
        return None


def test_cgtoaa_backmapping_restraints_are_read(tmp_path):
    """Backmapping TBL must appear in CNS output, confirming it was read."""
    # generate AA PSF
    topoaa_dir = tmp_path / "topoaa"
    topoaa_dir.mkdir()
    topoaa = TopoaaModule(
        order=0, path=topoaa_dir, initial_params=DEFAULT_TOPOAA_CONFIG
    )
    topoaa.params["molecules"] = [Path(GOLDEN_DATA, "hpr_ensemble.pdb")]
    topoaa.params["mol1"]["prot_segid"] = "B"
    topoaa.run()

    aa_psf_files = sorted(topoaa_dir.glob("*_haddock.psf"))
    assert aa_psf_files, "topoaa produced no AA PSF"
    aa_psf = aa_psf_files[0]
    aa_pdb = aa_psf.with_suffix(".pdb")
    assert aa_pdb.exists(), f"AA PDB not found: {aa_pdb}"

    # generate backmapping tbl
    martinize(aa_pdb, str(topoaa_dir), skipss=True)
    tbl_files = list(topoaa_dir.glob("*_cg_to_aa.tbl"))
    assert tbl_files, "martinize produced no backmapping TBL"
    cg_to_aa_tbl = tbl_files[0]

    # Append a unique marker comment to the TBL, the marker will appear in the
    # CNS output ONLY if the file was actually read.
    _MARKER = "! CGTOAA_TBL_WAS_READ"
    with open(cg_to_aa_tbl, "a") as fh:
        fh.write(f"\n{_MARKER}\n")

    # run cg-to-aa
    cgtoaa_dir = tmp_path / "cgtoaa"
    cgtoaa_dir.mkdir()
    cgtoaa = CgtoaaModule(
        order=1, path=cgtoaa_dir, initial_params=DEFAULT_CGTOAA_CONFIG
    )
    cgtoaa.previous_io = MockPreviousIO(
        path=cgtoaa_dir,
        cg_pdb=Path(GOLDEN_DATA, "hpr_haddock_cg.pdb"),
        cg_psf=Path(GOLDEN_DATA, "hpr_haddock_cg.psf"),
        aa_psf=aa_psf,
        cg_to_aa_tbl=cg_to_aa_tbl,
    )
    cgtoaa.params["debug"] = True
    cgtoaa.run()

    # check TBL was read
    output_pdb = cgtoaa_dir / "cgtoaa_1.pdb"
    assert output_pdb.exists(), "cgtoaa produced no output PDB"

    output_gz = cgtoaa_dir / "cgtoaa_1.out.gz"
    assert output_gz.exists(), "cgtoaa produced no CNS output"

    with gzip.open(output_gz, "rt") as fh:
        cns_output = fh.read()

    assert _MARKER in cns_output, (
        "Marker comment not found in CNS output — "
        "backmapping restraints were NOT read (Issue #1592 regression)"
    )
