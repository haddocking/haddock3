"""rnascan module."""

import os
import shutil
from pathlib import Path
from typing import List, Dict


from haddock import log
from haddock.core.typing import Any, Optional, Union
from haddock.libs.libalign import get_atoms, load_coords
from haddock.libs.libscan import (
    MutationResult,  # noqa: F401  re-exported for the module's public API
    calc_score,
    AddDeltaBFactor as _AddDeltaBFactor,
    ClusterOutputer as _ClusterOutputer,
    write_scan_out as _write_scan_out,
    group_scan_by_cluster as _group_scan_by_cluster,
    BaseInterfaceScanner as _BaseInterfaceScanner,
    ModelPointMutation as _ModelPointMutation,
    PURINE_BASE_ATOMS,  # noqa: F401  re-exported for the module's public API
    PYRIMIDINE_BASE_ATOMS,  # noqa: F401  re-exported for the module's public API
    norm_atom_name as _norm_atom_name,
    get_atoms_to_keep as _get_atoms_to_keep,
    validate_scan_bases as _validate_scan_bases,
    filter_interface,
    build_resname_dict,
)
from haddock.libs.libcapri import CAPRI

# Heavy atoms of the ribose-phosphate backbone that are common to all RNA
# nucleotides. These atoms are always preserved when mutating a base. The
# ribose 2'-hydroxyl oxygen (O2') is part of the backbone for RNA. Both the
# modern (OP1/OP2) and legacy (O1P/O2P) phosphate oxygen namings are accepted
# so that models from either convention are handled correctly.
BACKBONE_ATOMS = [
    "P",
    "OP1",
    "OP2",
    "O1P",
    "O2P",
    "O5'",
    "C5'",
    "C4'",
    "O4'",
    "C1'",
    "C2'",
    "O2'",
    "C3'",
    "O3'",
]

# RNA residues that can be scanned/mutated by this module.
RNA_RESIDUES = ("A", "C", "G", "U")

# Ring-type classification of the RNA bases.
PURINES = ("A", "G")
PYRIMIDINES = ("C", "U")

# Base ring atoms and glycosidic anchors are chemistry-invariant and shared with
# dnascan; they live in libscan (PURINE_BASE_ATOMS / PYRIMIDINE_BASE_ATOMS are
# re-exported above so this module's public names are preserved).


def get_atoms_to_keep(ori_resname: str, target_resname: str) -> Dict[str, str]:
    """Return the atoms to preserve when mutating one RNA base into another.

    Thin wrapper around :func:`haddock.libs.libscan.get_atoms_to_keep` bound to
    this module's ribose-phosphate backbone and RNA ring-type classification.
    """
    return _get_atoms_to_keep(
        ori_resname,
        target_resname,
        backbone_atoms=BACKBONE_ATOMS,
        purines=PURINES,
        pyrimidines=PYRIMIDINES,
    )


# Default set of target bases tested for each selected interface nucleotide.
DEFAULT_SCAN_BASES = list(RNA_RESIDUES)


def validate_scan_bases(scan_bases: List[str]) -> List[str]:
    """Validate and normalise the list of target RNA bases.

    Only the canonical RNA residue names (``A``, ``C``, ``G``, ``U``) are
    accepted, in a case-insensitive manner. Two-letter names prefixed with
    ``D`` (``DA``, ``DC``, ``DG``, ``DT``) denote DNA residues in HADDOCK and
    are therefore rejected by this RNA-specific module.

    Thin wrapper around :func:`haddock.libs.libscan.validate_scan_bases`.
    """
    return _validate_scan_bases(
        scan_bases,
        RNA_RESIDUES,
        base_kind="RNA",
        allowed_note=(
            f"RNA residue names {', '.join(RNA_RESIDUES)} (two-letter names "
            "such as DA, DC, DG, DT denote DNA and are not supported)."
        ),
    )


def mutate(pdb_f, target_chain, target_resid, mut_resname):
    """
    Mutate an RNA base in a PDB file into a different base.

    The ribose-phosphate backbone is always kept for the mutated nucleotide.
    When the original and target bases share a ring type (both purines or both
    pyrimidines) the common base ring atoms are kept as well to preserve the
    base orientation (see ``get_atoms_to_keep``). For cross-type mutations
    (purine <-> pyrimidine) the three glycosidic-region anchor atoms are kept and
    renamed to their counterpart in the target ring (pyrimidine N1/C2/C6 <->
    purine N9/C4/C8). The remaining base atoms are dropped and rebuilt by CNS
    during scoring.

    Parameters
    ----------
    pdb_f : str
        Path to the pdb file.

    target_chain : str
        Chain of the nucleotide to be mutated.

    target_resid : int
        Residue number of the nucleotide to be mutated.

    mut_resname : str
        Residue name of the target base (e.g. ``A``, ``C``, ``G``, ``U``).

    Returns
    -------
    mut_pdb_fname : str
        Path to the mutated pdb file.
    """
    mut_pdb_l = []
    resname = ""
    atoms_to_keep: Dict[str, str] = {}
    # RNA residue names are shorter than 3 characters, so they must be
    # right-justified to keep the PDB columns (18-20) aligned.
    resname_field = mut_resname.rjust(3)
    with open(pdb_f, "r") as fh:
        for line in fh:
            if line.startswith("ATOM"):
                chain = line[21]
                resid = int(line[22:26])
                atom_name = _norm_atom_name(line[12:16].strip())
                if target_chain == chain and target_resid == resid:
                    if not resname:
                        resname = line[17:20].strip()
                        # Determine which atoms to keep now that the original
                        # base is known (depends on the ori/target ring types).
                        atoms_to_keep = get_atoms_to_keep(resname, mut_resname)
                    if atom_name in atoms_to_keep:
                        # mutate the residue name
                        line = line[:17] + resname_field + line[20:]
                        # rename the atom if it maps to a different name in the
                        # target ring system (cross-type anchor atoms)
                        new_atom_name = atoms_to_keep[atom_name]
                        if new_atom_name != atom_name:
                            new_field = line[12:16].replace(atom_name, new_atom_name, 1)
                            line = line[:12] + new_field + line[16:]
                        mut_pdb_l.append(line)
                else:
                    mut_pdb_l.append(line)
    # RNA residue names are already the one-letter codes used in identifiers.
    if resname not in RNA_RESIDUES or mut_resname not in RNA_RESIDUES:
        raise KeyError(f"Could not mutate {resname} into {mut_resname}.")
    mut_id = f"{resname}{target_resid}{mut_resname}"
    mut_pdb_fname = Path(pdb_f.name.replace(".pdb", f"-{target_chain}_{mut_id}.pdb"))
    with open(mut_pdb_fname, "w") as fh:
        fh.write("".join(mut_pdb_l))
    return mut_pdb_fname


class ClusterOutputer(_ClusterOutputer):
    """Manage the generation of rnascan outputs for cluster-based analysis."""

    module_name = "rnascan"
    default_scan_residue = "RNA base"
    sort_columns = ["chain", "resid", "target_resname"]
    zscore_reference = "mutations"

    def _identity_columns(self):
        return ["chain", "resid", "resname", "target_resname", "full_resname"]

    def _identity_row(self, ident, clt_res_dt):
        # ident is "<chain>-<resid>-<ori_resname>-<target_resname>"
        chain, resid, resname, target_resname = ident.split("-")
        return [chain, int(resid), resname, target_resname, ident]


class AddDeltaBFactor(_AddDeltaBFactor):
    """Add rnascan delta score in the b-factor column of a PDB."""

    module_name = "rnascan"


def group_scan_by_cluster(models, results_by_model):
    """Group rnascan data per cluster, keyed by mutation (base included)."""
    return _group_scan_by_cluster(
        models,
        results_by_model,
        ident_builder=lambda r: (
            f"{r.chain}-{r.resid}-{r.ori_resname}-{r.target_resname}"
        ),
    )


def write_scan_out(results, model_id):
    """Save rnascan mutation results for one model to a tsv file."""
    _write_scan_out(
        results,
        model_id,
        module_name="rnascan",
        sort_columns=["chain", "res", "end_resname"],
        zscore_reference="mutations",
    )


class InterfaceScanner(_BaseInterfaceScanner):
    """Scan interface of a model to get target nucleotides and create
    corresponding mutation jobs.
    """

    def __init__(
        self,
        model: Union[str, Path, Any],
        scan_bases: Optional[List[str]] = None,
        params: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Initialize InterfaceScanner for a single model.

        Parameters
        ----------
        model : str, Path, or model object
            HADDOCK ``PDBFile`` model object or a path to a PDB file.
        scan_bases : list of str, optional
            Target bases tested for each interface nucleotide
            (default: A, C, G, T)
        params : dict, optional
            Additional parameters for interface detection
            (list on top of rnascan/__init__.py)
        """
        super().__init__(model, params)
        self.scan_bases = list(scan_bases) if scan_bases else list(DEFAULT_SCAN_BASES)

    @staticmethod
    def _iter_mutations(interface, resname_dict, scan_bases):
        """Yield each valid single-base mutation to perform.

        Flattens the interface (chain -> residues) and the requested target
        bases into a single stream of ``(chain, resid, ori_resname,
        target_resname)`` tuples, skipping non-RNA residues and no-op mutations
        (target base equal to the original, e.g. ``A -> A``).

        Parameters
        ----------
        interface : dict
            Mapping of chain id to the list of interface residue numbers.
        resname_dict : dict
            Mapping of ``"{chain}-{resid}"`` to the residue name.
        scan_bases : list of str
            Target bases to scan each nucleotide into.

        Yields
        ------
        tuple
            ``(chain, resid, ori_resname, target_resname)`` for each mutation.
        """
        for chain, residues in interface.items():
            for res in residues:
                ori_resname = resname_dict[f"{chain}-{res}"]
                # Only scan RNA nucleotides, skip protein/other residues
                if ori_resname not in RNA_RESIDUES:
                    continue
                for end_resname in scan_bases:
                    # Skip no-op mutation (e.g. A -> A)
                    if ori_resname == end_resname:
                        continue
                    yield chain, res, ori_resname, end_resname

    def run(self):
        """
        Get interface nucleotides and create the mutation jobs for this model.

        The jobs are returned (not executed): the caller hands them to a haddock
        Engine so that all mutations are scheduled together.

        Returns
        -------
        List[ModelPointMutation]
            The mutation jobs to perform for this model.
        """
        try:
            # Calculate native scores
            sc_dir = f"haddock3-score-{self.model_id}-{os.getpid()}"
            try:
                native_scores = calc_score(
                    self.model_path,
                    run_dir=sc_dir,
                    outputpdb=False,
                    ligand_param_fname=self.ligand_param_fname,
                    ligand_top_fname=self.ligand_top_fname,
                )
            finally:
                if os.path.exists(sc_dir):
                    shutil.rmtree(sc_dir)

            # Load coordinates
            atoms = get_atoms(self.model_path)
            coords, _chain_ranges = load_coords(
                self.model_path,
                atoms,
                add_resname=True,
            )

            # Determine target nucleotides: get interface, then apply user filters
            cutoff = self.params.get("int_cutoff", 5.0)
            interface = CAPRI.identify_interface(self.model_path, cutoff=cutoff)
            interface = filter_interface(
                interface, self.filter_resdic, self.params.get("chains", [])
            )

            # residue type lookup used to verify residue type down the line
            resname_dict = build_resname_dict(coords)

            # Create mutations
            output_mutants = self.params.get("output_mutants", False)
            for chain, res, ori_resname, end_resname in self._iter_mutations(
                interface, resname_dict, self.scan_bases
            ):
                job = ModelPointMutation(
                    model_path=self.model_path,
                    model_id=self.model_id,
                    chain=chain,
                    resid=res,
                    ori_resname=ori_resname,
                    target_resname=end_resname,
                    native_scores=native_scores,
                    output_mutants=output_mutants,
                    ligand_param_fname=self.ligand_param_fname,
                    ligand_top_fname=self.ligand_top_fname,
                )
                self.point_mutations_jobs.append(job)

            return self.point_mutations_jobs

        except Exception as e:
            log.error(f"Failed to scan model {self.model_id}: {e}")
            raise


class ModelPointMutation(_ModelPointMutation):
    """Execute a single rnascan (RNA base) point mutation.

    Shares its scoring flow with :class:`haddock.libs.libscan.ModelPointMutation`
    and only forwards this module's ``mutate``/``calc_score`` (which stay
    patchable in the module namespace).
    """

    def run(self):
        """Execute the point mutation."""
        return self._run(mutate, calc_score)
