"""rnascan module."""

import os
import shutil
from pathlib import Path
from typing import List, Tuple, Dict


from haddock import log
from haddock.core.exceptions import ConfigurationError
from haddock.core.typing import Any, Optional, Union
from haddock.libs.libalign import get_atoms, load_coords
from haddock.libs.libontology import PDBFile
from haddock.libs.libscan import (
    MutationResult,
    calc_score,
    AddDeltaBFactor as _AddDeltaBFactor,
    ClusterOutputer as _ClusterOutputer,
    write_scan_out as _write_scan_out,
    group_scan_by_cluster as _group_scan_by_cluster,
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

# Base ring atoms shared by both purines (A and G). Preserving them when
# mutating one purine into the other keeps the base plane and glycosidic
# orientation, so CNS only has to rebuild the differing substituents.
PURINE_BASE_ATOMS = ["N9", "C8", "N7", "C5", "C4", "N3", "C2", "N1", "C6"]

# Base ring atoms shared by both pyrimidines (C and U). Same rationale as for
# the purines above.
PYRIMIDINE_BASE_ATOMS = ["N1", "C2", "O2", "N3", "C4", "C5", "C6"]

# RNA residues that can be scanned/mutated by this module.
RNA_RESIDUES = ("A", "C", "G", "U")

# Ring-type classification of the RNA bases.
PURINES = ("A", "G")
PYRIMIDINES = ("C", "U")

# Glycosidic-region anchor atoms kept (and renamed) on cross-type mutations
# (purine <-> pyrimidine). Because the purine and pyrimidine ring systems do
# not share atom names, these three atoms are preserved and renamed to their
# counterpart in the target ring so that CNS rebuilds the new base with the
# correct glycosidic orientation (and syn/anti conformation). The
# correspondence is:
#   pyrimidine N1 <-> purine N9   (glycosidic nitrogen)
#   pyrimidine C2 <-> purine C4
#   pyrimidine C6 <-> purine C8
PYRIMIDINE_TO_PURINE_ANCHORS = {"N1": "N9", "C2": "C4", "C6": "C8"}
PURINE_TO_PYRIMIDINE_ANCHORS = {"N9": "N1", "C4": "C2", "C8": "C6"}


def get_atoms_to_keep(ori_resname: str, target_resname: str) -> Dict[str, str]:
    """Return the atoms to preserve when mutating one base into another.

    The result maps each atom name to keep (as found in the original residue)
    to the atom name it must be written with in the mutated residue. Backbone
    atoms and same-ring-type base atoms keep their name (mapped to themselves).

    The ribose-phosphate backbone is always kept. In addition, when the
    original and target bases are of the same ring type (both purines or both
    pyrimidines), the base ring atoms common to that type are kept as well, so
    that the base orientation is preserved and CNS only rebuilds the differing
    substituents. For cross-type mutations (purine <-> pyrimidine) the three
    glycosidic-region anchor atoms are kept and renamed to their counterpart in
    the target ring system (pyrimidine N1/C2/C6 <-> purine N9/C4/C8), so that the
    base orientation is preserved and CNS rebuilds the rest of the new base.

    Parameters
    ----------
    ori_resname : str
        Original (wild-type) residue name.
    target_resname : str
        Target base residue name.

    Returns
    -------
    dict of str -> str
        Mapping of original atom name -> atom name to write for the mutated
        nucleotide.
    """
    # Backbone is always kept, with unchanged atom names.
    atoms_to_keep: Dict[str, str] = {atom: atom for atom in BACKBONE_ATOMS}
    if ori_resname in PURINES and target_resname in PURINES:
        atoms_to_keep.update({atom: atom for atom in PURINE_BASE_ATOMS})
    elif ori_resname in PYRIMIDINES and target_resname in PYRIMIDINES:
        atoms_to_keep.update({atom: atom for atom in PYRIMIDINE_BASE_ATOMS})
    elif ori_resname in PYRIMIDINES and target_resname in PURINES:
        atoms_to_keep.update(PYRIMIDINE_TO_PURINE_ANCHORS)
    elif ori_resname in PURINES and target_resname in PYRIMIDINES:
        atoms_to_keep.update(PURINE_TO_PYRIMIDINE_ANCHORS)
    return atoms_to_keep


# Default set of target bases tested for each selected interface nucleotide.
DEFAULT_SCAN_BASES = ["A", "C", "G", "U"]


def validate_scan_bases(scan_bases: List[str]) -> List[str]:
    """Validate and normalise the list of target RNA bases.

    Only the canonical RNA residue names (``A``, ``C``, ``G``, ``U``) are
    accepted, in a case-insensitive manner. Two-letter names prefixed with
    ``D`` (``DA``, ``DC``, ``DG``, ``DT``) denote DNA residues in HADDOCK and
    are therefore rejected by this RNA-specific module.

    Parameters
    ----------
    scan_bases : list of str
        Target bases requested by the user.

    Returns
    -------
    list of str
        Normalised, de-duplicated list of canonical RNA residue names, in the
        order first seen.

    Raises
    ------
    ConfigurationError
        If ``scan_bases`` is empty or contains an unrecognised base.
    """
    if not scan_bases:
        raise ConfigurationError(
            "'scan_bases' must contain at least one RNA base "
            f"among {', '.join(RNA_RESIDUES)}."
        )
    normalised: List[str] = []
    for base in scan_bases:
        canonical = str(base).strip().upper()
        if canonical not in RNA_RESIDUES:
            raise ConfigurationError(
                f"Invalid 'scan_bases' entry {base!r}. Allowed values are the "
                f"RNA residue names {', '.join(RNA_RESIDUES)} (two-letter names "
                "such as DA, DC, DG, DT denote DNA and are not supported)."
            )
        if canonical not in normalised:
            normalised.append(canonical)
    return normalised


RES_CODES = dict(
    [
        # RNA
        ("A", "A"),
        ("C", "C"),
        ("G", "G"),
        ("U", "U"),
    ]
)


def _norm_atom_name(atom_name: str) -> str:
    """Normalise atom names so that primes are always written as `'`."""
    return atom_name.replace("*", "'")


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
    atoms_to_keep: List[str] = []
    # RNA residue names are shorter than 3 characters, so they must be
    # right-justified to keep the PDB columns (18-20) aligned.
    resname_field = mut_resname.rjust(3)
    with open(pdb_f, "r") as fh:
        for line in fh.readlines():
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
    try:
        mut_id = f"{RES_CODES[resname]}{target_resid}{RES_CODES[mut_resname]}"
    except KeyError:
        raise KeyError(f"Could not mutate {resname} into {mut_resname}.")
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


class InterfaceScanner:
    """Scan interface of a model to get target nucleotides and create
    corresponding mutation jobs.
    """

    def __init__(
        self,
        model: Union[str, Path, Any],
        scan_bases: Optional[List[str]] = None,
        params: Optional[Dict[str, Any]] = None,
        library_mode: bool = True,
    ) -> None:
        """
        Initialize InterfaceScanner for a single model.

        Parameters
        ----------
        model : str, Path, or model object
            HADDOCK model object (if library_mode = False) or
            Path to PDB file (if library mode = True)
        scan_bases : list of str, optional
            Target bases tested for each interface nucleotide
            (default: A, C, G, T)
        params : dict, optional
            Additional parameters for interface detection
            (list on top of rnascan/__init__.py)
        library_mode : bool
            If True, execute mutations sequentially inside InterfaceScanner.run()
            If False, just prepare jobs - execution will be taken care of in
            __init__.py of rnascan module by haddock Engine
        """
        self.model = model
        self.scan_bases = list(scan_bases) if scan_bases else list(DEFAULT_SCAN_BASES)
        self.library_mode = library_mode
        self.params = params or {}
        self.point_mutations_jobs = []
        self.ligand_param_fname = self.params.get("ligand_param_fname", "")
        self.ligand_top_fname = self.params.get("ligand_top_fname", "")
        self.filter_resdic = {
            key[-1]: value
            for key, value in self.params.items()
            if key.startswith("resdic")
        }
        if isinstance(model, PDBFile):
            self.model_path = model.rel_path
            self.model_id = model.file_name.removesuffix(".pdb")
        else:
            # for library mode
            self.model_path = Path(model)
            self.model_id = self.model_path.stem

    def run(self):
        """
        Get interface nucleotides and create mutation jobs.
        If library_mode=True, also execute the mutations sequentially.

        Returns
        -------
        Optional[List[ModelPointMutation]]
            List of mutation jobs if library_mode=False, None if library_mode=True
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
            # Get all interface residues
            cutoff = self.params.get("int_cutoff", 5.0)
            interface = CAPRI.identify_interface(self.model_path, cutoff=cutoff)

            # get user_chains for the check down the line
            user_chains = self.params.get("chains", [])

            # if user defined target residues, check they are in the interface
            if self.filter_resdic != {"_": []}:
                filtered_interface = {}
                for chain in self.filter_resdic:
                    if chain in interface:
                        # Search for the intersection of user queried residues and interface residues
                        user_res_valid = list(
                            set(self.filter_resdic[chain]).intersection(
                                set(interface[chain])
                            )
                        )
                        # If at least one residue must be analyzed, add it to residues to be scanned
                        if user_res_valid:
                            filtered_interface[chain] = user_res_valid
                interface = filtered_interface

            # if (user defined target chains) & (no user target residues) - do use user chains
            elif user_chains:
                interface = {
                    chain: res
                    for chain, res in interface.items()
                    if chain in user_chains
                }

            # get all atoms of the model to verify residue type down the line
            resname_dict = {}
            for chain, resid, _atom, resname in coords.keys():
                key = f"{chain}-{resid}"
                if key not in resname_dict:
                    resname_dict[key] = resname

            # Create mutations
            output_mutants = self.params.get("output_mutants", False)
            for chain in interface:
                for res in interface[chain]:
                    ori_resname = resname_dict[f"{chain}-{res}"]
                    # Only scan RNA nucleotides, skip protein/other residues
                    if ori_resname not in RNA_RESIDUES:
                        continue
                    # Test each requested target base
                    for end_resname in self.scan_bases:
                        # Skip if target base is the same as original
                        # (e.g. skip A->A)
                        if ori_resname == end_resname:
                            continue
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

            # Execute jobs if in library mode
            if self.library_mode:
                log.info(
                    f"Executing {len(self.point_mutations_jobs)} "
                    f"mutations for {self.model_id}"
                )
                results = []
                total = len(self.point_mutations_jobs)

                for i, job in enumerate(self.point_mutations_jobs, 1):
                    log.info(
                        f"Processing mutation {i}/{total}: "
                        f"{job.chain}:{job.resid} {job.ori_resname}"
                        f"->{job.target_resname}"
                    )
                    result = job.run()
                    results.append(result)
                    if result.success:
                        log.info(f"Delta score: {result.delta_scores[0]:.2f}")
                    else:
                        log.warning(f"Failed: {result.error_msg}")

                write_scan_out(results, self.model_id)
            else:
                # return point_mutations_jobs back to rnascan/__init__.py
                return self.point_mutations_jobs

        except Exception as e:
            log.error(f"Failed to scan model {self.model_id}: {e}")
            raise


class ModelPointMutation:
    """Executes a single base point mutation."""

    def __init__(
        self,
        model_path: Path,
        model_id: str,
        chain: str,
        resid: int,
        ori_resname: str,
        target_resname: str,
        native_scores: Tuple[float, float, float, float, float],
        output_mutants: bool = False,
        ligand_param_fname: Union[Path, str] = "",
        ligand_top_fname: Union[Path, str] = "",
    ) -> None:
        """
        Initialize a single point mutation job.

        Parameters
        ----------
        model_path : Path
            Path to the PDB file
        model_id : str
            Identifier for the model
        chain : str
            Chain identifier
        resid : int
            Residue number
        ori_resname : str
            Original residue name
        target_resname : str
            Target base name for mutation
        native_scores : tuple
            Native model scores (score, vdw, elec, desolv, bsa)
        output_mutants : bool
            Whether to keep mutant PDB files
        ligand_param_fname : Union[Path, str]
            Path to additional parameter file used by CNS
        ligand_top_fname : Union[Path, str]
            Path to additional topology file used by CNS
        """
        self.model_path = Path(model_path)
        self.model_id = model_id
        self.chain = chain
        self.resid = resid
        self.ori_resname = ori_resname
        self.target_resname = target_resname
        self.native_scores = native_scores
        self.output_mutants = output_mutants
        self.ligand_param_fname = ligand_param_fname
        self.ligand_top_fname = ligand_top_fname

    def run(self):
        """Execute the point mutation."""
        mutation_id = f"{self.model_id}_{self.chain}{self.resid}{self.target_resname}"

        try:
            # Setup working directory
            sc_dir = f"haddock3-score-{mutation_id}"
            os.makedirs(sc_dir, exist_ok=True)

            # Perform point mutation on pdb file
            mut_pdb = mutate(
                self.model_path, self.chain, self.resid, self.target_resname
            )

            # Calculate mutant scores
            mutant_scores = calc_score(
                mut_pdb,
                run_dir=sc_dir,
                outputpdb=self.output_mutants,
                ligand_param_fname=self.ligand_param_fname,
                ligand_top_fname=self.ligand_top_fname,
            )

            # Calculate deltas (native - mutant)
            n_score, n_vdw, n_elec, n_des, n_bsa = self.native_scores
            m_score, m_vdw, m_elec, m_des, m_bsa = mutant_scores
            delta_scores = (
                n_score - m_score,
                n_vdw - m_vdw,
                n_elec - m_elec,
                n_des - m_des,
                n_bsa - m_bsa,
            )

            # Handle output files
            em_mut_pdb = Path(f"{mut_pdb.stem}_hs.pdb")
            if not self.output_mutants:
                # if output_mutants = False, then remove both files
                if os.path.exists(mut_pdb):
                    os.remove(mut_pdb)
                if em_mut_pdb.exists():
                    os.remove(em_mut_pdb)
            else:
                # otherwise keep energy-minimized pdb
                if os.path.exists(em_mut_pdb):
                    shutil.move(em_mut_pdb, mut_pdb)
            # clean up scoring dir
            if os.path.exists(sc_dir):
                shutil.rmtree(sc_dir)

            return MutationResult(
                model_id=self.model_id,
                chain=self.chain,
                resid=self.resid,
                ori_resname=self.ori_resname,
                target_resname=self.target_resname,
                mutant_scores=mutant_scores,
                delta_scores=delta_scores,
                success=True,
            )

        except Exception as e:
            return MutationResult(
                model_id=self.model_id,
                chain=self.chain,
                resid=self.resid,
                ori_resname=self.ori_resname,
                target_resname=self.target_resname,
                mutant_scores=(0, 0, 0, 0, 0),
                delta_scores=(0, 0, 0, 0, 0),
                success=False,
                error_msg=str(e),
            )
