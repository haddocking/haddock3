"""alascan module."""

import os
import shutil
from pathlib import Path
from typing import Tuple, Dict


from haddock import log
from haddock.core.typing import Any, Optional, Union
from haddock.libs.libalign import get_atoms, load_coords
from haddock.libs.libontology import PDBFile
from haddock.libs.libscan import (
    MutationResult,
    calc_score,
    execute_scan_jobs,
    AddDeltaBFactor as _AddDeltaBFactor,
    ClusterOutputer as _ClusterOutputer,
    write_scan_out as _write_scan_out,
    group_scan_by_cluster as _group_scan_by_cluster,
)
from haddock.libs.libcapri import CAPRI

ATOMS_TO_BE_MUTATED = ["C", "N", "CA", "O", "CB"]

RES_CODES = dict(
    [
        ("CYS", "C"),
        ("ASP", "D"),
        ("SER", "S"),
        ("GLN", "Q"),
        ("LYS", "K"),
        ("ILE", "I"),
        ("PRO", "P"),
        ("THR", "T"),
        ("PHE", "F"),
        ("ASN", "N"),
        ("GLY", "G"),
        ("HIS", "H"),
        ("LEU", "L"),
        ("ARG", "R"),
        ("TRP", "W"),
        ("ALA", "A"),
        ("VAL", "V"),
        ("GLU", "E"),
        ("TYR", "Y"),
        ("MET", "M"),
        ("ALY", "K"),
        ("ASH", "D"),
        ("CFE", "C"),
        ("CSP", "C"),
        ("CYC", "C"),
        ("CYF", "C"),
        ("CYM", "C"),
        ("DDZ", "A"),
        ("GLH", "E"),
        ("HLY", "P"),
        ("HY3", "P"),
        ("HYP", "P"),
        ("M3L", "K"),
        ("MLY", "K"),
        ("MLZ", "K"),
        ("MSE", "M"),
        ("NEP", "H"),
        ("PNS", "S"),
        ("PTR", "Y"),
        ("SEP", "S"),
        ("TOP", "T"),
        ("TYP", "Y"),
        ("TYS", "Y"),
        ("CIR", "R"),
    ]
)


def mutate(pdb_f, target_chain, target_resid, mut_resname):
    """
    Mutate a residue in a PDB file into a different residue.

    Parameters
    ----------
    pdb_f : str
        Path to the pdb file.

    target_chain : str
        Chain of the residue to be mutated.

    target_resid : int
        Residue number of the residue to be mutated.

    mut_resname : str
        Residue name of the residue to be mutated.

    Returns
    -------
    mut_pdb_fname : str
        Path to the mutated pdb file.
    """
    mut_pdb_l = []
    resname = ""
    with open(pdb_f, "r") as fh:
        for line in fh.readlines():
            if line.startswith("ATOM"):
                chain = line[21]
                resid = int(line[22:26])
                atom_name = line[12:16].strip()
                if target_chain == chain and target_resid == resid:
                    if not resname:
                        resname = line[17:20].strip()
                    if atom_name in ATOMS_TO_BE_MUTATED:
                        # mutate
                        line = line[:17] + mut_resname + line[20:]
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
    """Manage the generation of alascan outputs for cluster-based analysis."""

    module_name = "alascan"
    default_scan_residue = "ALA"
    sort_columns = ["chain", "resid"]
    zscore_reference = "residues"

    def _identity_columns(self):
        return ["chain", "resid", "resname", "full_resname"]

    def _identity_row(self, ident, clt_res_dt):
        parts = ident.split("-")
        return [parts[0], int(parts[1]), parts[2], ident]


class AddDeltaBFactor(_AddDeltaBFactor):
    """Add alascan delta score in the b-factor column of a PDB."""

    module_name = "alascan"


def group_scan_by_cluster(models, results_by_model):
    """Group alascan data per cluster, keyed by residue."""
    return _group_scan_by_cluster(
        models,
        results_by_model,
        ident_builder=lambda r: f"{r.chain}-{r.resid}-{r.ori_resname}",
    )


def write_scan_out(results, model_id):
    """Save alascan mutation results for one model to a tsv file."""
    _write_scan_out(
        results,
        model_id,
        module_name="alascan",
        sort_columns=["chain", "res"],
        zscore_reference="residues",
    )


class InterfaceScanner:
    """Scan interface of a model to get tartget residues and create
    corresponding mutation jobs.
    """

    def __init__(
        self,
        model: Union[str, Path, Any],
        mutation_res: str = "ALA",
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
        mutation_res : str
            Target residue for mutation (default: "ALA")
        params : dict, optional
            Additional parameters for interface detection
            (list on top of alascan/__inint__.py)
        library_mode : bool
            If True, execute mutations sequentially inside InterfaceScanner.run()
            If False, just prepare jobs - execution will be taken care of in init.py
            of alascan module by haddock Engine
        """
        self.model = model
        self.mutation_res = mutation_res
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

    @staticmethod
    def _iter_mutations(interface, resname_dict, mutation_res):
        """Yield each mutation to perform for the interface residues.

        Flattens the interface (chain -> residues) into a single stream of
        ``(chain, resid, ori_resname, target_resname)`` tuples, skipping no-op
        mutations (target residue equal to the original, e.g. ``ALA -> ALA``).

        Parameters
        ----------
        interface : dict
            Mapping of chain id to the list of interface residue numbers.
        resname_dict : dict
            Mapping of ``"{chain}-{resid}"`` to the residue name.
        mutation_res : str
            Residue to mutate every interface residue into (e.g. ``ALA``).

        Yields
        ------
        tuple
            ``(chain, resid, ori_resname, target_resname)`` for each mutation.
        """
        for chain, residues in interface.items():
            for res in residues:
                ori_resname = resname_dict[f"{chain}-{res}"]
                # Skip no-op mutation (e.g. ALA -> ALA)
                if ori_resname == mutation_res:
                    continue
                yield chain, res, ori_resname, mutation_res

    def run(self):
        """
        Get interface residues and create mutation jobs.
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

            # Determine target residues: get interface, then apply user filers, if given
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

            # get all atoms of the model to verifiy residue type down the line
            resname_dict = {}
            for chain, resid, _atom, resname in coords.keys():
                key = f"{chain}-{resid}"
                if key not in resname_dict:
                    resname_dict[key] = resname

            # Create mutation
            output_mutants = self.params.get("output_mutants", False)
            for chain, res, ori_resname, end_resname in self._iter_mutations(
                interface, resname_dict, self.mutation_res
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

            if not self.library_mode:
                # return point_mutations_jobs back to alascan/__init__.py, which
                # hands them to a haddock Engine
                return self.point_mutations_jobs

            # Library mode: run the mutation jobs through a haddock Engine
            # (exactly as the module does in Step 2) instead of executing each
            # job inline, so they go through the normal job lifecycle.
            log.info(
                f"Executing {len(self.point_mutations_jobs)} "
                f"mutations for {self.model_id}"
            )
            results = execute_scan_jobs(self.point_mutations_jobs, self.params)
            write_scan_out(results, self.model_id)

        except Exception as e:
            log.error(f"Failed to scan model {self.model_id}: {e}")
            raise


class ModelPointMutation:
    """Executes a single point mutation."""

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
            Target residue name for mutation
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
                # othervise keep energy-minimized pdb
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
