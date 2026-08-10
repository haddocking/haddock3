"""alascan module."""

import os
import shutil
from pathlib import Path
from typing import Dict


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
    filter_interface,
    build_resname_dict,
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


class InterfaceScanner(_BaseInterfaceScanner):
    """Scan interface of a model to get tartget residues and create
    corresponding mutation jobs.
    """

    def __init__(
        self,
        model: Union[str, Path, Any],
        mutation_res: str = "ALA",
        params: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Initialize InterfaceScanner for a single model.

        Parameters
        ----------
        model : str, Path, or model object
            HADDOCK ``PDBFile`` model object or a path to a PDB file.
        mutation_res : str
            Target residue for mutation (default: "ALA")
        params : dict, optional
            Additional parameters for interface detection
            (list on top of alascan/__inint__.py)
        """
        super().__init__(model, params)
        self.mutation_res = mutation_res

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
        Get interface residues and create the mutation jobs for this model.

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

            # Determine target residues: get interface, then apply user filters
            cutoff = self.params.get("int_cutoff", 5.0)
            interface = CAPRI.identify_interface(self.model_path, cutoff=cutoff)
            interface = filter_interface(
                interface, self.filter_resdic, self.params.get("chains", [])
            )

            # residue type lookup used to verify residue type down the line
            resname_dict = build_resname_dict(coords)

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

            return self.point_mutations_jobs

        except Exception as e:
            log.error(f"Failed to scan model {self.model_id}: {e}")
            raise


class ModelPointMutation(_ModelPointMutation):
    """Execute a single alascan (protein) point mutation.

    Shares its scoring flow with :class:`haddock.libs.libscan.ModelPointMutation`
    and only forwards this module's ``mutate``/``calc_score`` (which stay
    patchable in the module namespace).
    """

    def run(self):
        """Execute the point mutation."""
        return self._run(mutate, calc_score)
