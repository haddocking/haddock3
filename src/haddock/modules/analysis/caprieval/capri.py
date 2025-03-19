"""CAPRI module."""

import copy
import json
import os
import shutil
import tempfile
from itertools import combinations
from math import isnan
from pathlib import Path


os.environ["OPENBLAS_NUM_THREADS"] = "1"

import numpy as np
from pdbtools import pdb_segxchain
from scipy.spatial.distance import cdist

from haddock import log
from haddock.core.defaults import CNS_MODULES
from haddock.core.typing import (
    Any,
    AtomsDict,
    FilePath,
    Iterable,
    NDFloat,
    Optional,
    ParamDict,
    ParamMap,
    Union,
    )
from haddock.gear.config import load as read_config
from haddock.libs.libalign import (
    ALIGNError,
    calc_rmsd,
    centroid,
    check_chains,
    get_align,
    get_atoms,
    kabsch,
    load_coords,
    make_range,
    )
from haddock.libs.libio import write_dic_to_file, write_nested_dic_to_file
from haddock.libs.libontology import PDBFile, PDBPath
from haddock.modules import get_module_steps_folders


WEIGHTS = ["w_elec", "w_vdw", "w_desolv", "w_bsa", "w_air"]


def get_previous_cns_step(sel_steps: list, st_order: int) -> Union[str, None]:
    """
    Get the previous CNS step.

    Parameters
    ----------
    run_path : Path
        Path to the run folder.

    Returns
    -------
    cns_step : str
        Name of the CNS step.
    """
    # get the previous CNS step
    cns_step = None
    # just to be careful, remove steps with more than one underscore
    sel_steps = [step for step in sel_steps if step.count("_") == 1]
    mod = min(st_order - 1, len(sel_steps) - 1)
    # loop
    while mod > -1:
        st_name = sel_steps[mod].split("_")[1]
        if st_name in CNS_MODULES:
            cns_step = sel_steps[mod]
            break
        mod -= 1

    return cns_step


def save_scoring_weights(cns_step: str) -> Path:
    """Save the scoring weights in a json file.

    Parameters
    ----------
    cns_step : str
        Name of the CNS step.

    Returns
    -------
    scoring_params_fname : Path
        Path to the json file.
    """
    cns_params = read_config(Path("..", cns_step, "params.cfg"))
    key = list(cns_params["final_cfg"].keys())[0]
    scoring_pars = {kv: cns_params["final_cfg"][key][kv] for kv in WEIGHTS}

    scoring_params_fname = Path("weights_params.json")
    # write json file
    with open(scoring_params_fname, "w", encoding="utf-8") as jsonf:
        json.dump(
            scoring_pars,
            jsonf,
            indent=4,
        )
    return scoring_params_fname


def load_contacts(
    pdb_f: Union[Path, PDBFile],
    cutoff: float = 5.0,
    numbering_dic: Optional[dict[str, dict[int, int]]] = None,
    model2ref_chain_dict: Optional[dict[str, str]] = None,
) -> set[tuple]:
    """Load residue-based contacts.

    Parameters
    ----------
    pdb_f : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`
        PDB file of the model to have its atoms identified
    cutoff : float, optional
        Cutoff distance for the interface identification.

    Returns
    -------
    set(con_list) : set
        set of unique contacts
    """
    if isinstance(pdb_f, PDBFile):
        pdb_f = pdb_f.rel_path
    # get also side chains atoms
    atoms = get_atoms(pdb_f, full=True)
    ref_coord_dic, _ = load_coords(
        pdb_f,
        atoms,
        numbering_dic=numbering_dic,
        model2ref_chain_dict=model2ref_chain_dict,
    )
    # create coordinate arrays
    coord_arrays: dict[str, NDFloat] = {}
    coord_ids: dict[str, list[int]] = {}
    for atom in ref_coord_dic.keys():
        chain = atom[0]
        if chain not in coord_arrays.keys():  # initialize lists
            coord_arrays[chain], coord_ids[chain] = [], []  # type: ignore
        coord_arrays[chain].append(ref_coord_dic[atom])  # type: ignore
        coord_ids[chain].append(atom[1])  # only the resid is appended
    for chain in coord_arrays.keys():
        coord_arrays[chain] = np.array(coord_arrays[chain])

    # combinations of chains
    unique_chain_combs = list(combinations(sorted(coord_arrays.keys()), 2))

    # calculating contacts
    con_list: list[tuple] = []
    for pair in unique_chain_combs:
        # cycling over each coordinate of the first chain
        for s in range(coord_arrays[pair[0]].shape[0]):
            s_xyz = coord_arrays[pair[0]][s].reshape(1, 3)
            s_cid = coord_ids[pair[0]][s]
            dist = cdist(s_xyz, coord_arrays[pair[1]])
            npw = np.where(dist < cutoff)
            del dist
            for k in range(npw[0].shape[0]):
                con = (pair[0], s_cid, pair[1], coord_ids[pair[1]][npw[1][k]])
                con_list.append(con)
    return set(con_list)


class CAPRI:
    """CAPRI class."""

    def __init__(
        self,
        identificator: int,
        model: PDBPath,
        path: Path,
        reference: PDBPath,
        params: ParamMap,
        ref_id: int = 1,
    ) -> None:
        """
        Initialize the class.

        Parameters
        ----------
        identificator : int
            The identificator of the object.
        model : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`
            The model to be evaluated.
        path : Path
            Reference that defines where output should be saved.
        reference : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`
            The reference structure.
        params : dict
            The parameters for the CAPRI evaluation.
        """
        self.reference = reference
        if not isinstance(model, PDBFile):
            self.model = PDBFile(model)
            self.md5 = ""
            self.score = float("nan")
        else:
            self.model = model
            self.md5 = model.md5
            self.score = model.score
        self.path = path
        self.params = params
        self.irmsd = float("nan")
        self.lrmsd = float("nan")
        self.ilrmsd = float("nan")
        self.fnat = float("nan")
        self.dockq = float("nan")
        self.rmsd = float("nan")
        self.allatoms = params["allatoms"]
        self.atoms = self._load_atoms(model, reference, full=self.allatoms)
        self.r_chain = params["receptor_chain"]
        self.l_chains = params["ligand_chains"]
        self.model2ref_numbering = None
        self.model2ref_chain_dict = None
        self.output_ss_fname = Path(f"capri_ss_{identificator}.tsv")
        self.output_clt_fname = Path(f"capri_clt_{identificator}.tsv")
        self.output = self.output_ss_fname
        self.identificator = identificator
        self.core_model_idx = identificator
        self.ref_id = ref_id

    def calc_irmsd(self, cutoff: float = 5.0) -> None:
        """Calculate the I-RMSD.

        Parameters
        ----------
        cutoff : float
            The cutoff distance for the intermolecular contacts.
        """
        # Identify reference interface
        ref_interface_resdic = self.identify_interface(self.reference, cutoff)

        if len(ref_interface_resdic) == 0:
            log.warning("No reference interface found")
        else:
            # Load interface coordinates
            ref_coord_dic, _ = load_coords(
                self.reference, self.atoms, ref_interface_resdic
            )
            try:
                mod_coord_dic, _ = load_coords(
                    self.model,
                    self.atoms,
                    ref_interface_resdic,
                    numbering_dic=self.model2ref_numbering,
                    model2ref_chain_dict=self.model2ref_chain_dict,
                )
            except ALIGNError as alignerror:
                log.warning(alignerror)
                return

            # Here _coord_dic keys are matched
            #  and formatted as (chain, resnum, atom)
            #  we will use atoms that are present in both
            P = []
            Q = []

            for k in ref_coord_dic.keys() & mod_coord_dic.keys():
                ref_xyz = ref_coord_dic[k]
                mod_xyz = mod_coord_dic[k]

                Q.append(ref_xyz)
                P.append(mod_xyz)

            Q = np.asarray(Q)
            P = np.asarray(P)
            # write_coords("model.pdb", P)
            # write_coords("ref.pdb", Q)

            Q = Q - centroid(Q)
            P = P - centroid(P)
            U = kabsch(P, Q)
            P = np.dot(P, U)

            self.irmsd = calc_rmsd(P, Q)
            # write_coords("model_aln.pdb", P)
            # write_coords("ref_aln.pdb", Q)

    def calc_lrmsd(self) -> None:
        """Calculate the L-RMSD."""
        ref_coord_dic, _ = load_coords(self.reference, self.atoms)
        try:
            mod_coord_dic, _ = load_coords(
                self.model,
                self.atoms,
                numbering_dic=self.model2ref_numbering,
                model2ref_chain_dict=self.model2ref_chain_dict,
            )
        except ALIGNError as alignerror:
            log.warning(alignerror)
            return

        Q = []
        P = []
        # Note: this MUST be sorted since we will use the indexes to
        #  separate between receptor and ligand coordinates
        intersection = sorted(ref_coord_dic.keys() & mod_coord_dic.keys())

        chain_ranges: dict[Any, Any] = {}
        for i, segment in enumerate(intersection):
            chain, _, _ = segment
            if chain not in chain_ranges:
                chain_ranges[chain] = []
            chain_ranges[chain].append(i)

        chain_ranges = make_range(chain_ranges)
        obs_chains = list(chain_ranges.keys())  # observed chains
        if len(obs_chains) < 2:
            log.warning("Not enough chains for calculating lrmsd")
        else:
            r_chain, l_chains = check_chains(obs_chains, self.r_chain, self.l_chains)
            r_start, r_end = chain_ranges[r_chain]
            l_starts = [chain_ranges[_l][0] for _l in l_chains]
            l_ends = [chain_ranges[_l][1] for _l in l_chains]

            for k in intersection:
                ref_xyz = ref_coord_dic[k]
                mod_xyz = mod_coord_dic[k]

                Q.append(ref_xyz)
                P.append(mod_xyz)

            Q = np.asarray(Q)
            P = np.asarray(P)

            # write_coords("ref_first.pdb", Q)
            # write_coords("model_first.pdb", P)

            # get receptor and ligand coordinates
            Q_r_first = Q[r_start : r_end + 1]
            P_r_first = P[r_start : r_end + 1]
            # write_coords("ref_r_first.pdb", Q_r_first)
            # write_coords("model_r_first.pdb", P_r_first)
            # Q_l_first = Q[l_start: l_end + 1]
            # P_l_first = P[l_start: l_end + 1]
            # write_coords("ref_l_first.pdb", Q_l_first)
            # write_coords("model_l_first.pdb", P_l_first)

            # move to the origin of the receptor

            Q = Q - centroid(Q_r_first)
            P = P - centroid(P_r_first)

            # get receptor coordinates
            Q_r = Q[r_start : r_end + 1]
            P_r = P[r_start : r_end + 1]
            # Center receptors and get rotation matrix
            # Q_r = Q_r - centroid(Q_r)
            # P_r = P_r - centroid(P_r)
            # write_coords("ref_r_centr.pdb", Q_r)
            # write_coords("model_r_centr.pdb", P_r)

            U_r = kabsch(P_r, Q_r)

            # Apply rotation to complex
            #  - complex are now aligned by the receptor
            P = np.dot(P, U_r)

            # write_coords("ref.pdb", Q)
            # write_coords("model.pdb", P)

            # Identify ligand coordinates concatenating all the ligand chains
            Q_l = np.empty((0, 3))
            P_l = np.empty((0, 3))
            for l_start, l_end in zip(l_starts, l_ends):
                Q_l = np.concatenate((Q_l, Q[l_start : l_end + 1]))
                P_l = np.concatenate((P_l, P[l_start : l_end + 1]))
            # Q_l = Q[l_start: l_end + 1]
            # P_l = P[l_start: l_end + 1]

            # write_coords("ref_l.pdb", Q_l)
            # write_coords("model_l.pdb", P_l)

            # Calculate the RMSD of the ligands
            self.lrmsd = calc_rmsd(P_l, Q_l)

    def calc_ilrmsd(self, cutoff: float = 10.0) -> None:
        """Calculate the Interface Ligand RMSD.

        Parameters
        ----------
        cutoff : float
            The cutoff distance for the intermolecular contacts.
        """
        # Identify interface
        ref_interface_resdic = self.identify_interface(self.reference, cutoff)
        # Load interface coordinates

        ref_int_coord_dic, _ = load_coords(
            self.reference, self.atoms, ref_interface_resdic
        )
        try:
            mod_int_coord_dic, _ = load_coords(
                self.model,
                self.atoms,
                ref_interface_resdic,
                numbering_dic=self.model2ref_numbering,
                model2ref_chain_dict=self.model2ref_chain_dict,
            )
        except ALIGNError as alignerror:
            log.warning(alignerror)
            return

        # write_coord_dic("ref.pdb", ref_int_coord_dic)
        # write_coord_dic("model.pdb", mod_int_coord_dic)

        # find atoms present in both interfaces
        Q_int = []
        P_int = []
        common_keys = ref_int_coord_dic.keys() & mod_int_coord_dic.keys()
        for k in sorted(common_keys):
            ref_xyz = ref_int_coord_dic[k]
            mod_xyz = mod_int_coord_dic[k]

            Q_int.append(ref_xyz)
            P_int.append(mod_xyz)

        Q_int = np.asarray(Q_int)
        P_int = np.asarray(P_int)

        # write_coords("ref.pdb", Q_int)
        # write_coords("model.pdb", P_int)

        chain_ranges: dict[Any, Any] = {}
        for i, segment in enumerate(sorted(common_keys)):
            chain, _, _ = segment
            if chain not in chain_ranges:
                chain_ranges[chain] = []
            chain_ranges[chain].append(i)

        chain_ranges = make_range(chain_ranges)
        obs_chains = list(chain_ranges.keys())  # observed chains
        if len(obs_chains) < 2:
            log.warning("Not enough chains for calculating ilrmsd")
        else:
            r_chain, l_chains = check_chains(obs_chains, self.r_chain, self.l_chains)
            r_start, r_end = chain_ranges[r_chain]
            l_starts = [chain_ranges[l_chain][0] for l_chain in l_chains]
            l_ends = [chain_ranges[l_chain][1] for l_chain in l_chains]

            # write_coords("ref.pdb", Q)
            # write_coords("model.pdb", P)

            # put system at origin of the receptor interface
            Q_r_int = Q_int[r_start : r_end + 1]
            P_r_int = P_int[r_start : r_end + 1]

            Q_int = Q_int - centroid(Q_r_int)
            P_int = P_int - centroid(P_r_int)
            # put interfaces at the origin

            # find the rotation that minimizes the receptor interface rmsd
            Q_r_int = Q_int[r_start : r_end + 1]
            P_r_int = P_int[r_start : r_end + 1]

            U_int = kabsch(P_r_int, Q_r_int)
            P_int = np.dot(P_int, U_int)
            # just for checks.
            # the interface rmsd for the rec interface should be almost zero
            # Q_r_int = Q_int[r_start: r_end + 1]
            # P_r_int = P_int[r_start: r_end + 1]
            # r_rmsd = calc_rmsd(Q_r_int, P_int[r_start: r_end + 1])
            # print(r_rmsd)
            # Identify ligand coordinates concatenating all the ligand chains
            Q_l_int = np.empty((0, 3))
            P_l_int = np.empty((0, 3))
            for l_start, l_end in zip(l_starts, l_ends):
                Q_l_int = np.concatenate((Q_l_int, Q_int[l_start : l_end + 1]))
                P_l_int = np.concatenate((P_l_int, P_int[l_start : l_end + 1]))
            # prior to multibody:
            # Q_l_int = Q_int[l_start: l_end + 1]
            # P_l_int = P_int[l_start: l_end + 1]
            # write_coords("ref_l_int_fin.pdb", Q_l_int)
            # write_coords("mod_l_int_fin.pdb", P_l_int)

            # # this will be the interface-ligand-rmsd
            self.ilrmsd = calc_rmsd(P_l_int, Q_l_int)

    def calc_fnat(self, cutoff: float = 5.0) -> None:
        """Calculate the frequency of native contacts.

        Parameters
        ----------
        cutoff : float
            The cutoff distance for the intermolecular contacts.
        """
        ref_contacts = load_contacts(self.reference, cutoff)
        if len(ref_contacts) != 0:
            try:
                model_contacts = load_contacts(
                    self.model,
                    cutoff,
                    numbering_dic=self.model2ref_numbering,  # type: ignore
                    model2ref_chain_dict=self.model2ref_chain_dict,  # type: ignore
                )
            except ALIGNError as alignerror:
                log.warning(alignerror)
            else:
                intersection = ref_contacts & model_contacts
                self.fnat = len(intersection) / float(len(ref_contacts))
        else:
            log.warning("No reference contacts found")

    def calc_global_rmsd(self) -> None:
        """Calculate the full structure RMSD."""
        # Load reference atomic coordinates
        ref_coord_dic, _ = load_coords(self.reference, self.atoms)
        # Load model atomic coordinates
        try:
            model_coord_dic, _ = load_coords(
                self.model,
                self.atoms,
                numbering_dic=self.model2ref_numbering,
                model2ref_chain_dict=self.model2ref_chain_dict,
            )
        except ALIGNError as alignerror:
            log.warning(alignerror)
            return
        # Obtain list of coordinates
        Q = []
        P = []
        for k in ref_coord_dic.keys() & model_coord_dic.keys():
            ref_xyz = ref_coord_dic[k]
            mod_xyz = model_coord_dic[k]
            Q.append(ref_xyz)
            P.append(mod_xyz)
        # Cast indo array
        Q = np.asarray(Q)
        P = np.asarray(P)
        # Center to 0
        Q = Q - centroid(Q)
        P = P - centroid(P)
        # Obtain rotation matrix
        U = kabsch(P, Q)
        # Rotate model (the actual superimposition)
        P = np.dot(P, U)
        # Compute full RMSD
        self.rmsd = calc_rmsd(P, Q)

    def calc_dockq(self) -> None:
        """Calculate the DockQ metric."""
        self.dockq = 0.0
        if self.fnat:
            self.dockq += float(self.fnat) / 3
        if self.irmsd:
            irmsd_denom = 1 + (self.irmsd / 1.5) * (self.irmsd / 1.5)
            self.dockq += (1 / irmsd_denom) / 3
        if self.lrmsd:
            lrmsd_denom = 1 + (self.lrmsd / 8.5) * (self.lrmsd / 8.5)
            self.dockq += (1 / lrmsd_denom) / 3

    def has_cluster_info(self) -> bool:
        """
        Check wether this object contains cluster information.

        Returns
        -------
        bool
            True if this object contains cluster information.
        """
        has_cluster_info = False
        if self.model.clt_id:
            has_cluster_info = True
        return has_cluster_info

    def run(self) -> Union[None, "CAPRI"]:
        """Get the CAPRI metrics."""
        try:
            align_func = get_align(
                method=self.params["alignment_method"],
                lovoalign_exec=self.params["lovoalign_exec"],
            )
            self.model2ref_numbering, self.model2ref_chain_dict = align_func(
                self.reference, self.model, self.path
            )
        except ALIGNError:
            log.warning(
                f"Alignment failed between {self.reference} "
                f"and {self.model}, skipping..."
            )
            return
        # print(f"model2ref_numbering {self.model2ref_numbering}")
        # print(f"model2ref_chain_dict {self.model2ref_chain_dict}")
        if self.params["fnat"]:
            fnat_cutoff = self.params["fnat_cutoff"]
            self.calc_fnat(cutoff=fnat_cutoff)

        if self.params["irmsd"]:
            irmsd_cutoff = self.params["irmsd_cutoff"]
            self.calc_irmsd(cutoff=irmsd_cutoff)

        if self.params["lrmsd"]:
            self.calc_lrmsd()

        if self.params["ilrmsd"]:
            ilrmsd_cutoff = self.params["irmsd_cutoff"]
            self.calc_ilrmsd(cutoff=ilrmsd_cutoff)

        if self.params["dockq"]:
            self.calc_dockq()

        if self.params["global_rmsd"]:
            self.calc_global_rmsd()

        # The scheduler will use the return of the `run` method as the output of the tasks
        return copy.deepcopy(self)

    def __eq__(self, other):
        if self.params["dockq"] and \
                not (isnan(self.dockq) or isnan(other.dockq)):
            return self.dockq == other.dockq
        elif self.params["fnat"] and \
                not (isnan(self.fnat) or isnan(other.fnat)):
            return self.fnat == other.fnat
        elif self.params["ilrmsd"] and \
                not (isnan(self.ilrmsd) or isnan(other.ilrmsd)):
            return self.ilrmsd == other.ilrmsd
        elif self.params["lrmsd"] and \
                not (isnan(self.lrmsd) or isnan(other.lrmsd)):
            return self.lrmsd == other.lrmsd
        elif self.params["irmsd"] and \
                not (isnan(self.irmsd) or isnan(other.irmsd)):
            return self.irmsd == other.irmsd
        elif self.params["global_rmsd"] and \
                not (isnan(self.rmsd) or isnan(other.rmsd)):
            return self.rmsd == other.rmsd
        return True

    def __lt__(self, other):
        if self.params["dockq"] and \
                not (isnan(self.dockq) or isnan(other.dockq)):
            return self.dockq > other.dockq
        elif self.params["fnat"] and \
                not (isnan(self.fnat) or isnan(other.fnat)):
            return self.fnat > other.fnat
        elif self.params["ilrmsd"] and \
                not (isnan(self.ilrmsd) or isnan(other.ilrmsd)):
            return self.ilrmsd < other.ilrmsd
        elif self.params["lrmsd"] and \
                not (isnan(self.lrmsd) or isnan(other.lrmsd)):
            return self.lrmsd < other.lrmsd
        elif self.params["irmsd"] and \
                not (isnan(self.irmsd) or isnan(other.irmsd)):
            return self.irmsd < other.irmsd
        elif self.params["global_rmsd"] and \
                not (isnan(self.rmsd) or isnan(other.rmsd)):
            return self.rmsd < other.rmsd
        return False

    @staticmethod
    def _load_atoms(
        model: PDBPath,
        reference: PDBPath,
        full: bool = False,
    ) -> AtomsDict:
        """
        Load atoms from a model and reference.

        Parameters
        ----------
        model : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`
            PDB file of the model to have its atoms identified
        reference : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`
            PDB file of the model to have its atoms identified
        full : bool
            If False, only backbone atoms will be retrieved, otherwise all atoms

        Returns
        -------
        atom_dic : dict
            Dictionary containing atoms observed in model and reference
        """
        model_atoms = get_atoms(model, full=full)
        reference_atoms = get_atoms(reference, full=full)
        atoms_dict: AtomsDict = {}
        atoms_dict.update(model_atoms)
        atoms_dict.update(reference_atoms)
        return atoms_dict

    @staticmethod
    def identify_interface(
        pdb_f: PDBPath,
        cutoff: float = 5.0,
    ) -> dict[str, list[int]]:
        """Identify the interface.

        Parameters
        ----------
        pdb_f : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`
            PDB file of the model to have its atoms identified
        cutoff : float, optional
            Cutoff distance for the interface identification.

        Returns
        -------
        interface_resdic : dict[str, list[int]]
            Dictionary holding list of interface residues ids for each chains.
        """
        if isinstance(pdb_f, PDBFile):
            pdb_f = pdb_f.rel_path

        interface_resdic: dict[str, list[int]] = {}
        contacts = load_contacts(pdb_f, cutoff)

        for contact in contacts:
            first_chain, first_resid, sec_chain, sec_resid = contact

            if first_chain not in interface_resdic:
                interface_resdic[first_chain] = []
            if sec_chain not in interface_resdic:
                interface_resdic[sec_chain] = []

            if first_resid not in interface_resdic[first_chain]:
                interface_resdic[first_chain].append(first_resid)
            if sec_resid not in interface_resdic[sec_chain]:
                interface_resdic[sec_chain].append(sec_resid)

        return interface_resdic

    @staticmethod
    def add_chain_from_segid(pdb_path: PDBPath) -> Path:
        """
        Replace the chainID with the segID.

        Parameters
        ----------
        pdb_path : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`
            PDB file to be replaced
        """
        if isinstance(pdb_path, PDBFile):
            pdb_path = pdb_path.rel_path
        temp_f = tempfile.NamedTemporaryFile(delete=False, mode="w+t")
        with open(pdb_path) as fh:
            for line in list(pdb_segxchain.run(fh)):
                temp_f.writelines(line)
        temp_f.close()
        # REPLACE!
        new_pdb_path = shutil.move(temp_f.name, pdb_path)
        return new_pdb_path


def rank_according_to_score(
    data: dict[int, ParamDict], sort_key: str, sort_ascending: bool
) -> dict[int, ParamDict]:
    """
    Ranks a dictionary of data based on a specified sort key and sort order,
    and assigns a rank to each entry based on its 'score' attribute.

    Args:
        data (dict[int, ParamDict]): Dictionary where each key is an index and each
                                     value is a ParamDict containing data attributes.
        sort_key (str): Key by which to sort the data within the ParamDict.
                        Must correspond to a valid attribute in ParamDict.
        sort_ascending (bool): If True, sorts the data in ascending order based on
                               the sort_key; if False, sorts in descending order.

    Returns:
        dict[int, ParamDict]: A new dictionary where entries are sorted according
                              to the sort_key and optionally sorted order. Each entry
                              also includes a 'caprieval_rank' attribute indicating
                              its rank based on the 'score'.
    """
    score_rankkey_values = [(k, v["score"]) for k, v in data.items()]
    score_rankkey_values.sort(key=lambda x: x[1])

    for i, k in enumerate(score_rankkey_values):
        data_idx, _ = k
        data[data_idx]["caprieval_rank"] = i + 1

    # Sort according to the sort key
    rankkey_values = [(k, v[sort_key]) for k, v in data.items()]
    rankkey_values.sort(
        key=lambda x: x[1],
        reverse=True if not sort_ascending else False,
    )

    _data = {}
    for i, (data_idx, _) in enumerate(rankkey_values):
        _data[i + 1] = data[data_idx]
    data = _data

    return _data


def extract_models_best_references(capri_objects: list[CAPRI]) -> list[CAPRI]:
    """Extract best reference for each input model.

    Same input models are combined and best reference is later found by
    sorting the CAPRI objects. Only the best performing CAPRI object is
    kept and returned.
    This step was implemented to handle comparisons against multiple refs.

    Parameters
    ----------
    capri_objects : list[CAPRI]
        List of CAPRI object.

    Returns
    -------
    selected_capri_objects : list[CAPRI]
        List of selected best CAPIR object for each model.
    """
    # Group results by models
    by_model_data: dict[int, dict[Path, CAPRI]] = {}
    for capri in capri_objects:
        model_data = by_model_data.setdefault(capri.identificator, [])
        model_data.append(capri)
    # Finds best performances for each model
    selected_capri_objects: list[CAPRI] = []
    # Loop over model referneces performances
    for references_perfs in by_model_data.values():
        # Sort them using built in __eq__ and __lt__ CAPRI methods
        sorted_perfs = sorted(references_perfs)
        # Select first on (best)
        best_perf = sorted_perfs[0]
        # Hold that guy
        selected_capri_objects.append(best_perf)
    return selected_capri_objects


def extract_data_from_capri_class(
    capri_objects: list[CAPRI],
    sort_key: str,
    sort_ascending: bool,
    output_fname: Path,
    add_reference_id: bool = False,
) -> Union[dict[int, ParamDict], None]:
    """Extracts data attributes from a list of CAPRI objects into a structured
    dictionary, optionally sorts the data based on a specified key, and writes
    the sorted data to a file.

    Parameters
    ----------
    capri_objects : list[CAPRI]
        List of CAPRI objects containing data attributes to be extracted.
    sort_key : str
        Key by which to sort the extracted data. Must correspond to a valid
        attribute in the CAPRI object (e.g., 'score', 'irmsd').
    sort_ascending : bool
        If True, sorts the data in ascending order based on the sort_key;
        if False, sorts in descending order.
    output_fname : Path
        Path to the output file where the sorted data will be written.
    add_reference_id : bool, optional
        Should the reference id be added to the capri table?, by default False

    Returns
    -------
    ranked_data : Union[dict[int, ParamDict], None]
        The sorted and structured data dictionary if successful,
        None if no data was processed.
    """
    # Retrieve data for each model
    data: dict[int, ParamDict] = {}
    for i, c in enumerate(capri_objects, start=1):
        data[i] = {
            "model": c.model,
            "md5": c.md5,
            "caprieval_rank": None,
            "score": c.score,
            "irmsd": c.irmsd,
            "fnat": c.fnat,
            "lrmsd": c.lrmsd,
            "ilrmsd": c.ilrmsd,
            "dockq": c.dockq,
            "rmsd": c.rmsd,
            "cluster_id": c.model.clt_id if c.model.clt_id else None,
            "cluster_ranking": c.model.clt_rank if c.model.clt_rank else None,
            "model-cluster_ranking": (
                c.model.clt_model_rank if c.model.clt_model_rank else None
            ),
        }
        if c.model.unw_energies is not None:
            data[i].update(c.model.unw_energies)
        if add_reference_id:
            data[i]["ref_id"] = c.ref_id
    # Sort models
    ranked_data = rank_according_to_score(
        data, sort_key=sort_key, sort_ascending=sort_ascending
    )
    if not ranked_data:
        # This means no files have been collected
        return
    else:
        write_nested_dic_to_file(data_dict=ranked_data, output_fname=output_fname)
        return ranked_data


def calc_stats(data: list) -> tuple[float, float]:
    """
    Calculate the mean and stdev.

    Parameters
    ----------
    data : list
        List of values.

    Returns
    -------
    mean : float
        Mean of the values.
    stdev : float
        Standard deviation of the values.
    """
    mean = np.mean(data)
    stdev = np.std(data)
    return mean, stdev


# Define dict types
CltData = dict[tuple[Optional[int], Union[int, str, None]], list[tuple[CAPRI, PDBFile]]]


def capri_cluster_analysis(
    capri_list: Iterable[CAPRI],
    model_list: Iterable[PDBFile],
    output_fname: FilePath,
    clt_threshold: int,
    sort_key: str,
    sort_ascending: bool,
    path: FilePath,
) -> None:
    """Consider the cluster results for the CAPRI evaluation."""
    capri_keys = ["irmsd", "fnat", "lrmsd", "dockq", "ilrmsd", "rmsd"]
    model_keys = ["air", "bsa", "desolv", "elec", "total", "vdw"]
    log.info(f"Rearranging cluster information into {output_fname}")
    # get the cluster data
    clt_data: CltData = dict(((m.clt_rank, m.clt_id), []) for m in model_list)

    # add models to each cluster
    for capri, model in zip(capri_list, model_list):
        clt_data[(model.clt_rank, model.clt_id)].append((capri, model))

    output_dic: dict[int, ParamDict] = {}

    for i, element in enumerate(clt_data):
        data: ParamDict = {}
        number_of_models_in_cluster = len(clt_data[element])
        # rank, cluster id, number of models in cluster
        data["cluster_rank"] = element[0]
        data["cluster_id"] = element[1]
        data["n"] = number_of_models_in_cluster
        if number_of_models_in_cluster < clt_threshold:
            # under-evaluated, the mean was divided by a value
            #  larger than the total number of models in the cluster
            data["under_eval"] = "yes"
        else:
            data["under_eval"] = "-"
        # score
        try:
            score_array = [e[1].score for e in clt_data[element][:clt_threshold]]
            data["score"], data["score_std"] = calc_stats(score_array)
        except KeyError:
            data["score"] = float("nan")
            data["score_std"] = float("nan")

        # capri keys
        for key in capri_keys:
            std_key = f"{key}_std"
            try:
                key_array = [vars(e[0])[key] for e in clt_data[element][:clt_threshold]]
                data[key], data[std_key] = calc_stats(key_array)
            except KeyError:
                data[key] = float("nan")
                data[std_key] = float("nan")

        # model keys
        for key in model_keys:
            std_key = f"{key}_std"
            if clt_data[element][0][1].unw_energies:
                try:
                    key_array = [
                        vars(e[1])["unw_energies"][key]
                        for e in clt_data[element][:clt_threshold]
                    ]
                    data[key], data[std_key] = calc_stats(key_array)
                except KeyError:
                    data[key] = float("nan")
                    data[std_key] = float("nan")

        output_dic[i] = data

    # Rank according to the score
    score_rankkey_values = [(key, v["score"]) for key, v in output_dic.items()]
    score_rankkey_values.sort(key=lambda x: x[1])
    for i, k in enumerate(score_rankkey_values):
        idx, _ = k
        output_dic[idx]["caprieval_rank"] = i + 1

    # Rank according to the sorting key
    rankkey_values = [(key, v[sort_key]) for key, v in output_dic.items()]
    rankkey_values.sort(
        key=lambda x: x[1],
        reverse=True if not sort_ascending else False,
    )

    _output_dic = {}
    for i, k in enumerate(rankkey_values):
        idx, _ = k
        _output_dic[i + 1] = output_dic[idx]
    output_dic = _output_dic

    output_fname = Path(path, output_fname)

    info_header = "#" * 40 + os.linesep
    info_header += "# `caprieval` cluster-based analysis" + os.linesep
    info_header += "#" + os.linesep
    info_header += f"# > sortby_key={sort_key}" + os.linesep
    info_header += f"# > sort_ascending={sort_ascending}" + os.linesep
    info_header += f"# > clt_threshold={clt_threshold}" + os.linesep
    info_header += "#" + os.linesep
    info_header += (
        "# NOTE: if under_eval=yes, it means that there were less models in"
        " a cluster than" + os.linesep
    )
    info_header += (
        "#    clt_threshold, thus these values were under " "evaluated." + os.linesep
    )
    info_header += (
        "#   You might need to tweak the value of clt_threshold or change"
        " some parameters" + os.linesep
    )
    info_header += "#    in `clustfcc` depending on your " "analysis." + os.linesep
    info_header += "#" + os.linesep
    info_header += "#" * 40

    if not data:
        # This means there were only "dummy" values
        return
    else:
        write_nested_dic_to_file(
            output_dic,
            output_fname,
            info_header=info_header,
        )


class CAPRIError(Exception):
    """Raised when something goes wrong with the CAPRI class."""

    def __init__(self, msg: str = "") -> None:
        self.msg = msg
        super().__init__(self.msg)


def dump_weights(order: int) -> None:
    sel_steps = get_module_steps_folders(Path(".."))
    cns_step = get_previous_cns_step(sel_steps=sel_steps, st_order=order)
    if cns_step:
        log.info(f"Found previous CNS step: {cns_step}")
        scoring_params_fname = save_scoring_weights(cns_step)
        log.info(f"Saved scoring weights to: {scoring_params_fname}")
    else:
        log.info("No previous CNS step found. Cannot save scoring weights.")


# # debug only
# def write_coord_dic(output_name, coord_dic):
#     """Add a dummy atom to a PDB file according to a list of coordinates."""
#     with open(output_name, "w") as fh:
#         for i, k in enumerate(coord_dic):
#             atom_num = f"{i+1}".rjust(4, " ")
#             chain, resnum, atom = k
#             resnum = int(resnum)
#             resnum = f"{resnum}".rjust(3, " ")
#             atom_name = f"{atom}".rjust(3, " ")
#             x, y, z = coord_dic[k]
#             dum_x = f"{x:.3f}".rjust(7, " ")
#             dum_y = f"{y:.3f}".rjust(7, " ")
#             dum_z = f"{z:.3f}".rjust(7, " ")
#             dummy_line = (
#                 f"ATOM   {atom_num} {atom_name}  DUM {chain} {resnum}   "
#                 f"  {dum_x} {dum_y} {dum_z}  1.00  1.00   "
#                 "        H  " + os.linesep
#                 )
#             fh.write(dummy_line)


# # debug only
# def write_coords(output_name, coor_list):
#     """Add a dummy atom to a PDB file according to a list of coordinates."""
#     with open(output_name, "w") as fh:
#         for i, dummy_coord in enumerate(coor_list):
#             atom_num = f"{i}".rjust(4, " ")
#             resnum = f"{i}".rjust(3, " ")
#             dum_x = f"{dummy_coord[0]:.3f}".rjust(7, " ")
#             dum_y = f"{dummy_coord[1]:.3f}".rjust(7, " ")
#             dum_z = f"{dummy_coord[2]:.3f}".rjust(7, " ")
#             dummy_line = (
#                 f"ATOM   {atom_num}  H   DUM X {resnum}   "
#                 f"  {dum_x} {dum_y} {dum_z}  1.00  1.00   "
#                 "        H  " + os.linesep
#                 )
#             fh.write(dummy_line)


# # debug only
# def write_pymol_viz(resdic):
#     """Write PyMol vizualitation."""
#     for k in resdic:
#         reslist = "+".join(map(str, resdic[k]))
#         cmd = f"sele {k}, chain {k} and resid {reslist}"
#         print(cmd)
