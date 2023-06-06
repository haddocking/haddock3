"""CAPRI module."""
import os
import shutil
import tempfile
from itertools import combinations
from pathlib import Path

import numpy as np
from pdbtools import pdb_segxchain
from scipy.spatial.distance import cdist

from haddock import log
from haddock.libs.libalign import (
    AlignError,
    calc_rmsd,
    centroid,
    get_align,
    get_atoms,
    kabsch,
    load_coords,
    make_range,
    )
from haddock.libs.libio import write_dic_to_file, write_nested_dic_to_file
from haddock.libs.libontology import PDBFile


def load_contacts(pdb_f, cutoff=5.0):
    """
    Load residue-based contacts.

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
    con_list = []
    if isinstance(pdb_f, PDBFile):
        pdb_f = pdb_f.rel_path
    # get also side chains atoms
    atoms = get_atoms(pdb_f, full=True)
    ref_coord_dic, _ = load_coords(pdb_f, atoms)
    # create coordinate arrays
    coord_arrays, coord_ids = {}, {}
    for atom in ref_coord_dic.keys():
        chain = atom[0]
        if chain not in coord_arrays.keys():  # initialize lists
            coord_arrays[chain], coord_ids[chain] = [], []
        coord_arrays[chain].append(ref_coord_dic[atom])
        coord_ids[chain].append(atom[1])  # only the resid is appended
    for chain in coord_arrays.keys():
        coord_arrays[chain] = np.array(coord_arrays[chain])

    # combinations of chains
    unique_chain_combs = list(combinations(coord_arrays.keys(), 2))

    # calculating contacts
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
            identificator,
            model,
            path,
            reference,
            params,
            ):
        """
        Initialize the class.

        Parameters
        ----------
        identificator : str
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
        self.model = model
        self.path = path
        self.params = params
        self.irmsd = float('nan')
        self.lrmsd = float('nan')
        self.ilrmsd = float('nan')
        self.fnat = float('nan')
        self.dockq = float('nan')
        self.atoms = self._load_atoms(model, reference)
        self.r_chain = params["receptor_chain"]
        self.l_chain = params["ligand_chain"]
        self.model2ref_numbering = None
        self.output_ss_fname = Path(f"capri_ss_{identificator}.tsv")
        self.output_clt_fname = Path(f"capri_clt_{identificator}.tsv")
        # for parallelisation
        self.output = self.output_ss_fname
        self.identificator = identificator
        self.core_model_idx = identificator

    def calc_irmsd(self, cutoff=5.0):
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

            mod_coord_dic, _ = load_coords(
                self.model,
                self.atoms,
                ref_interface_resdic,
                numbering_dic=self.model2ref_numbering
                )

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

    def calc_lrmsd(self):
        """Calculate the L-RMSD."""
        ref_coord_dic, _ = load_coords(self.reference, self.atoms)

        mod_coord_dic, _ = load_coords(
            self.model,
            self.atoms,
            numbering_dic=self.model2ref_numbering
            )

        Q = []
        P = []
        # Note: this MUST be sorted since we will use the indexes to
        #  separate between receptor and ligand coordinates
        intersection = sorted(ref_coord_dic.keys() & mod_coord_dic.keys())

        chain_ranges = {}
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
            r_chain, l_chain = self.check_chains(obs_chains)
            r_start, r_end = chain_ranges[r_chain]
            l_start, l_end = chain_ranges[l_chain]

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
            Q_r_first = Q[r_start: r_end + 1]
            P_r_first = P[r_start: r_end + 1]
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
            Q_r = Q[r_start: r_end + 1]
            P_r = P[r_start: r_end + 1]
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

            # Identify the ligand coordinates
            Q_l = Q[l_start: l_end + 1]
            P_l = P[l_start: l_end + 1]

            # write_coords("ref_l.pdb", Q_l)
            # write_coords("model_l.pdb", P_l)

            # Calculate the RMSD of the ligands
            self.lrmsd = calc_rmsd(P_l, Q_l)

    def calc_ilrmsd(self, cutoff=10.0):
        """
        Calculate the Interface Ligand RMSD.

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

        mod_int_coord_dic, _ = load_coords(
            self.model,
            self.atoms,
            ref_interface_resdic,
            numbering_dic=self.model2ref_numbering
            )

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

        chain_ranges = {}
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
            r_chain, l_chain = self.check_chains(obs_chains)

            l_start, l_end = chain_ranges[l_chain]
            r_start, r_end = chain_ranges[r_chain]

            # write_coords("ref.pdb", Q)
            # write_coords("model.pdb", P)

            # put system at origin of the receptor interface
            Q_r_int = Q_int[r_start: r_end + 1]
            P_r_int = P_int[r_start: r_end + 1]

            Q_int = Q_int - centroid(Q_r_int)
            P_int = P_int - centroid(P_r_int)
            # put interfaces at the origin

            # find the rotation that minimizes the receptor interface rmsd
            Q_r_int = Q_int[r_start: r_end + 1]
            P_r_int = P_int[r_start: r_end + 1]

            U_int = kabsch(P_r_int, Q_r_int)
            P_int = np.dot(P_int, U_int)
            # just for checks.
            # the interface rmsd for the rec interface should be almost zero
            # Q_r_int = Q_int[r_start: r_end + 1]
            # P_r_int = P_int[r_start: r_end + 1]
            # r_rmsd = calc_rmsd(Q_r_int, P_int[r_start: r_end + 1])
            # print(r_rmsd)

            Q_l_int = Q_int[l_start: l_end + 1]
            P_l_int = P_int[l_start: l_end + 1]
            # write_coords("ref_l_int_fin.pdb", Q_l_int)
            # write_coords("mod_l_int_fin.pdb", P_l_int)

            # # this will be the interface-ligand-rmsd
            self.ilrmsd = calc_rmsd(P_l_int, Q_l_int)

    def calc_fnat(self, cutoff=5.0):
        """
        Calculate the frequency of native contacts.

        Parameters
        ----------
        cutoff : float
            The cutoff distance for the intermolecular contacts.
        """
        ref_contacts = load_contacts(self.reference, cutoff)
        if len(ref_contacts) != 0:
            model_contacts = load_contacts(self.model, cutoff)
            intersection = ref_contacts & model_contacts
            self.fnat = len(intersection) / float(len(ref_contacts))
        else:
            log.warning("No reference contacts found")

    def calc_dockq(self):
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

    def has_cluster_info(self):
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

    def make_output(self):
        """Output the CAPRI results to a .tsv file."""
        data = {}
        # keep always "model" the first key
        data["model"] = self.model
        data["md5"] = self.model.md5
        # create the empty rank here so that it will appear
        #  as the second column
        data["caprieval_rank"] = None
        data["score"] = self.model.score
        data["irmsd"] = self.irmsd
        data["fnat"] = self.fnat
        data["lrmsd"] = self.lrmsd
        data["ilrmsd"] = self.ilrmsd
        data["dockq"] = self.dockq

        if self.has_cluster_info():
            data["cluster-id"] = self.model.clt_id
            data["cluster-ranking"] = self.model.clt_rank
            data["model-cluster-ranking"] = self.model.clt_model_rank
        else:
            data["cluster-id"] = None
            data["cluster-ranking"] = None
            data["model-cluster-ranking"] = None

        # energies
        if self.model.unw_energies:
            for key in self.model.unw_energies:
                data[key] = self.model.unw_energies[key]

        output_fname = Path(self.path, self.output_ss_fname)

        write_dic_to_file(data, output_fname)

    def run(self):
        """Get the CAPRI metrics."""
        try:
            align_func = get_align(
                method=self.params["alignment_method"],
                lovoalign_exec=self.params["lovoalign_exec"]
                )
            self.model2ref_numbering = align_func(
                self.reference,
                self.model,
                self.path
                )
        except AlignError:
            log.warning(
                f"Alignment failed between {self.reference} "
                f"and {self.model}, skipping..."
                )
            return

        if self.params["fnat"]:
            log.debug(f"id {self.identificator}, calculating FNAT")
            fnat_cutoff = self.params["fnat_cutoff"]
            log.debug(f" cutoff: {fnat_cutoff}A")
            self.calc_fnat(cutoff=fnat_cutoff)

        if self.params["irmsd"]:
            log.debug(f"id {self.identificator}, calculating I-RMSD")
            irmsd_cutoff = self.params["irmsd_cutoff"]
            log.debug(f" cutoff: {irmsd_cutoff}A")
            self.calc_irmsd(cutoff=irmsd_cutoff)

        if self.params["lrmsd"]:
            log.debug(f"id {self.identificator}, calculating L-RMSD")
            self.calc_lrmsd()

        if self.params["ilrmsd"]:
            log.debug(f"id {self.identificator}, calculating I-L-RMSD")
            ilrmsd_cutoff = self.params["irmsd_cutoff"]
            log.debug(f" cutoff: {ilrmsd_cutoff}A")
            self.calc_ilrmsd(cutoff=ilrmsd_cutoff)

        if self.params["dockq"]:
            log.debug(f"id {self.identificator}, calculating DockQ metric")
            self.calc_dockq()

        self.make_output()

    def check_chains(self, obs_chains):
        """Check observed chains against the expected ones."""
        r_found, l_found = False, False
        obs_chains_cp = obs_chains.copy()
        if self.r_chain in obs_chains:
            r_chain = self.r_chain
            obs_chains.remove(r_chain)
            r_found = True
        if self.l_chain in obs_chains:
            l_chain = self.l_chain
            obs_chains.remove(l_chain)
            l_found = True
        # if one or both exp chains are not observed, use the observed chains
        if obs_chains != []:
            exps = {self.r_chain, self.l_chain}
            log.warning(f"observed chains != expected chains {exps}.")
            log.info(f"Sticking to observed chains {obs_chains_cp}")
        if not r_found:
            r_chain = obs_chains[0]
            obs_chains.remove(r_chain)
        if not l_found:
            l_chain = obs_chains[0]

        return r_chain, l_chain

    @staticmethod
    def _load_atoms(model, reference):
        """
        Load atoms from a model and reference.

        Parameters
        ----------
        model : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`
            PDB file of the model to have its atoms identified
        reference : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`
            PDB file of the model to have its atoms identified

        Returns
        -------
        atom_dic : dict
            Dictionary containing atoms observed in model and reference
        """
        model_atoms = get_atoms(model)
        reference_atoms = get_atoms(reference)
        atoms_dict = {}
        atoms_dict.update(model_atoms)
        atoms_dict.update(reference_atoms)
        return atoms_dict

    @staticmethod
    def identify_interface(pdb_f, cutoff=5.0):
        """Identify the interface.

        Parameters
        ----------
        pdb_f : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`
            PDB file of the model to have its atoms identified
        cutoff : float, optional
            Cutoff distance for the interface identification.
        """
        if isinstance(pdb_f, PDBFile):
            pdb_f = pdb_f.rel_path

        interface_resdic = {}
        contacts = load_contacts(pdb_f, cutoff)

        for contact in contacts:

            first_chain, first_resid = contact[0], contact[1]
            sec_chain, sec_resid = contact[2], contact[3]

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
    def add_chain_from_segid(pdb_path):
        """
        Replace the chainID with the segID.

        Parameters
        ----------
        pdb_path : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`
            PDB file to be replaced
        """
        temp_f = tempfile.NamedTemporaryFile(delete=False, mode="w+t")
        with open(pdb_path) as fh:
            for line in list(pdb_segxchain.run(fh)):
                temp_f.writelines(line)
        temp_f.close()
        # REPLACE!
        new_pdb_path = shutil.move(temp_f.name, pdb_path)
        return new_pdb_path


def merge_data(capri_jobs):
    """Merge CAPRI data."""
    capri_dic = {}
    for ident in range(1, len(capri_jobs) + 1):
        out_file = Path(f"capri_ss_{ident}.tsv")
        if not out_file.exists():
            continue
        header, content = out_file.read_text().split(os.linesep, 1)

        header_data = header.split('\t')
        content_data = content.split('\t')

        model_name = Path(content_data[header_data.index("model")]).name
        capri_dic[model_name] = {}
        target_keys = ['irmsd', 'fnat', 'ilrmsd', 'lrmsd', 'dockq']
        for key in target_keys:
            capri_dic[model_name][key] = float(
                content_data[header_data.index(key)])

    for j in capri_jobs:
        for m in capri_dic:
            if m == j.model.file_name:
                # add the data
                j.irmsd = capri_dic[m]['irmsd']
                j.fnat = capri_dic[m]['fnat']
                j.lrmsd = capri_dic[m]['lrmsd']
                j.ilrmsd = capri_dic[m]['ilrmsd']
                j.dockq = capri_dic[m]['dockq']

    return capri_jobs


def rearrange_ss_capri_output(
        output_name,
        output_count,
        sort_key,
        sort_ascending,
        path
        ):
    """
    Combine different capri outputs in a single file.

    Parameters
    ----------
    output_name : str
        Name of the output file.
    output_count : int
        Number of output files to combine.
    sort_key : str
        Key to sort the output files.
    path : Path
        Path to the output directory.
    """
    # this would be easier and more readable with pandas (:
    output_fname = Path(path, output_name)
    log.info(f"Rearranging output files into {output_fname}")
    keyword = output_name.split(".")[0]
    split_dict = {
        "capri_ss": "model-cluster-ranking",
        "capri_clt": "caprieval_rank"
        }
    if keyword not in split_dict.keys():
        raise Exception(f'Keyword {keyword} does not exist.')

    # Load the information of each intermediate file
    data = {}
    for ident in range(1, output_count + 1):
        out_file = Path(path, f"{keyword}_{ident}.tsv")

        # raise a warning if file does not exist.
        if not out_file.exists():
            log.warning((f"Output file {out_file} does not exist. "
                        "Caprieval will not be exhaustive..."))
            continue

        data[ident] = {}
        header, content = out_file.read_text().split(os.linesep, 1)

        header_data = header.split('\t')
        content_data = content.split('\t')

        # find out the data type of each field
        for key, value in zip(header_data, content_data):
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    value = str(value).strip(os.linesep)
            data[ident][key] = value

        out_file.unlink()

    # Rank according to the score
    score_rankkey_values = [(k, data[k]['score']) for k in data.keys()]
    score_rankkey_values.sort(key=lambda x: x[1])

    for i, k in enumerate(score_rankkey_values):
        data_idx, _ = k
        data[data_idx]["caprieval_rank"] = i + 1

    # Sort according to the sort key
    rankkey_values = [(k, data[k][sort_key]) for k in data.keys()]
    rankkey_values.sort(
        key=lambda x: x[1],
        reverse=True if not sort_ascending else False
        )

    _data = {}
    for i, (data_idx, _) in enumerate(rankkey_values):
        _data[i + 1] = data[data_idx]
    data = _data

    if not data:
        # This means no files have been collected
        return
    else:
        write_nested_dic_to_file(data, output_name)


def calc_stats(data):
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


def capri_cluster_analysis(
        capri_list,
        model_list,
        output_fname,
        clt_threshold,
        sort_key,
        sort_ascending,
        path
        ):
    """Consider the cluster results for the CAPRI evaluation."""
    capri_keys = ["irmsd", "fnat", "lrmsd", "dockq"]
    model_keys = ["air", "bsa", "desolv", "elec", "total", "vdw"]
    log.info(f"Rearranging cluster information into {output_fname}")
    # get the cluster data
    clt_data = dict(((m.clt_rank, m.clt_id), []) for m in model_list)

    # add models to each cluster
    for capri, model in zip(capri_list, model_list):
        clt_data[(model.clt_rank, model.clt_id)].append((capri, model))

    output_dic = {}

    for i, element in enumerate(clt_data):
        data = {}
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
            score_array = [
                e[1].score for e in clt_data[element][:clt_threshold]]
            data["score"], data["score_std"] = calc_stats(score_array)
        except KeyError:
            data["score"] = float("nan")
            data["score_std"] = float("nan")

        # capri keys
        for key in capri_keys:
            std_key = f"{key}_std"
            try:
                key_array = [
                    vars(e[0])[key] for e in clt_data[element][:clt_threshold]]
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
                        vars(e[1])["unw_energies"][key] for e in clt_data[element][:clt_threshold]] # noqa
                    data[key], data[std_key] = calc_stats(key_array)
                except KeyError:
                    data[key] = float("nan")
                    data[std_key] = float("nan")

        output_dic[i] = data

    # Rank according to the score
    score_rankkey_values = [(key, output_dic[key]['score'])
                            for key in output_dic.keys()]
    score_rankkey_values.sort(key=lambda x: x[1])
    for i, k in enumerate(score_rankkey_values):
        idx, _ = k
        output_dic[idx]["caprieval_rank"] = i + 1

    # Rank according to the sorting key
    rankkey_values = [(key, output_dic[key][sort_key])
                      for key in output_dic.keys()]
    rankkey_values.sort(
        key=lambda x: x[1],
        reverse=True if not sort_ascending else False
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
        "#    clt_threshold, thus these values were under evaluated."
        + os.linesep
        )
    info_header += (
        "#   You might need to tweak the value of clt_threshold or change"
        " some parameters" + os.linesep
        )
    info_header += (
        "#    in `clustfcc` depending on your analysis." + os.linesep
        )
    info_header += "#" + os.linesep
    info_header += "#" * 40

    if not data:
        # This means there were only "dummy" values
        return
    else:
        write_nested_dic_to_file(
            output_dic,
            output_fname,
            info_header=info_header)


class CAPRIError(Exception):
    """Raised when something goes wrong with the CAPRI class."""

    def __init__(self, msg=""):
        self.msg = msg
        super().__init__(self.msg)


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