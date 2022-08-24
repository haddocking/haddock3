"""CAPRI module."""
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
from fccpy import read_pdb
from fccpy.contacts import get_intermolecular_contacts
from pdbtools import pdb_segxchain

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

        # write_coord_dic("ref.pdb", ref_coord_dic)
        # write_coord_dic("model.pdb", mod_coord_dic)

        # write_coords("ref.pdb", Q)
        # write_coords("model.pdb", P)

        # move to the origin
        Q = Q - centroid(Q)
        P = P - centroid(P)

        # get receptor coordinates
        Q_r = Q[r_start: r_end - 1]
        P_r = P[r_start: r_end - 1]

        # Center receptors and get rotation matrix
        # Q_r = Q_r - centroid(Q_r)
        # P_r = P_r - centroid(P_r)

        U_r = kabsch(P_r, Q_r)

        # Center complexes at receptor centroids
        Q = Q - centroid(Q_r)
        P = P - centroid(P_r)

        # Apply rotation to complex
        #  - complex are now aligned by the receptor
        P = np.dot(P, U_r)

        # write_coords("ref.pdb", Q)
        # write_coords("model.pdb", P)

        # Identify the ligand coordinates
        Q_l = Q[l_start: l_end - 1]
        P_l = P[l_start: l_end - 1]

        # write_coords("ref_l.pdb", Q_l)
        # write_coords("model_l.pdb", P_l)

        # Calculate the RMSD of the ligands
        self.lrmsd = calc_rmsd(P_l, Q_l)

        # write_coords("ref.pdb", Q)
        # write_coords("model.pdb", P)

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
        ref_coord_dic, _ = load_coords(self.reference, self.atoms)

        ref_int_coord_dic, _ = load_coords(
            self.reference, self.atoms, ref_interface_resdic
            )

        mod_coord_dic, _ = load_coords(
            self.model,
            self.atoms,
            numbering_dic=self.model2ref_numbering
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

        # find atoms present in both molecules
        Q = []
        P = []
        intersection = sorted(ref_coord_dic.keys() & mod_coord_dic.keys())
        chain_ranges = {}
        for i, segment in enumerate(intersection):
            chain, _, _ = segment
            if chain not in chain_ranges:
                chain_ranges[chain] = []
            chain_ranges[chain].append(i)

        chain_ranges = make_range(chain_ranges)

        obs_chains = list(chain_ranges.keys())  # observed chains
        r_chain, l_chain = self.check_chains(obs_chains)

        l_start, l_end = chain_ranges[l_chain]

        for k in sorted(ref_coord_dic.keys() & mod_coord_dic.keys()):
            ref_xyz = ref_coord_dic[k]
            mod_xyz = mod_coord_dic[k]

            Q.append(ref_xyz)
            P.append(mod_xyz)

        Q = np.asarray(Q)
        P = np.asarray(P)

        # write_coords("ref.pdb", Q)
        # write_coords("model.pdb", P)

        # put system at origin
        Q_int = Q_int - centroid(Q_int)
        P_int = P_int - centroid(P_int)

        # put interfaces at the origin
        # Q_int = Q_int - centroid(Q_int)
        # P_int = P_int - centroid(P_int)

        # find the rotation that minimizes the interface rmsd
        U_int = kabsch(P_int, Q_int)
        P_int = np.dot(P_int, U_int)

        # write_coords("ref.pdb", Q_int)
        # write_coords("model.pdb", P_int)

        # # move the system to the centroid of the interfaces
        Q = Q - centroid(Q)
        P = P - centroid(P)

        # write_coords("ref_1.pdb", Q)
        # write_coords("model_1.pdb", P)

        Q = Q - centroid(Q_int)
        P = P - centroid(P_int)

        # write_coords("ref_2.pdb", Q)
        # write_coords("model_2.pdb", P)

        # apply this rotation to the model
        #  - complexes are now aligned by the interfaces
        P = np.dot(P, U_int)

        # write_coords("ref_i.pdb", Q)
        # write_coords("model_i.pdb", P)

        # Calculate the rmsd of the ligand
        Q_l = Q[l_start: l_end - 1]
        P_l = P[l_start: l_end - 1]

        # write_coords("ref_l.pdb", Q_l)
        # write_coords("model_l.pdb", P_l)

        # this will be the interface-ligand-rmsd
        self.ilrmsd = calc_rmsd(P_l, Q_l)

    def calc_fnat(self, cutoff=5.0):
        """
        Calculate the frequency of native contacts.

        Parameters
        ----------
        cutoff : float
            The cutoff distance for the intermolecular contacts.
        """
        ref_contacts = self.load_contacts(self.reference, cutoff)
        model_contacts = self.load_contacts(self.model, cutoff)
        intersection = ref_contacts & model_contacts
        self.fnat = len(intersection) / float(len(ref_contacts))

    def calc_dockq(self):
        """Calculate the DockQ metric."""
        if self.fnat and self.irmsd and self.lrmsd:
            self.dockq = (
                float(self.fnat)
                + 1 / (1 + (self.irmsd / 1.5) * (self.irmsd / 1.5))
                + 1 / (1 + (self.lrmsd / 8.5) * (self.lrmsd / 8.5))
                ) / 3
        else:
            self.dockq = float("nan")

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
            data["self.model-cluster-ranking"] = None

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
        pdb = read_pdb(pdb_f)

        interface_resdic = {}
        for atom_i, atom_j in get_intermolecular_contacts(pdb, cutoff):

            if atom_i.chain not in interface_resdic:
                interface_resdic[atom_i.chain] = []
            if atom_j.chain not in interface_resdic:
                interface_resdic[atom_j.chain] = []

            if atom_i.resid not in interface_resdic[atom_i.chain]:
                interface_resdic[atom_i.chain].append(atom_i.resid)
            if atom_j.resid not in interface_resdic[atom_j.chain]:
                interface_resdic[atom_j.chain].append(atom_j.resid)

        return interface_resdic

    @staticmethod
    def load_contacts(pdb_f, cutoff=5.0):
        """
        Load residue-based contacts.

        Parameters
        ----------
        pdb_f : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`
            PDB file of the model to have its atoms identified
        cutoff : float, optional
            Cutoff distance for the interface identification.
        """
        con_list = []
        if isinstance(pdb_f, PDBFile):
            pdb_f = pdb_f.rel_path
        structure = read_pdb(pdb_f)
        for atom_i, atom_j in get_intermolecular_contacts(structure, cutoff):
            con = (atom_i.chain, atom_i.resid, atom_j.chain, atom_j.resid)
            con_list.append(con)
        return set(con_list)

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
        data[ident] = {}

        # add this dummy data so we can rank it later
        #  without messing up the order of the output
        if not out_file.exists():
            data[ident]['score'] = 99999.9
            continue

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
    score_rankkey_values = [(i, data[k]['score']) for i, k in enumerate(data)]
    score_rankkey_values.sort(key=lambda x: x[1])

    for i, k in enumerate(score_rankkey_values):
        idx, _ = k
        data[idx + 1]["caprieval_rank"] = i + 1

        if data[idx + 1]['score'] == 99999.9:
            del data[idx + 1]

    # Sort according to the sort key
    rankkey_values = [(i, data[k][sort_key]) for i, k in enumerate(data)]
    rankkey_values.sort(
        key=lambda x: x[1],
        reverse=True if not sort_ascending else False
        )

    _data = {}
    for i, (k, _) in enumerate(rankkey_values):
        _data[i + 1] = data[k + 1]
    data = _data

    if not data:
        # This means there were only "dummy" values
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
    # get the cluster data
    clt_data = dict(((m.clt_rank, m.clt_id), []) for m in model_list)

    # add models to each cluster
    for capri, model in zip(capri_list, model_list):
        clt_data[(model.clt_rank, model.clt_id)].append((capri, model))

    output_dic = {}
    for i, element in enumerate(clt_data):
        data = {}
        number_of_models_in_cluster = len(clt_data[element])
        # TODO: Refactor these ugly try/excepts
        try:
            score_array = [
                e[1].score for e in clt_data[element][:clt_threshold]]
            score_mean, score_stdev = calc_stats(score_array)
        except KeyError:
            score_mean = float("nan")
            score_stdev = float("nan")

        try:
            irmsd_array = [
                e[0].irmsd for e in clt_data[element][:clt_threshold]]
            irmsd_mean, irmsd_stdev = calc_stats(irmsd_array)
        except KeyError:
            irmsd_mean = float("nan")
            irmsd_stdev = float("nan")

        try:
            fnat_array = [e[0].fnat for e in clt_data[element][:clt_threshold]]
            fnat_mean, fnat_stdev = calc_stats(fnat_array)
        except KeyError:
            fnat_mean = float("nan")
            fnat_stdev = float("nan")

        try:
            lrmsd_array = [
                e[0].lrmsd for e in clt_data[element][:clt_threshold]]
            lrmsd_mean, lrmsd_stdev = calc_stats(lrmsd_array)
        except KeyError:
            lrmsd_mean = float("nan")
            lrmsd_stdev = float("nan")

        try:
            dockq_array = [
                e[0].dockq for e in clt_data[element][:clt_threshold]]
            dockq_mean, dockq_stdev = calc_stats(dockq_array)
        except KeyError:
            dockq_mean = float("nan")
            dockq_stdev = float("nan")

        data["cluster_rank"] = element[0]
        data["cluster_id"] = element[1]
        data["n"] = number_of_models_in_cluster
        if number_of_models_in_cluster < clt_threshold:
            # under-evaluated, the mean was divided by a value
            #  larger than the total number of models in the cluster
            data["under_eval"] = "yes"
        else:
            data["under_eval"] = "-"

        data["score"] = score_mean
        data["score_std"] = score_stdev
        data["irmsd"] = irmsd_mean
        data["irmsd_std"] = irmsd_stdev
        data["fnat"] = fnat_mean
        data["fnat_std"] = fnat_stdev
        data["lrmsd"] = lrmsd_mean
        data["lrmsd_std"] = lrmsd_stdev
        data["dockqn"] = dockq_mean
        data["dockq_std"] = dockq_stdev

        output_dic[i] = data

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
