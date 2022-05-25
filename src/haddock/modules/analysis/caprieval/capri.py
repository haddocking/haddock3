"""CAPRI module."""
import os
import shutil
import tempfile
from math import sqrt
from pathlib import Path

import numpy as np
from fccpy import read_pdb
from fccpy.contacts import get_intermolecular_contacts
from pdbtools import pdb_segxchain

from haddock.libs.libalign import (
    calc_rmsd,
    centroid,
    get_align,
    get_atoms,
    kabsch,
    load_coords,
    make_range,
    )
from haddock.libs.libontology import PDBFile


class CAPRI:
    """CAPRI class."""

    def __init__(
            self,
            reference,
            model_list,
            receptor_chain,
            ligand_chain,
            aln_method,
            path,
            core=None,
            core_model_idx=0,
            **params,
            ):
        self.reference = reference
        self.model_list = model_list
        self.irmsd_dic = {}
        self.lrmsd_dic = {}
        self.ilrmsd_dic = {}
        self.fnat_dic = {}
        self.dockq_dic = {}
        self.atoms = get_atoms(model_list + [reference])
        self.r_chain = receptor_chain
        self.l_chain = ligand_chain
        self.path = path
        # for parallelisation
        self.core = core
        self.core_model_idx = core_model_idx

        # TODO: For scoring we might need to get one alignment per model
        reference = str(reference)
        model = model_list[0].rel_path
        align_func = get_align(aln_method, **params)
        self.model2ref_numbering = align_func(reference, model, path)
        if not self.model2ref_numbering:
            raise CAPRIError("Could not align reference and model")

    def irmsd(self, cutoff=5.0):
        """Calculate the I-RMSD."""
        # Identify reference interface
        ref_interface_resdic = self.identify_interface(self.reference, cutoff)

        # Load interface coordinates
        ref_coord_dic, _ = load_coords(
            self.reference, self.atoms, ref_interface_resdic
            )

        for model in self.model_list:

            mod_coord_dic, _ = load_coords(
                model,
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
            i_rmsd = calc_rmsd(P, Q)
            # write_coords("model_aln.pdb", P)
            # write_coords("ref_aln.pdb", Q)

            self.irmsd_dic[model] = i_rmsd

        return self.irmsd_dic

    def lrmsd(self):
        """Calculate the L-RMSD."""
        ref_coord_dic, _ = load_coords(self.reference, self.atoms)

        for model in self.model_list:

            mod_coord_dic, _ = load_coords(
                model,
                self.atoms,
                numbering_dic=self.model2ref_numbering
                )

            Q = []
            P = []
            # Note: this MUST be sorted since we will use the indexes to
            #  separate between receptor and ligand coordinates
            instersection = sorted(ref_coord_dic.keys() & mod_coord_dic.keys())

            chain_ranges = {}
            for i, segment in enumerate(instersection):
                chain, _, _ = segment
                if chain not in chain_ranges:
                    chain_ranges[chain] = []
                chain_ranges[chain].append(i)

            chain_ranges = make_range(chain_ranges)
            r_start, r_end = chain_ranges[self.r_chain]
            l_start, l_end = chain_ranges[self.l_chain]

            for k in instersection:
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

            # # move to the origin
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
            l_rmsd = calc_rmsd(P_l, Q_l)

            # write_coords("ref.pdb", Q)
            # write_coords("model.pdb", P)

            self.lrmsd_dic[model] = l_rmsd

        return self.lrmsd_dic

    def ilrmsd(self, cutoff=10.0):
        """Calculate the Interface Ligand RMSD."""
        # Identify interface
        ref_interface_resdic = self.identify_interface(self.reference, cutoff)

        # Load interface coordinates
        ref_coord_dic, _ = load_coords(self.reference, self.atoms)

        ref_int_coord_dic, _ = load_coords(
            self.reference, self.atoms, ref_interface_resdic
            )

        for model in self.model_list:

            mod_coord_dic, _ = load_coords(
                model,
                self.atoms,
                numbering_dic=self.model2ref_numbering
                )

            mod_int_coord_dic, _ = load_coords(
                model,
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
            l_start, l_end = chain_ranges[self.l_chain]

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

            # # put interfaces at the origin
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
            i_l_rmsd = calc_rmsd(P_l, Q_l)
            self.ilrmsd_dic[model] = i_l_rmsd

        return self.ilrmsd_dic

    def fnat(self, cutoff=5.0):
        """Calculate the frequency of native contacts."""
        ref_contacts = self.load_contacts(self.reference, cutoff)
        for model in self.model_list:
            model_contacts = self.load_contacts(model, cutoff)
            intersection = ref_contacts & model_contacts
            fnat = len(intersection) / float(len(ref_contacts))
            self.fnat_dic[model] = fnat
        return self.fnat_dic

    def dockq(self):
        """Calculate the DockQ metric."""
        for model in self.model_list:
            irmsd = self.irmsd_dic[model]
            fnat = self.fnat_dic[model]
            lrmsd = self.lrmsd_dic[model]
            dockq = (
                float(fnat)
                + 1 / (1 + (irmsd / 1.5) * (irmsd / 1.5))
                + 1 / (1 + (lrmsd / 8.5) * (lrmsd / 8.5))
                ) / 3
            self.dockq_dic[model] = dockq

        return self.dockq_dic

    def output(
            self,
            clt_threshold,
            sortby_key,
            sort_ascending,
            ):
        """Output the CAPRI results to a .tsv file."""
        self._output_ss(sortby_key, sort_ascending)
        self._output_clt(clt_threshold, sortby_key, sort_ascending,)

    def _output_ss(self, sortby_key, sort_ascending):
        output_l = []
        for model in self.model_list:
            data = {}
            # keep always "model" the first key
            data["model"] = model
            # create the empty rank here so that it will appear
            #  as the second column
            data["caprieval_rank"] = None
            data["score"] = model.score
            if model in self.irmsd_dic:
                data["irmsd"] = self.irmsd_dic[model]
            if model in self.fnat_dic:
                data["fnat"] = self.fnat_dic[model]
            if model in self.lrmsd_dic:
                data["lrmsd"] = self.lrmsd_dic[model]
            if model in self.ilrmsd_dic:
                data["ilrmsd"] = self.ilrmsd_dic[model]
            if model in self.dockq_dic:
                data["dockq"] = self.dockq_dic[model]
            # add cluster data
            data["cluster-id"] = model.clt_id
            data["cluster-ranking"] = model.clt_rank
            data["model-cluster-ranking"] = model.clt_model_rank
            # list of dictionaries
            output_l.append(data)

        if self.core is None:
            capri_name = "capri_ss.tsv"
        else:
            capri_name = "capri_ss_" + str(self.core) + ".tsv"
        output_fname = Path(self.path, capri_name)

        self._dump_file(
            output_l,
            output_fname,
            sortby_key,
            sort_ascending,
            )

    def _output_clt(
            self,
            clt_threshold,
            sortby_key,
            sort_ascending,
            ):
        """Output cluster-based results."""
        has_cluster_info = any(m.clt_id for m in self.model_list)
        if not has_cluster_info:
            return

        # get the cluster data
        clt_data = dict(((m.clt_rank, m.clt_id), []) for m in self.model_list)

        # add models to each cluster
        for model in self.model_list:
            clt_data[(model.clt_rank, model.clt_id)].append(model)

        output_l = []
        for element in clt_data:
            data = {}
            number_of_models_in_cluster = len(clt_data[element])
            # TODO: Refactor these ugly try/excepts
            try:
                score_array = [v.score
                               for v in clt_data[element][: clt_threshold]]
                score_mean, score_stdev = self._calc_stats(
                    score_array, clt_threshold)
            except KeyError:
                score_mean = float("nan")
                score_stdev = float("nan")

            try:
                irmsd_array = [
                    self.irmsd_dic[v]
                    for v in clt_data[element]
                    [: clt_threshold]]
                irmsd_mean, irmsd_stdev = self._calc_stats(
                    irmsd_array, clt_threshold)
            except KeyError:
                irmsd_mean = float("nan")
                irmsd_stdev = float("nan")

            try:
                fnat_array = [self.fnat_dic[v]
                              for v in clt_data[element][:clt_threshold]]
                fnat_mean, fnat_stdev = self._calc_stats(
                    fnat_array, clt_threshold)
            except KeyError:
                fnat_mean = float("nan")
                fnat_stdev = float("nan")

            try:
                lrmsd_array = [self.lrmsd_dic[v]
                               for v in clt_data[element]
                               [: clt_threshold]]
                lrmsd_mean, lrmsd_stdev = self._calc_stats(
                    lrmsd_array, clt_threshold)
            except KeyError:
                lrmsd_mean = float("nan")
                lrmsd_stdev = float("nan")
            try:
                dockq_array = [self.dockq_dic[v]
                               for v in clt_data[element]
                               [: clt_threshold]]
                dockq_mean, dockq_stdev = self._calc_stats(
                    dockq_array, clt_threshold)
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

            output_l.append(data)
        if self.core is None:
            capri_name = "capri_clt.tsv"
        else:
            capri_name = "capri_clt_" + str(self.core) + ".tsv"
        output_fname = Path(self.path, capri_name)

        info_header = "#" * 40 + os.linesep
        info_header += "# `caprieval` cluster-based analysis" + os.linesep
        info_header += "#" + os.linesep
        info_header += f"# > sortby_key={sortby_key}" + os.linesep
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

        self._dump_file(
            output_l,
            output_fname,
            sortby_key,
            sort_ascending,
            info_header=info_header,
            )

    def _dump_file(
            self,
            container,
            output_fname,
            sortby_key,
            sort_ascending,
            info_header="",
            ):

        # rank
        ranked_output_l = self._rank(
            container, key='score',
            ascending=True, core_model_idx=self.core_model_idx
            )

        # sort
        sorted_keys = self._sort(
            ranked_output_l, key=sortby_key, ascending=sort_ascending
            )

        header = "\t".join(ranked_output_l[0].keys())

        if info_header:
            header = info_header + os.linesep + header

        with open(output_fname, "w") as out_fh:
            out_fh.write(header + os.linesep)
            for idx, _ in sorted_keys:
                row_l = []
                for value in ranked_output_l[idx].values():
                    if isinstance(value, Path):
                        row_l.append(str(value))
                    elif isinstance(value, PDBFile):
                        row_l.append(str(value.rel_path))
                    elif isinstance(value, int):
                        row_l.append(f"{value}")
                    elif isinstance(value, str):
                        row_l.append(f"{value}")
                    elif value is None:
                        row_l.append("-")
                    else:
                        row_l.append(f"{value:.3f}")
                out_fh.write("\t".join(row_l) + os.linesep)

    @staticmethod
    def _sort(container, key, ascending):
        # Sort the column
        key_values = [(i, k[key]) for i, k in enumerate(container)]
        key_values.sort(key=lambda x: x[1], reverse=not ascending)
        return key_values

    @staticmethod
    def _rank(container, key, ascending, core_model_idx):
        rankkey_values = [(i, k[key]) for i, k in enumerate(container)]
        rankkey_values.sort(key=lambda x: x[1], reverse=not ascending)
        for i, k in enumerate(rankkey_values, start=1):
            idx, _ = k
            container[idx]["caprieval_rank"] = core_model_idx + i
        return container

    @staticmethod
    def _calc_stats(data, n):
        """Calculate the mean and stdev."""
        mean = sum(data) / n
        var = sum((x - mean) ** 2 for x in data) / n
        stdev = sqrt(var)
        return mean, stdev

    @staticmethod
    def identify_interface(pdb_f, cutoff=5.0):
        """Identify the interface."""
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
        """Load residue-based contacts."""
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
        """Replace the chainID with the segID."""
        temp_f = tempfile.NamedTemporaryFile(delete=False, mode="w+t")
        with open(pdb_path) as fh:
            for line in list(pdb_segxchain.run(fh)):
                temp_f.writelines(line)
        temp_f.close()
        # REPLACE!
        new_pdb_path = shutil.move(temp_f.name, pdb_path)
        return new_pdb_path


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
