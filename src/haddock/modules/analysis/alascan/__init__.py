"""HADDOCK3 module for alanine scan."""
from pathlib import Path
import os
from haddock.gear.haddockmodel import HaddockModel
from haddock.libs.libcns import prepare_cns_input, prepare_expected_pdb, PDBFile
from haddock.libs.libontology import ModuleIO, TopologyFile
from haddock.libs.libsubprocess import CNSJob
from haddock.gear.config_reader import read_config
from haddock.modules import get_engine
from haddock.modules.base_cns_module import BaseCNSModule
from haddock.modules.analysis.caprieval.capri import CAPRI
from haddock.modules.analysis.alascan.scan import mutate, add_delta_to_bfactor
from haddock.modules.topology.topoaa import generate_topology
from haddock.modules.topology.topoaa import DEFAULT_CONFIG as TOPO_PARAMS
from haddock.modules.topology.topoaa import HaddockModule as Topoaa
from pdbtools import pdb_tidy
import operator
from functools import partial

RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.cfg")


class HaddockModule(BaseCNSModule):
    """HADDOCK3 module for alanine scan."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        cns_script = Path(RECIPE_PATH, "cns", "emref.cns")
        super().__init__(order, path, initial_params, cns_script=cns_script)
        self.top_script = Path(
            RECIPE_PATH,
            '../../topology/topoaa/cns/generate-topology.cns'
            ).read_text()

        self.top_params = read_config(TOPO_PARAMS)

    @classmethod
    def confirm_installation(cls):
        """Confirm if module is installed."""
        return

    def _run(self):
        """Execute module."""

        # Get the models generated in previous step
        try:
            models_to_scan = \
                self.previous_io.retrieve_models(individualize=True)
        except Exception as e:
            self.finish_with_error(e)

        # FIXME: this is definitively wrong
        topoaa = Topoaa(1, path=self.path)
        mol_params = topoaa.params.pop('input')
        mol_params_keys = list(mol_params.keys())[::-1]
        if topoaa.params['limit']:
            mol_params_get = mol_params_keys.pop
        else:
            mol_params_get = partial(operator.getitem, mol_params_keys, -1)
        parameters_for_this_molecule = mol_params[mol_params_get()]

        # generate topology
        input_dic = {}
        topo_jobs = []
        for native in models_to_scan:
            input_dic[native] = []
            interface = CAPRI.identify_interface(native.rel_path)
            self.log(f'Mutating interface of {native.file_name}...')
            for chain in interface:
                for res in interface[chain]:
                    # Mutate and make it look nice
                    mut_pdb_f = mutate(native.rel_path, chain, res, 'ALA')
                    mut_id = \
                        str(mut_pdb_f).split('-')[-1].split('.pdb')[0]
                    # self.debug(f'> Mutating {chain}.{mut_id}')
                    with open(mut_pdb_f, 'r') as fin:
                        lines = list(
                            pdb_tidy.run(fin, strict=False)
                            )  # be explicit in the `strict`

                    with open(mut_pdb_f, 'w') as fout:
                        fout.write(''.join(lines))

                    topology_filename = generate_topology(
                        Path(mut_pdb_f),
                        self.top_script,
                        topoaa.params,
                        parameters_for_this_molecule,
                        default_params_path=topoaa.toppar_path,
                        )
                    # Add new job to the pool
                    output_filename = \
                        Path(mut_pdb_f).with_suffix('.out')

                    job = CNSJob(
                        topology_filename,
                        output_filename,
                        envvars=topoaa.default_envvars(),
                        cns_exec=self.params["cns_exec"],
                        )

                    topo_jobs.append(job)

                    # processed pdb
                    pdb_name = str(mut_pdb_f)
                    mut_pdb = PDBFile(
                        pdb_name.replace('.pdb', '_haddock.pdb')
                        )
                    top_name = pdb_name.replace('.pdb', '_haddock.psf')
                    mut_pdb.topology = TopologyFile(
                        Path(top_name),
                        path=".")
                    # save the mutation id as ori_name
                    mut_pdb.ori_name = mut_id
                    input_dic[native].append(mut_pdb)

        # get the topologies
        self.log(
            f"Generating topologies of mutated PDBs n={len(topo_jobs)}"
            )
        Engine = get_engine(self.params['mode'], self.params)
        engine = Engine(topo_jobs)
        engine.run()
        self.log("Topologies done.")

        alascan_dic = {}
        idx = 1
        mut_jobs = []
        expected = []
        for native in input_dic:
            # prepare this job
            inp_file = prepare_cns_input(
                idx,
                native,
                self.path,
                self.recipe_str,
                self.params,
                "alascan",
                )
            out_file = f"alascan_{idx}.out"
            native_pdb = prepare_expected_pdb(native, idx, ".", "alascan")
            alascan_dic[native_pdb] = []
            expected.append(native_pdb)
            job = CNSJob(inp_file, out_file, envvars=self.envvars)
            mut_jobs.append(job)
            idx += 1
            for mutant in input_dic[native]:
                inp_file = prepare_cns_input(
                    idx,
                    mutant,
                    self.path,
                    self.recipe_str,
                    self.params,
                    "alascan",
                    )
                out_file = f"alascan_{idx}.out"
                expected_mutated_pdb = prepare_expected_pdb(
                    mutant, idx, ".", "alascan"
                    )
                expected_mutated_pdb.ori_name = mutant.ori_name
                expected.append(expected_mutated_pdb)
                alascan_dic[native_pdb].append(expected_mutated_pdb)
                job = CNSJob(inp_file, out_file, envvars=self.envvars)
                mut_jobs.append(job)
                idx += 1

        # Run CNS Jobs
        self.log(f"Calculating energies n={len(mut_jobs)}")
        Engine = get_engine(self.params['mode'], self.params)
        engine = Engine(mut_jobs)
        engine.run()
        self.log("Energy calculations done")

        # Get the weights from the defaults
        _weight_keys = ("w_vdw", "w_elec", "w_desolv", "w_air", "w_bsa")
        weights = {e: self.params[e] for e in _weight_keys}

        for native in alascan_dic:
            if native.is_present():
                native_score = \
                    HaddockModel(native.file_name).calc_haddock_score(
                        **weights
                        )
                native.score = native_score
                for mutated_pdb in alascan_dic[native]:
                    if mutated_pdb.is_present():
                        mutated_pdb.score = \
                            HaddockModel(
                                mutated_pdb.file_name).calc_haddock_score(
                                    **weights
                                    )

        # write output
        output_f = Path('alascan.txt')
        with open(output_f, 'w') as out_fh:
            out_fh.write(f'model\ttype\tscore\tdelta{os.linesep}')
            for native in alascan_dic:
                out_fh.write(
                    f'{native.file_name}\tWT\t{native.score:.2f}\t'
                    f'{float("nan")}{os.linesep}')
                for mutant in alascan_dic[native]:
                    delta_score = mutant.score - native.score
                    mut_id = mutant.ori_name
                    out_fh.write(
                        f'{mutant.file_name}\t{mut_id}\t{mutant.score:.2f}'
                        f'\t{delta_score:.2f}{os.linesep}'
                        )

        # Put the delta in the b-factor
        bfactor_l = []
        for native in alascan_dic:
            for mutant in alascan_dic[native]:
                chain, mut_id = mutant.ori_name.split('_')
                resnum = int(mut_id[1:-1])
                delta = mutant.score - native.score
                ident = f'{chain}.{resnum}'
                # if chain not in bfactor_dic:
                #     bfactor_dic[chain] = {}
                # bfactor_dic[chain][resnum] = delta
                bfactor_l.append((ident, delta))

            _l = [e[1] for e in bfactor_l]
            norm_l = [-1 + 2 * (x - min(_l)) / (max(_l) - min(_l)) for x in _l]
            norm_bfactor_dic = dict(
                (e[0], j) for e, j in zip(bfactor_l, norm_l))

            bfactor_pdb = \
                add_delta_to_bfactor(native.file_name, norm_bfactor_dic)
            self.log(f'Saved visual output to {bfactor_pdb}')

        # Save module information
        io = ModuleIO()
        io.add(expected, "o")
        faulty = io.check_faulty()
        tolerance = self.params["tolerance"]
        if faulty > tolerance:
            _msg = (
                f"{faulty:.2f}% of output was not generated for this module "
                f"and tolerance was set to {tolerance:.2f}%.")
            self.finish_with_error(_msg)
        io.save()
