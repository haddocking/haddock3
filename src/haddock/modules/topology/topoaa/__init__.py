"""Create and manage CNS all-atom topology."""
import shutil
from pathlib import Path

from haddock import log
from haddock.libs import libpdb
from haddock.libs.libcns import (
    generate_default_header,
    load_workflow_params,
    prepare_output,
    prepare_single_input,
    )
from haddock.libs.libontology import Format, ModuleIO, PDBFile, TopologyFile
from haddock.libs.libparallel import Scheduler
from haddock.libs.libstructure import make_molecules_from_pathlist
from haddock.libs.libsubprocess import CNSJob
from haddock.modules import BaseHaddockModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.cfg")


def generate_topology(input_pdb, step_path, recipe_str, defaults, mol_params,
                      protonation=None):
    """Generate a HADDOCK topology file from input_pdb."""
    # this is a special cases that only applies to topolyaa.
    general_param = load_workflow_params(defaults)

    input_mols_params = load_workflow_params(mol_params, param_header='')

    general_param = general_param + input_mols_params

    param, top, link, topology_protonation, \
        trans_vec, tensor, scatter, \
        axis, water_box = generate_default_header(protonation)

    abs_path = input_pdb.resolve().parent.absolute()
    output_pdb_filename = abs_path / (f'{input_pdb.stem}_'
                                      f'haddock{input_pdb.suffix}')
    output_psf_filename = abs_path / (f'{input_pdb.stem}_'
                                      f'haddock.{Format.TOPOLOGY}')
    output = prepare_output(output_psf_filename, output_pdb_filename)

    input_str = prepare_single_input(str(input_pdb.resolve().absolute()))

    inp = general_param + param + top + input_str + output + link \
        + topology_protonation + trans_vec + tensor + scatter + axis \
        + water_box + recipe_str

    output_inp_filename = abs_path / f'{input_pdb.stem}.{Format.CNS_INPUT}'
    with open(output_inp_filename, 'w') as output_handler:
        output_handler.write(inp)

    return output_inp_filename


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module to create CNS all-atom topologies."""

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        cns_script = RECIPE_PATH / "cns" / "generate-topology.cns"
        super().__init__(order, path, initial_params, cns_script)

    @classmethod
    def confirm_installation(cls):
        """Confirm if module is installed."""
        return

    def run(self, molecules, **params):
        """Execute module."""
        log.info("Running [allatom] module")
        log.info("Generating topologies")

        super().run(params)

        molecules = make_molecules_from_pathlist(molecules)
        # extracts `input` key from params. The `input` keyword needs to
        # be treated separately
        mol_params = self.params.pop('input')
        # to facilite the for loop down the line, we create a list with the keys
        # of `mol_params` with inverted order (we will use .pop)
        mol_params_keys = list(mol_params.keys())[::-1]

        # Pool of jobs to be executed by the CNS engine
        jobs = []

        models_dic = {}
        for i, molecule in enumerate(molecules, start=1):
            log.info(f"{i} - {molecule.file_name}")
            models_dic[i] = {}
            # Copy the molecule to the step folder
            step_molecule_path = Path(self.path, molecule.file_name.name)
            shutil.copyfile(molecule.file_name, step_molecule_path)

            # Split models
            log.info(f"Split models if needed for {step_molecule_path}")
            splited_models = libpdb.split_ensemble(step_molecule_path)

            # molecule parameters are shared among models of the same molecule
            parameters_for_this_molecule = mol_params[mol_params_keys.pop()]

            # Sanitize the different PDB files
            for j, model in enumerate(splited_models):
                models_dic[i][j] = []
                chains = libpdb.split_by_chain(model)

                # Each chain must become a model
                for chained_model in chains:

                    log.info(f"Sanitizing molecule {chained_model.name}")
                    models_dic[i][j].append(chained_model)

                    if self.params['ligand_top_fname']:
                        custom_top = self.params['ligand_top_fname']
                        log.info(f'Using custom topology {custom_top}')
                        libpdb.sanitize(chained_model,
                                        overwrite=True,
                                        custom_topology=custom_top)

                    else:
                        libpdb.sanitize(chained_model, overwrite=True)

                    # Prepare generation of topologies jobs
                    topology_filename = generate_topology(
                        chained_model,
                        self.path,
                        self.recipe_str,
                        self.params,
                        parameters_for_this_molecule,
                        )
                    log.info("Topology CNS input created"
                             f" in {topology_filename}")

                    # Add new job to the pool
                    output_filename = Path(
                        chained_model.resolve().parent,
                        f"{chained_model.stem}.{Format.CNS_OUTPUT}",
                        )

                    job = CNSJob(
                        topology_filename,
                        output_filename,
                        cns_folder=self.cns_folder_path,
                        modpath=self.path,
                        config_path=self.params['config_path'],
                        cns_exec=self.params['cns_exec'],
                        )

                    jobs.append(job)

        # Run CNS engine
        log.info(f"Running CNS engine with {len(jobs)} jobs")
        engine = Scheduler(jobs, ncores=self.params['ncores'])
        engine.run()
        log.info("CNS engine has finished")

        # Check for generated output, fail it not all expected files
        #  are found
        expected = {}
        not_found = []
        to_be_cleaned = []
        for i, mol_id in enumerate(models_dic):
            expected[i] = {}
            for j, chained_model_id in enumerate(models_dic[mol_id]):

                chain_model_list = models_dic[mol_id][chained_model_id]
                topology_list = []
                processed_chained_model_list = []
                for sub_model in chain_model_list:

                    processed_pdb = Path(
                        self.path,
                        f"{sub_model.stem}_haddock.{Format.PDB}"
                        )

                    processed_topology = Path(
                        self.path,
                        f"{sub_model.stem}_haddock.{Format.TOPOLOGY}"
                        )

                    if not processed_pdb.exists():
                        not_found.append(processed_pdb.name)

                    if not processed_topology.exists():
                        not_found.append(processed_topology.name)

                    topology_list.append(processed_topology)
                    processed_chained_model_list.append(processed_pdb)

                    to_be_cleaned.append(sub_model)

                m_stem = '_'.join(chain_model_list[0].stem.split('_')[:-1])
                model_fname = Path(self.path, f"{m_stem}_haddock.{Format.PDB}")

                libpdb.merge(processed_chained_model_list, model_fname)

                to_be_cleaned.extend(processed_chained_model_list)

                topology = []
                for top in topology_list:
                    topology.append(TopologyFile(top, path=self.path))

                pdb = PDBFile(model_fname, topology, path=self.path)

                expected[i][j] = pdb

        if not_found:
            self.finish_with_error("Several files were not generated:"
                                   f" {not_found}")

        for temp_pdb in to_be_cleaned:
            temp_pdb.unlink()

        # Save module information
        io = ModuleIO()
        for i in expected:
            io.add(expected[i], "o")
        io.save(self.path)
