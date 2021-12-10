"""Create and manage CNS all-atom topology."""
import shutil
from pathlib import Path

from haddock.libs import libpdb
from haddock.libs.libcns import (
    generate_default_header,
    load_workflow_params,
    prepare_output,
    prepare_single_input,
    )
from haddock.libs.libontology import Format, ModuleIO, PDBFile, TopologyFile
from haddock.libs.libstructure import make_molecules
from haddock.libs.libsubprocess import CNSJob
from haddock.modules import BaseHaddockModule, get_engine


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.cfg")


def generate_topology(input_pdb, step_path, recipe_str, defaults, mol_params,
                      protonation=None):
    """Generate a HADDOCK topology file from input_pdb."""
    # this is a special cases that only applies to topolyaa.
    general_param = load_workflow_params(defaults)

    input_mols_params = load_workflow_params(mol_params, param_header='')

    general_param = general_param + input_mols_params

    link, topology_protonation, \
        trans_vec, tensor, scatter, \
        axis, water_box = generate_default_header(protonation)

    abs_path = input_pdb.resolve().parent.absolute()
    output_pdb_filename = abs_path / (f'{input_pdb.stem}_'
                                      f'haddock{input_pdb.suffix}')
    output_psf_filename = abs_path / (f'{input_pdb.stem}_'
                                      f'haddock.{Format.TOPOLOGY}')
    output = prepare_output(output_psf_filename, output_pdb_filename)

    input_str = prepare_single_input(str(input_pdb.resolve().absolute()))

    inp = general_param + input_str + output + link \
        + topology_protonation + trans_vec + tensor + scatter + axis \
        + water_box + recipe_str

    output_inp_filename = abs_path / f'{input_pdb.stem}.{Format.CNS_INPUT}'
    with open(output_inp_filename, 'w') as output_handler:
        output_handler.write(inp)

    return output_inp_filename


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module to create CNS all-atom topologies."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        cns_script = RECIPE_PATH / "cns" / "generate-topology.cns"
        super().__init__(order, path, initial_params, cns_script)

    @classmethod
    def confirm_installation(cls):
        """Confirm if module is installed."""
        return

    def _run(self):
        """Execute module."""
        molecules = make_molecules(self.params.pop('molecules'))

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
            self.log(f"Molecule {i}: {molecule.file_name.name}")
            models_dic[i] = []
            # Copy the molecule to the step folder
            step_molecule_path = Path(self.path, molecule.file_name.name)
            shutil.copyfile(molecule.file_name, step_molecule_path)

            # Split models
            self.log(
                f"Split models if needed for {step_molecule_path}",
                level='debug',
                )
            # these come already sorted
            splited_models = libpdb.split_ensemble(step_molecule_path)

            # nice variable name, isn't it? :-)
            # molecule parameters are shared among models of the same molecule
            parameters_for_this_molecule = mol_params[mol_params_keys.pop()]

            # Sanitize the different PDB files
            for model in splited_models:
                self.log(f"Sanitizing molecule {model.name}")
                models_dic[i].append(model)

                if self.params['ligand_top_fname']:
                    custom_top = self.params['ligand_top_fname']
                    self.log(f'Using custom topology {custom_top}')
                    libpdb.sanitize(model,
                                    overwrite=True,
                                    custom_topology=custom_top)

                else:
                    libpdb.sanitize(model, overwrite=True)

                # Prepare generation of topologies jobs
                topology_filename = generate_topology(
                    model,
                    self.path,
                    self.recipe_str,
                    self.params,
                    parameters_for_this_molecule,
                    )
                self.log(
                    f"Topology CNS input created in {topology_filename.name}"
                    )

                # Add new job to the pool
                output_filename = Path(
                    model.resolve().parent,
                    f"{model.stem}.{Format.CNS_OUTPUT}",
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

        # Run CNS Jobs
        self.log(f"Running CNS Jobs n={len(jobs)}")
        Engine = get_engine(self.params['mode'], self.params)
        engine = Engine(jobs)
        engine.run()
        self.log("CNS jobs have finished")

        # Check for generated output, fail it not all expected files
        #  are found
        expected = {}
        not_found = []
        for i in models_dic:
            expected[i] = {}
            for j, model in enumerate(models_dic[i]):
                model_name = model.stem
                processed_pdb = Path(
                    self.path,
                    f"{model_name}_haddock.{Format.PDB}"
                    )
                if not processed_pdb.is_file():
                    not_found.append(processed_pdb.name)
                processed_topology = Path(
                    self.path,
                    f"{model_name}_haddock.{Format.TOPOLOGY}"
                    )
                if not processed_topology.is_file():
                    not_found.append(processed_topology.name)

                topology = TopologyFile(processed_topology,
                                        path=self.path)
                pdb = PDBFile(processed_pdb,
                              topology,
                              path=self.path)
                expected[i][j] = pdb

        if not_found:
            self.finish_with_error(
                "Several files were not generated:"
                f" {', '.join(not_found)}"
                )

        # Save module information
        io = ModuleIO()
        for i in expected:
            io.add(expected[i], "o")
        io.save(self.path)
