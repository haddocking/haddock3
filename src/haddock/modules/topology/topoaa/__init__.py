"""Create and manage CNS all-atom topology."""
import operator
import os
import re
from functools import partial
from pathlib import Path

from haddock.libs import libpdb
from haddock.libs.libcns import (
    generate_default_header,
    load_workflow_params,
    prepare_output,
    prepare_single_input,
    )
from haddock.libs.libontology import Format, PDBFile, TopologyFile
from haddock.libs.libstructure import make_molecules
from haddock.libs.libsubprocess import CNSJob
from haddock.modules import get_engine
from haddock.modules.base_cns_module import BaseCNSModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")


def generate_topology(
        input_pdb,
        recipe_str,
        defaults,
        mol_params,
        default_params_path=None,
        ):
    """Generate a HADDOCK topology file from input_pdb."""
    # generate params headers
    general_param = load_workflow_params(**defaults)
    input_mols_params = load_workflow_params(param_header='', **mol_params)
    general_param = general_param + input_mols_params

    # generate default headers
    link, trans_vec, tensor, scatter, axis, water_box = \
        generate_default_header(path=default_params_path)

    output = prepare_output(
        output_pdb_filename=f'{input_pdb.stem}_haddock{input_pdb.suffix}',
        output_psf_filename=f'{input_pdb.stem}_haddock.{Format.TOPOLOGY}',
        )

    input_str = prepare_single_input(str(input_pdb))

    inp_parts = (
        general_param,
        input_str,
        output,
        link,
        trans_vec,
        tensor,
        scatter,
        axis,
        water_box,
        recipe_str,
        )

    inp = "".join(inp_parts)

    output_inp_filename = Path(f'{input_pdb.stem}.{Format.CNS_INPUT}')
    output_inp_filename.write_text(inp)

    return output_inp_filename


class HaddockModule(BaseCNSModule):
    """HADDOCK3 module to create CNS all-atom topologies."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        cns_script = RECIPE_PATH / "cns" / "generate-topology.cns"
        super().__init__(order, path, initial_params, cns_script=cns_script)

    @classmethod
    def confirm_installation(cls):
        """Confirm if module is installed."""
        return

    @staticmethod
    def get_md5(ensemble_f):
        """Get MD5 hash of a multi-model PDB file."""
        md5_dic = {}
        text = Path(ensemble_f).read_text()
        lines = text.split(os.linesep)
        REMARK_lines = (line for line in lines if line.startswith('REMARK'))
        remd5 = re.compile(r"^[a-f0-9]{32}$")
        for line in REMARK_lines:
            parts = line.strip().split()

            try:
                idx = parts.index('MODEL')
            except ValueError:  # MODEL not in parts, this line can be ignored
                continue

            # check if there's a md5 hash in line
            for part in parts:
                group = remd5.fullmatch(part)
                if group:
                    # the model num comes after the MODEL
                    model_num = int(parts[idx + 1])
                    md5_dic[model_num] = group.string  # md5 hash
                    break

        return md5_dic

    def _run(self):
        """Execute module."""
        if self.order == 0:
            # topoaa is the first step in the workflow
            molecules = make_molecules(self.params.pop('molecules'))

        else:
            # in case topoaa is not the first step it reads the input
            # molecules from the previous step but it only takes the
            # first model because all models will have the same
            # identity, as is the case for all PDBs created by rigidbody.
            # In these cases, it is likely that each PDB contains the
            # complex structure instead of the separated chains.
            # Currently, we have not implemented a way to split chains
            # and recreate the topology of the chains separately.
            _molecules = self.previous_io.retrieve_models()[0]
            molecules = make_molecules([_molecules.rel_path], no_parent=True)

        # extracts `input` key from params. The `input` keyword needs to
        # be treated separately
        mol_params = {}
        for k in list(self.params.keys()):
            if k.startswith("mol") and k[3:].isdigit():
                mol_params[k] = self.params.pop(k)

        # to facilite the for loop down the line, we create a list with the keys
        # of `mol_params` with inverted order (we will use .pop)
        mol_params_keys = list(mol_params.keys())[::-1]

        if self.params['limit']:
            mol_params_get = mol_params_keys.pop
        else:
            mol_params_get = partial(operator.getitem, mol_params_keys, -1)

        # Pool of jobs to be executed by the CNS engine
        jobs = []

        models_dic = {}
        ens_dic = {}
        for i, molecule in enumerate(molecules, start=1):
            self.log(f"Molecule {i}: {molecule.file_name.name}")
            models_dic[i] = []
            # Copy the molecule to the step folder

            # Split models
            self.log(
                f"Split models if needed for {molecule.with_parent}",
                level='debug',
                )
            # these come already sorted
            splited_models = libpdb.split_ensemble(
                molecule.with_parent,
                dest=Path.cwd(),
                )

            # get the MD5 hash of each model
            ens_dic[i] = self.get_md5(molecule.with_parent)

            # nice variable name, isn't it? :-)
            # molecule parameters are shared among models of the same molecule
            parameters_for_this_molecule = mol_params[mol_params_get()]

            # Sanitize the different PDB files
            relative_paths_models = (
                Path(_p.name) if _p.is_absolute() else _p
                for _p in splited_models
                )

            for model in relative_paths_models:
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
                    self.recipe_str,
                    self.params,
                    parameters_for_this_molecule,
                    default_params_path=self.toppar_path,
                    )

                self.log(
                    f"Topology CNS input created in {topology_filename.name}"
                    )

                # Add new job to the pool
                output_filename = Path(f"{model.stem}.{Format.CNS_OUTPUT}")

                job = CNSJob(
                    topology_filename,
                    output_filename,
                    envvars=self.envvars,
                    cns_exec=self.params["cns_exec"],
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
        for i in models_dic:
            expected[i] = {}
            md5_dic = ens_dic[i]
            for j, model in enumerate(models_dic[i]):
                md5_hash = None
                try:
                    model_id = int(model.stem.split('_')[-1])
                except ValueError:
                    model_id = ''

                if model_id in md5_dic:
                    md5_hash = md5_dic[model_id]

                model_name = model.stem
                processed_pdb = Path(f"{model_name}_haddock.{Format.PDB}")
                processed_topology = \
                    Path(f"{model_name}_haddock.{Format.TOPOLOGY}")

                topology = TopologyFile(processed_topology, path=".")
                pdb = PDBFile(
                    file_name=processed_pdb,
                    topology=topology,
                    path=".",
                    md5=md5_hash)
                pdb.ori_name = model.stem
                expected[i][j] = pdb

        # Save module information
        self.output_models = list(expected.values())
        self.export_output_models(faulty_tolerance=self.params["tolerance"])
