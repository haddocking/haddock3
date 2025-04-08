"""Create and manage CNS all-atom topology.

The ``[topoaa]`` module is dedicated to the generation of CNS compatible
parameters (.param) and topologies (.psf) for each of the input structures.

It will:
- Detect missing atoms, including hydrogens
- Re-build them when missing
- Build and write out topologies (psf) and coordinates (pdb) files

This module is a pre-requisite to run any downstream modules using CNS.
Having access to parameters and topology is mandatory for any kind
of EM/MD related tasks.
Therefore this is the reason why the module [topoaa] is often used as first
module in a workflow.

Note that for non-standard bio-molecules
(apart from standard amino-acids, some modified ones, DNA, RNA, ions
and carbohydrates ... see `detailed list of supported molecules 
<https://wenmr.science.uu.nl/haddock2.4/library>`_),
such as small-molecules, parameters and topology must be obtained and provided
by the user, as there is currently no built-in solution to generate
them on the fly.
"""

import operator
import os
import re
from functools import partial
from pathlib import Path

from haddock.core.defaults import MODULE_DEFAULT_YAML, cns_exec
from haddock.core.typing import FilePath, Optional, ParamDict, ParamMap, Union
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
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


def generate_topology(
    input_pdb: Path,
    recipe_str: str,
    defaults: ParamMap,
    mol_params: ParamMap,
    default_params_path: Optional[FilePath] = None,
    write_to_disk: Optional[bool] = True,
) -> Union[Path, str]:
    """Generate a HADDOCK topology file from input_pdb."""
    # generate params headers
    general_param = load_workflow_params(**defaults)
    input_mols_params = load_workflow_params(param_header="", **mol_params)
    general_param = general_param + input_mols_params

    # generate default headers
    link, trans_vec, tensor, scatter, axis, water_box = generate_default_header(
        path=default_params_path
    )

    output = prepare_output(
        output_pdb_filename=f"{input_pdb.stem}_haddock{input_pdb.suffix}",
        output_psf_filename=f"{input_pdb.stem}_haddock.{Format.TOPOLOGY}",
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

    if write_to_disk:
        output_inp_filename = Path(f"{input_pdb.stem}.{Format.CNS_INPUT}")
        output_inp_filename.write_text(inp)
        return output_inp_filename
    else:
        return inp


class HaddockModule(BaseCNSModule):
    """HADDOCK3 module to create CNS all-atom topologies."""

    name = RECIPE_PATH.name

    def __init__(
        self, order: int, path: Path, initial_params: FilePath = DEFAULT_CONFIG
    ) -> None:
        cns_script = RECIPE_PATH / "cns" / "generate-topology.cns"
        super().__init__(order, path, initial_params, cns_script=cns_script)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm if module is installed."""
        return

    @staticmethod
    def get_md5(ensemble_f: FilePath) -> dict[int, str]:
        """Get MD5 hash of a multi-model PDB file."""
        md5_dic: dict[int, str] = {}
        text = Path(ensemble_f).read_text()
        lines = text.split(os.linesep)
        REMARK_lines = (line for line in lines if line.startswith("REMARK"))
        remd5 = re.compile(r"^[a-f0-9]{32}$")
        for line in REMARK_lines:
            parts = line.strip().split()

            try:
                idx = parts.index("MODEL")
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

    @staticmethod
    def get_ensemble_origin(ensemble_f: FilePath) -> dict[int, str]:
        """Try to find origin for each model in ensemble.

        Parameters
        ----------
        ensemble_f : FilePath
            Path to a pdb file containing an ensemble.

        Returns
        -------
        origin_dic : dict[int, str]
            Dictionary holding as keys the modelID and values its origin.
        """
        origin_dic: dict[int, str] = {}
        text = Path(ensemble_f).read_text()
        lines = text.split(os.linesep)
        REMARK_lines = (line for line in lines if line.startswith("REMARK"))
        re_origin = re.compile(
            r"REMARK\s+MODEL\s+(\d+)\s+(FROM|from|From)\s+(([\w_-]+\.?)+)"
        )
        for line in REMARK_lines:
            if match := re_origin.search(line):
                model_num = int(match.group(1).strip())
                original_path = match.group(3).strip()
                original_name = Path(original_path).stem
                origin_dic[model_num] = original_name
        return origin_dic

    def _run(self) -> None:
        """Execute module."""
        if self.order == 0:
            # topoaa is the first step in the workflow
            molecules = make_molecules(self.params.pop("molecules"))

        else:
            # in case topoaa is not the first step, the topology is rebuilt for
            # each retrieved model
            _molecules = self.previous_io.retrieve_models()
            molecules_paths: list[Path] = [mol.rel_path for mol in _molecules]  # type: ignore
            molecules = make_molecules(molecules_paths, no_parent=True)

        # extracts `input` key from params. The `input` keyword needs to
        # be treated separately
        mol_params: ParamDict = {}
        for k in list(self.params.keys()):
            if k.startswith("mol") and k[3:].isdigit():
                mol_params[k] = self.params.pop(k)

        # to facilitate the for loop down the line, we create a list with the
        #  keys of `mol_params` with inverted order (we will use .pop)
        mol_params_keys = list(mol_params.keys())[::-1]

        # limit is only useful when order == 0
        if self.order == 0 and self.params["limit"]:
            mol_params_get = mol_params_keys.pop

        # `else` is used in any case where limit is False.
        else:
            mol_params_get = partial(operator.getitem, mol_params_keys, -1)

        # Pool of jobs to be executed by the CNS engine
        jobs: list[CNSJob] = []

        models_dic: dict[int, list[Path]] = {}
        ens_dic: dict[int, dict[int, str]] = {}
        origi_ens_dic: dict[int, dict[int, str]] = {}
        for i, molecule in enumerate(molecules, start=1):
            self.log(f"Molecule {i}: {molecule.file_name.name}")
            models_dic[i] = []
            # Copy the molecule to the step folder

            # Split models
            self.log(
                f"Split models if needed for {molecule.with_parent}",
                level="debug",
            )
            # these come already sorted
            splited_models = libpdb.split_ensemble(
                molecule.with_parent,
                dest=Path.cwd(),
            )

            # get the MD5 hash of each model
            ens_dic[i] = self.get_md5(molecule.with_parent)
            origi_ens_dic[i] = self.get_ensemble_origin(molecule.with_parent)
            # nice variable name, isn't it? :-)
            # molecule parameters are shared among models of the same molecule
            parameters_for_this_molecule = mol_params[mol_params_get()]

            for task_id, model in enumerate(splited_models):
                self.log(f"Sanitizing molecule {model.name}")
                models_dic[i].append(model)

                if self.params["ligand_top_fname"]:
                    custom_top = self.params["ligand_top_fname"]
                    self.log(f"Using custom topology {custom_top}")
                    libpdb.sanitize(
                        model,
                        overwrite=True,
                        custom_topology=custom_top,
                    )

                else:
                    libpdb.sanitize(model, overwrite=True)

                # Prepare generation of topologies jobs
                topoaa_input = generate_topology(
                    model,
                    self.recipe_str,
                    self.params,
                    parameters_for_this_molecule,
                    default_params_path=self.toppar_path,
                    write_to_disk=self.params["debug"],
                )

                self.log("Topology CNS input created")

                # Add new job to the pool
                output_filename = Path(f"{model.stem}.{Format.CNS_OUTPUT}")
                err_fname = f"{model.stem}.cnserr"
                job = CNSJob(
                    topoaa_input,
                    output_filename,
                    err_fname,
                    envvars=self.envvars,
                    cns_exec=cns_exec,
                )

                jobs.append(job)

        # Run CNS Jobs
        self.log(f"Running CNS Jobs n={len(jobs)}")
        Engine = get_engine(self.params["mode"], self.params)
        engine = Engine(jobs)
        engine.run()
        self.log("CNS jobs have finished")

        # Check for generated output, fail it not all expected files
        #  are found
        expected: dict[int, dict[int, PDBFile]] = {}
        for i in models_dic:
            expected[i] = {}
            md5_dic = ens_dic[i]
            origin_names = origi_ens_dic[i]
            for j, model in enumerate(models_dic[i]):
                md5_hash = None
                try:
                    model_id = int(model.stem.split("_")[-1])
                except ValueError:
                    model_id = 0

                if model_id in md5_dic:
                    md5_hash = md5_dic[model_id]

                model_name = model.stem
                processed_pdb = Path(f"{model_name}_haddock.{Format.PDB}")
                processed_topology = Path(f"{model_name}_haddock.{Format.TOPOLOGY}")

                # Check if origin or md5 is available
                if md5_hash or model_id in origin_names.keys():
                    # Select prefix
                    if md5_hash:
                        prefix_name = md5_hash
                    else:
                        prefix_name = origin_names[model_id]
                    # Check if topology and file created
                    if processed_pdb.exists() and processed_topology.exists():
                        # Build new filename
                        model_name = f"{prefix_name}_from_{model_name}"
                        # Rename files
                        processed_pdb = processed_pdb.rename(
                            f"{model_name}_haddock.{Format.PDB}"
                        )
                        processed_topology = processed_topology.rename(
                            f"{model_name}_haddock.{Format.TOPOLOGY}"
                        )

                topology = TopologyFile(processed_topology, path=".")
                pdb = PDBFile(
                    file_name=processed_pdb,
                    topology=topology,
                    path=".",
                    md5=md5_hash,
                )
                pdb.ori_name = model.stem
                expected[i][j] = pdb

        # Save module information
        self.output_models = list(expected.values())  # type: ignore
        self.export_io_models(faulty_tolerance=self.params["tolerance"])
