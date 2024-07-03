"""Create and manage CNS all-atom topology."""
import operator
from functools import partial
from pathlib import Path

from haddock.core.typing import FilePath, Optional, ParamDict, ParamMap
from haddock.libs import libpdb
from haddock.libs.libcns import (
    generate_default_header,
    load_workflow_params,
    prepare_output,
    prepare_single_input,
    )
from haddock.libs.libontology import Format, PDBFile, TopologyFile
from haddock.libs.libsubprocess import CNSJob
from haddock.modules import get_engine
from haddock.modules.base_cns_module import BaseCNSModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")


def generate_topology(
        input_pdb: Path,
        recipe_str: str,
        defaults: ParamMap,
        mol_params: ParamMap,
        default_params_path: Optional[FilePath] = None,
        ) -> Path:
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

    output_inp_filename = Path(f"{input_pdb.stem}.{Format.CNS_INPUT}")
    output_inp_filename.write_text(inp)

    return output_inp_filename


class HaddockModule(BaseCNSModule):
    """HADDOCK3 module to create CNS all-atom topologies."""

    name = RECIPE_PATH.name

    def __init__(
            self,
            order: int,
            path: Path,
            initial_params: FilePath = DEFAULT_CONFIG,
            ) -> None:
        cns_script = RECIPE_PATH / "cns" / "generate-topology.cns"
        super().__init__(order, path, initial_params, cns_script=cns_script)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm if module is installed."""
        return

    def _run(self) -> None:
        """Execute module."""
        md5_dic: dict[int, dict[int, str]] = {}
        self.params.pop("molecules")
        molecules: list[list[Path]] = []
        if self.order == 0:
            _molecules = self.previous_io.output
            for i, models in enumerate(_molecules, start=1):
                molecules.append([model.rel_path for model in models.values()])
                md5_dic[i] = {
                    j: model.md5
                    for j, model in enumerate(models.values())
                    }
        else:
            # in case topoaa is not the first step, the topology is rebuilt for
            # each retrieved model
            _molecules = self.previous_io.retrieve_models()
            molecules = [[mol.rel_path] for mol in _molecules]

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

        for i, models in enumerate(molecules, start=1):
            self.log(f"Molecule {i}")
            models_dic[i] = []
            # nice variable name, isn't it? :-)
            # molecule parameters are shared among models of the same molecule
            parameters_for_this_molecule = mol_params[mol_params_get()]

            for model in models:
                self.log(f"Sanitizing model {model.name}")
                custom_top: Optional[FilePath] = None
                if self.params["ligand_top_fname"]:
                    custom_top = self.params["ligand_top_fname"]
                    self.log(f"Using custom topology {custom_top}")
                libpdb.sanitize(
                    model,
                    overwrite=True,
                    custom_topology=custom_top,
                    )
                models_dic[i].append(model)

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
        Engine = get_engine(self.params["mode"], self.params)
        engine = Engine(jobs)
        engine.run()
        self.log("CNS jobs have finished")

        # Check for generated output, fail if not all expected files
        #  are found
        expected: dict[int, dict[int, PDBFile]] = {}
        for i, models in models_dic.items():
            expected.setdefault(i, {})
            for j, model in enumerate(models):
                model_name = model.stem
                processed_pdb = Path(f"{model_name}_haddock.{Format.PDB}")
                processed_topology = Path(
                    f"{model_name}_haddock.{Format.TOPOLOGY}"
                    )
                topology = TopologyFile(processed_topology, path=".")
                try:
                    md5 = md5_dic[i][j]
                except KeyError:
                    md5 = None
                pdb = PDBFile(
                    file_name=processed_pdb,
                    topology=topology,
                    path=".",
                    md5=md5,
                    )
                pdb.ori_name = model.stem
                expected[i][j] = pdb

        # Save module information
        self.output_models = list(expected.values())  # type: ignore
        self.export_io_models(faulty_tolerance=self.params["tolerance"])
