"""HADDOCK3 scoring module."""
from os import linesep
from pathlib import Path

from haddock.gear.haddockmodel import HaddockModel
from haddock.libs.libcns import prepare_cns_input, prepare_expected_pdb
from haddock.libs.libontology import ModuleIO
from haddock.libs.libsubprocess import CNSJob
from haddock.modules import BaseHaddockModule, get_engine


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.cfg")


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module to perform energy minimization scoring."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        cns_script = Path(RECIPE_PATH, "cns", "scoring.cns")
        super().__init__(order, path, initial_params, cns_script)

    @classmethod
    def confirm_installation(cls):
        """Confirm module is installed."""
        return

    def _run(self):
        """Execute module."""
        # Pool of jobs to be executed by the CNS engine
        jobs = []

        try:
            models_to_score = self.previous_io.retrieve_models(
                individualize=True
                )
        except Exception as e:
            self.finish_with_error(e)

        scored_structure_list = []
        for model_num, model in enumerate(models_to_score, start=1):
            scoring_inp = prepare_cns_input(
                model_num,
                model, self.path,
                self.recipe_str,
                self.params,
                "emscoring",
                native_segid=True)

            scoring_out = Path(self.path, f"emscoring_{model_num}.out")

            # create the expected PDBobject
            expected_pdb = prepare_expected_pdb(
                model, model_num, self.path, "emscoring"
                )
            scored_structure_list.append(expected_pdb)

            job = CNSJob(
                scoring_inp,
                scoring_out,
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

        # Get the weights from the defaults
        _weight_keys = ("w_vdw", "w_elec", "w_desolv", "w_air", "w_bsa")
        weights = {e: self.params[e] for e in _weight_keys}

        # Check for generated output, fail it not all expected files are found
        expected = []
        not_found = []
        for pdb in scored_structure_list:
            if not pdb.is_present():
                not_found.append(pdb.file_name)
            else:
                haddock_score = HaddockModel(pdb.full_name).calc_haddock_score(
                    **weights
                    )

                pdb.score = haddock_score
                expected.append(pdb)

        if not_found:
            # fail if any expected files are found
            self.finish_with_error("Several files were not generated:"
                                   f" {not_found}")

        output_fname = Path(self.path, "emscoring.tsv")
        self.log(f"Saving output to {output_fname.name}")
        with open(output_fname, "w") as fh:
            fh.write(f"structure\tscore{linesep}")
            for pdb in expected:
                # temporary solution to get the original name
                #  here one .pdb = one .psf and the .psf is not changed by the module
                #  so use its name as the original one
                original_name = pdb.topology.file_name.replace('.psf', '.pdb')
                fh.write(f"{pdb.file_name}\t{original_name}\t{pdb.score}{linesep}")
        fh.close()

        # Save module information
        io = ModuleIO()
        io.add(expected, "o")
        io.save(self.path)
