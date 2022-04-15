"""Calculate CAPRI metrics."""
import os
from pathlib import Path

from haddock import log
from haddock.libs.libparallel import Scheduler
from haddock.libs.libutil import parse_ncores
from haddock.modules import BaseHaddockModule
from haddock.modules.analysis.caprieval.capri import CAPRI
from haddock.modules.analysis.caprieval.caprijob import CapriJob


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module to calculate the CAPRI metrics."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, *ignore, init_params=DEFAULT_CONFIG,
                 **everything):
        super().__init__(order, path, init_params)

    @classmethod
    def confirm_installation(cls):
        """Confirm if contact executable is compiled."""
        return

    def _rearrange_output(self, output_name, path, ncores):
        """Combine different capri outputs in a single file."""
        output_fname = Path(path, output_name)
        log.info(f"rearranging output files into {output_fname}")
        keyword = output_name.split(".")[0]
        split_dict = {
            "capri_ss": "model-cluster-ranking",
            "capri_clt": "caprieval_rank"
            }
        if keyword not in split_dict.keys():
            raise Exception(f'Keyword {keyword} does not exist.')
        # Combine files
        with open(output_fname, 'w') as out_file:
            for core in range(ncores):
                tmp_file = Path(path, keyword + "_" + str(core) + ".tsv")
                with open(tmp_file) as infile:
                    if core == 0:
                        content = infile.read().rstrip(os.linesep)
                    else:
                        kw = split_dict[keyword]
                        content = infile.read().split(kw)[1].rstrip(os.linesep)
                    out_file.write(content)
                log.debug(f"File number {core} written")
                tmp_file.unlink()
            # adding linesep to output file
            out_file.write(os.linesep)
        log.info("Completed reconstruction of caprieval files.")
        log.info(f"{output_fname} created.")

    def _run(self):
        """Execute module."""
        # Get the models generated in previous step
        if type(self.previous_io) == iter:
            _e = "This module cannot come after one that produced an iterable."
            self.finish_with_error(_e)

        models = self.previous_io.retrieve_models()

        #  Sort by score
        models.sort()
        best_model_fname = Path(models[0].rel_path)

        if self.params["reference_fname"]:
            reference = Path(self.params["reference_fname"])
        else:
            self.log(
                "No reference was given. "
                "Using the structure with the lowest score from previous step")
            reference = best_model_fname

        # Parallelisation : optimal dispatching of models
        nmodels = len(models)
        ncores = parse_ncores(n=self.params['ncores'], njobs=nmodels)
        base_models = nmodels // ncores
        modulo = nmodels % ncores
        chain_of_idx = [0]
        for core in range(ncores):
            if core < modulo:
                chain_of_idx.append(chain_of_idx[core] + base_models + 1)
            else:
                chain_of_idx.append(chain_of_idx[core] + base_models)
        # starting jobs
        capri_jobs = []
        cluster_info = []
        log.info(f"running Capri Jobs with {ncores} cores")
        for core in range(ncores):
            # init Capri
            capri_obj = CAPRI(
                reference,
                models[chain_of_idx[core]:chain_of_idx[core + 1]],
                receptor_chain=self.params["receptor_chain"],
                ligand_chain=self.params["ligand_chain"],
                aln_method=self.params["alignment_method"],
                path=Path("."),
                lovoalign_exec=self.params["lovoalign_exec"],
                core=core,
                core_model_idx=chain_of_idx[core]
                )
            # Name job
            job_f = Path("capri_ss_" + str(core) + ".tsv")
            # Get cluster info
            has_cluster = any(m.clt_id for m in capri_obj.model_list)
            cluster_info.append(has_cluster)
            # init CapriJob
            job = CapriJob(
                job_f,
                self.params,
                capri_obj
                )
            capri_jobs.append(job)
        # Running parallel Capri Jobs
        capri_engine = Scheduler(capri_jobs, ncores=ncores)
        capri_engine.run()
        # Check correct execution of parallel Capri Jobs
        capri_file_l = []
        not_found = []
        for job in capri_jobs:
            if not job.output.exists():
                jobn = job.output.name
                not_found.append(jobn)
                log.warning(f'Capri results were not calculated for {jobn}')
            else:
                capri_file_l.append(str(job.output))
        if not_found:
            # Not all capri objects were calculated
            self.finish_with_error("Several capri files were not generated:"
                                   f" {not_found}")
        # Post-processing : single structure
        self._rearrange_output(
            "capri_ss.tsv",
            path=capri_obj.path,
            ncores=ncores
            )
        # Post-processing : clusters
        has_cluster_info = any(cluster_info)
        if not has_cluster_info:
            self.log("No cluster information")
        else:
            self._rearrange_output(
                "capri_clt.tsv",
                path=capri_obj.path,
                ncores=ncores
                )
        
        # Sending models to the next step of the workflow
        self.output_models = models
        self.export_output_models()
