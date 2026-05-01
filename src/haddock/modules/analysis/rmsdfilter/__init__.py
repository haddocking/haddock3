"""RMDS-based filtering module.

This module calculates RMSD for input models against a user-supplied reference structure(s)
and filters out models based on RMSD threshold. The idea behind this module is to simplify 
removal of a priori incorrect models, for example those generated with diffusion algorithms, 
such as (partially) unfolded antibodies. 

The following file is generated:
- **rmsdfilter_ss.tsv**: a table with the global RMSD and score (if available) for each model.
- **rmsdfilter_ss_multiref.tsv**: when multiple references are provided, a table with RMSD for every (model, reference) pair.

This module will terminate with an error message in the following cases:
* ``reference_fname`` is not set (opposite to `caprieval`).
* Alignment fails for all models, i.e. all RMSD values are NaN.
* All models exceed the RMSD threshold, i.e. no models pass through filtering.
* Models have already been clustered, i.e. models carry cluster metadata.

For more details about this module, please `refer to the haddock3 user manual
<https://www.bonvinlab.org/haddock3-user-manual/modules/analysis.html#rmsdfilter-module>`_
"""

from math import isnan
from pathlib import Path

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import FilePath, Union
from haddock.libs.libalign import get_align
from haddock.libs.libaa2cg import martinize
from haddock.libs.libontology import PDBFile
from haddock.libs.libparallel import Scheduler
from haddock.libs.libstructure import find_ff
from haddock.libs.libpdb import handle_input_reference
from haddock.modules import BaseHaddockModule
from haddock.modules.analysis.rmsdfilter.rmsdfilter import RMSDFilter


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


def _sort_key(row, sortby: str) -> float:
    model, rmsd = row
    if sortby == "score":
        val = model.score
        if val is None or (isinstance(val, float) and isnan(val)):
            return float("inf")
        return val
    return float("inf") if isnan(rmsd) else rmsd


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module to filter models by RMSD."""

    name = RECIPE_PATH.name

    def __init__(
        self,
        order: int,
        path: Path,
        init_params: FilePath = DEFAULT_CONFIG,
    ) -> None:
        super().__init__(order, path, init_params)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm if module is installed."""
        return

    @staticmethod
    def is_nested(models: list[Union[PDBFile, list[PDBFile]]]) -> bool:
        for model in models:
            if isinstance(model, list):
                return True
        return False

    def _run(self) -> None:
        """Execute module."""
        # Get the models generated in previous step
        if type(self.previous_io) == iter:
            self.finish_with_error(
                "[rmsdfilter] module cannot come after one that "
                "produced an iterable."
            )

        # Get the models generated in previous step
        models = self.previous_io.retrieve_models(individualize=True)
        if self.is_nested(models):
            raise ValueError(
                "[rmsdfilter] module cannot be executed after "
                "modules that produce a nested list of models."
            )
        
        # Check if cluster info is present
        if any(m.clt_id is not None for m in models):
            self.finish_with_error(
                "Models have been clustered!"
                "[rmsdfilter] cannot be performed after clustering - "
                "filtering individual models after clustering would leave "
                "remaining models with stale and inconsistent cluster assignments."
            )

        # load reference(s)
        if self.params["reference_fname"]:
            reference = Path(self.params["reference_fname"])
            references = handle_input_reference(reference)
        else:
            self.finish_with_error(
                "[rmsdfilter] No valid reference structure(s) provided!"
                "A reference structure is required for this module to work,"
                "please set the 'reference_fname' parameter to a valid PDB file."
            )

        # Detect force field (aa of cg)
        ff = find_ff(models)
        # if cg, convert reference to cg
        if ff == "martini2":
            references = [
                Path(martinize(ref, self.path.resolve().parent, False))
                for ref in references
            ]

        # Build alignment function (to be passed to each job)
        align_func = get_align(
            method=self.params["alignment_method"],
            lovoalign_exec=self.params["lovoalign_exec"],
            keep_hetatm=self.params["keep_hetatm"],
        )

        # Create one job per (model, reference) pair
        jobs: list[RMSDFilter] = []
        for i, model in enumerate(models, start=1):
            for ref_id, reference in enumerate(references, start=1):
                jobs.append(
                    RMSDFilter(
                        identificator=i,
                        model=model,
                        reference=reference,
                        path=Path("."),
                        params=self.params,
                        align_func=align_func,
                        ref_id=ref_id,
                    )
                )

        engine = Scheduler(
            tasks=jobs,
            ncores=self.params["ncores"],
            max_cpus=self.params["max_cpus"],
        )
        engine.run()
        results: list[RMSDFilter] = engine.results

        # For each model take the minimum RMSD across all references
        rmsd_map: dict[int, float] = {}
        for job in results:
            current = rmsd_map.get(job.identificator, float("nan"))
            job_rmsd = job.rmsd
            if isnan(current):
                rmsd_map[job.identificator] = job_rmsd
            elif not isnan(job_rmsd):
                rmsd_map[job.identificator] = min(current, job_rmsd)

        # Map identificator back to model
        id_to_model: dict[int, PDBFile] = {
            i: m for i, m in enumerate(models, start=1)
        }

        # Sort models to write in output tsv
        sortby = self.params["sortby"]
        sort_ascending = self.params["sort_ascending"]

        # Check if score is available
        has_score = any(
            m.score is not None and not (
                isinstance(m.score, float) and isnan(m.score)
            )
            for m in models
        )
        if sortby == "score" and not has_score:
            self.log(
                "Cannot sort models by score, falling back to sorting by RMSD."
            )
            sortby = "rmsd"

        rows = [
            (id_to_model[i], rmsd_map[i])
            for i in sorted(rmsd_map.keys())
        ]

        rows.sort(key=lambda row: _sort_key(row, sortby), reverse=not sort_ascending)

        # Compute filtering stats for the tsv header
        valid_rows = [(m, r) for m, r in rows if not isnan(r)]
        nan_count = len(rows) - len(valid_rows)
        filtered = [m for m, r in valid_rows if r <= self.params["threshold"]]
        percent_filtered = (1 - len(filtered) / len(models)) * 100

        with open("rmsdfilter_ss.tsv", "w") as fh:
            fh.write(
                f"# RMSD filtering threshold is set to {self.params['threshold']:.3f} Å; "
                f"{len(filtered)} model(s) were kept; {percent_filtered:.2f}% were filtered out. "
                "This file contains all models for user information.\n"
            )
            fh.write("model\tscore\trmsd\n")
            for model, rmsd in rows:
                score_str = (
                    f"{model.score:.3f}"
                    if model.score is not None and not (
                        isinstance(model.score, float) and isnan(model.score)
                    )
                    else "nan"
                )
                rmsd_str = f"{rmsd:.3f}" if not isnan(rmsd) else "nan"
                fh.write(f"{model.rel_path}\t{score_str}\t{rmsd_str}\n")

        # When multiple references were used, write multiref tsv,
        # so it is possible to trace from which reference RMSD came from
        if len(references) > 1:
            multiref_rows = sorted(
                results,
                key=lambda job: (str(id_to_model[job.identificator].rel_path), job.ref_id),
            )
            with open("rmsdfilter_ss_multiref.tsv", "w") as fh:
                fh.write("model\tref_id\trmsd\n")
                for job in multiref_rows:
                    model = id_to_model[job.identificator]
                    rmsd_str = f"{job.rmsd:.3f}" if not isnan(job.rmsd) else "nan"
                    fh.write(
                        f"{model.rel_path}\t{job.ref_id}\t{rmsd_str}\n"
                    )

        if nan_count > 0:
            self.log(
                f"{100 * nan_count / len(rows):6.2f}% of models had NaN RMSD "
                "(alignment failed) and will be excluded from filtering."
            )

        if not valid_rows: 
            self.finish_with_error(
                "[rmsdfilter] All models have NaN RMSD - alignment failed for every model."
            )

        if not filtered:
            self.finish_with_error(
                f"[rmsdfilter] With threshold {self.params['threshold']:.3f} Å, "
                "ALL models were filtered out!"
            )

        self.log(
            f"with threshold {self.params['threshold']:.3f} Å:"
            f"{percent_filtered:6.2f}% of models were filtered out, {len(filtered)} model(s) passed."
        )

        self.output_models = filtered
        self.export_io_models()
