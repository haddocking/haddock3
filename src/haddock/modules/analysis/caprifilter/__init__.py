"""Filter models by any combination of CAPRI metrics computed against a reference structure.

The user specifies a list of metrics via ``filter_by`` (e.g. ``['irmsd', 'fnat']``).
For each metric, a per-metric cutoff and filter direction are set independently:

- ``{metric}_filter_cutoff``: the numerical threshold.
- ``{metric}_filter_out``: ``above`` removes models where value > cutoff (keep ≤);
  ``below`` removes models where value < cutoff (keep ≥).

All active filters are applied simultaneously (AND logic), a model must pass
every filter to be kept. The following metrics are supported:
- RMSD
- iRMSD (interface RMSD)
- lRMSD (ligand RMSD)
- ilRMSD (interface ligand RMSD)
- FNAT (fraction of native contacts)
- DOCKQ

The following files are generated:
- **caprifilter.tsv**: kept models with score and the user-requested metric
  columns only.
- **caprifilter_all_models.tsv**: all input models, user-requested metrics,
  and a ``status`` column (``kept`` / ``filtered``). Written only when
  ``caprifilter_full = true``.
- **caprifilter_ss.tsv**: caprieval-style ranked table of all models.
- **caprifilter_multiref.tsv**: when multiple references are provided,
  one row per (model, reference) pair.

For more details about this module, please `refer to the haddock3 user manual
<https://www.bonvinlab.org/haddock3-user-manual/modules/analysis.html#caprifilter-module>`_
"""

from math import isnan
from pathlib import Path

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import FilePath, Union
from haddock.libs.libaa2cg import martinize
from haddock.libs.libcapri import (
    CAPRI,
    extract_data_from_capri_class,
    extract_models_best_references,
)
from haddock.libs.libontology import PDBFile
from haddock.libs.libparallel import Scheduler
from haddock.libs.libpdb import handle_input_reference
from haddock.libs.libstructure import find_ff
from haddock.modules import BaseHaddockModule
from haddock.modules.analysis.caprifilter.caprifilter import (
    VALID_FILTER_METRICS,
    collect_metrics,
    filter_models,
    get_capri_params,
    write_caprifilter,
    write_caprifilter_full,
)


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module to filter models by CAPRI metrics."""

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
        # Get models from previous step
        models = self.previous_io.retrieve_models(individualize=True)
        if self.is_nested(models):
            raise ValueError(
                "[caprifilter] cannot be executed after modules that produce "
                "a nested list of models."
            )

        if not self.params["reference_fname"]:
            self.finish_with_error(
                "[caprifilter] No reference structure provided! "
                "Please set 'reference_fname' to a valid PDB file."
            )

        filter_by: list[str] = self.params["filter_by"]
        if not filter_by:
            self.finish_with_error(
                "[caprifilter] filter_by is empty, please specify at least one metric."
            )

        unknown = [m for m in filter_by if m not in VALID_FILTER_METRICS]
        if unknown:
            self.finish_with_error(
                f"[caprifilter] Unknown metric(s) in filter_by: {unknown}. "
                f"Valid choices are: {sorted(VALID_FILTER_METRICS)}."
            )

        self.params.update(get_capri_params(filter_by))

        # Build filter_specs: {metric: (cutoff, filter_out)}
        filter_specs: dict[str, tuple[float, str]] = {
            metric: (
                self.params[f"{metric}_filter_cutoff"],
                self.params[f"{metric}_filter_out"],
            )
            for metric in filter_by
        }

        # Build reference list
        reference = Path(self.params["reference_fname"])
        references = handle_input_reference(reference)

        # Handle coarse-grain force field
        ff = find_ff(models)
        if ff == "martini2":
            references = [
                Path(martinize(ref, self.path.resolve().parent, False))
                for ref in references
            ]

        # Create one CAPRI job per (model, reference) pair
        jobs: list[CAPRI] = []
        for i, model in enumerate(models, start=1):
            for ref_id, ref in enumerate(references, start=1):
                jobs.append(
                    CAPRI(
                        identificator=i,
                        model=model,
                        path=Path("."),
                        reference=ref,
                        params=self.params,
                        ref_id=ref_id,
                        ff=ff,
                    )
                )

        engine = Scheduler(
            tasks=jobs,
            ncores=self.params["ncores"],
            max_cpus=self.params["max_cpus"],
        )
        engine.run()

        jobs = engine.results
        jobs = sorted(jobs, key=lambda c: c.identificator)

        # Best reference per model
        best_ref_jobs = extract_models_best_references(jobs)

        # Warn if models carry cluster assignments — filtering ignores them,
        # so surviving models may have stale clt_id/clt_rank values.
        if any(getattr(m, "clt_id", None) is not None for m in models):
            self.log(
                "Warning: input models have cluster assignments (clt_id/clt_rank). "
                "[caprifilter] filters on a per-model basis and does not update cluster info. "
                "Remaining models will retain their original cluster labels."
            )

        # Write full metrics table (same format as caprieval)
        extract_data_from_capri_class(
            capri_objects=best_ref_jobs,
            output_fname=Path(".", "caprifilter_ss.tsv"),
            sort_key=self.params["sortby"],
            sort_ascending=self.params["sort_ascending"],
            add_reference_id=len(references) > 1,
        )

        if len(references) > 1:
            extract_data_from_capri_class(
                capri_objects=jobs,
                output_fname=Path(".", "caprifilter_multiref.tsv"),
                sort_key=self.params["sortby"],
                sort_ascending=self.params["sort_ascending"],
                add_reference_id=True,
            )

        # Collect metrics and apply all filters
        metrics_data = collect_metrics(best_ref_jobs)
        kept, _ = filter_models(metrics_data, filter_specs)

        # Log NaN counts per filtered metric
        for metric in filter_by:
            n_nan = sum(1 for vals in metrics_data.values() if isnan(vals[metric]))
            if n_nan:
                self.log(
                    f"{100 * n_nan / len(models):6.2f}% of models had NaN {metric} "
                    "(alignment failed) and will be excluded."
                )

        if not metrics_data or all(
            all(isnan(v) for v in vals.values()) for vals in metrics_data.values()
        ):
            self.finish_with_error(
                "[caprifilter] All models have NaN metrics, i.e. alignment failed for every model."
            )

        if not kept:
            specs_str = ", ".join(
                f"{m} filter_out={fo} cutoff={c:.3f}"
                for m, (c, fo) in filter_specs.items()
            )
            self.finish_with_error(
                f"[caprifilter] With filters [{specs_str}], "
                "ALL models were filtered out!!"
            )

        pct_filtered = (1 - len(kept) / len(models)) * 100
        specs_str = ", ".join(
            f"{m} filter_out={fo} cutoff={c:.3f}" for m, (c, fo) in filter_specs.items()
        )
        self.log(
            f"Filters: [{specs_str}] — "
            f"{pct_filtered:.2f}% filtered out, {len(kept)} model(s) passed."
        )

        # Write clean output: kept models, requested metric columns only
        write_caprifilter(
            kept=kept,
            capri_objects=best_ref_jobs,
            filter_specs=filter_specs,
        )

        # Write full status table if requested
        if self.params["caprifilter_full"]:
            write_caprifilter_full(
                capri_objects=best_ref_jobs,
                metrics_data=metrics_data,
                kept=kept,
                filter_specs=filter_specs,
            )

        # Send models to the next step, no operation is done on them
        self.output_models = kept  # type: ignore # ignore this here only if we are checking the return type of `retrieve_models` is not nested!!
        self.export_io_models()
