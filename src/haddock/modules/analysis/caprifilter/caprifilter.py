"""Helper functions for the caprifilter module."""

import os
from math import isnan
from pathlib import Path

from haddock.libs.libontology import PDBFile


# All metric attribute names as they appear on the CAPRI class.
METRIC_NAMES = ("irmsd", "lrmsd", "ilrmsd", "fnat", "dockq", "rmsd")

# Mapping from filter_by to metric name, so relevant computations are performed in libcapri.
_METRIC_TO_PARAM: dict[str, str] = {
    "irmsd": "irmsd",
    "lrmsd": "lrmsd",
    "ilrmsd": "ilrmsd",
    "fnat": "fnat",
    "dockq": "dockq",
    "rmsd": "global_rmsd",
}

# Metrics that dockq computation depends on.
_DOCKQ_DEPS: tuple[str, ...] = ("fnat", "irmsd", "lrmsd")

# Public set of valid metric names accepted in filter_by.
VALID_FILTER_METRICS: frozenset[str] = frozenset(_METRIC_TO_PARAM)


def get_capri_params(filter_by: list[str]) -> dict[str, bool]:
    """Return CAPRI computation params required for the given filter metrics.

    Translates the user-facing ``filter_by`` metric names into the boolean
    computation flags expected by the CAPRI class. All six computation params
    are returned: those required by ``filter_by`` are set to True, the rest
    to False. If ``dockq`` is requested, its dependencies (fnat, irmsd, lrmsd)
    are also enabled automatically.

    Parameters
    ----------
    filter_by : list[str]
        List of metrics to filter on. Valid values: irmsd, lrmsd, ilrmsd,
        fnat, dockq, rmsd.

    Returns
    -------
    params : dict[str, bool]
        Computation param flags ready to be injected into ``self.params``
        before CAPRI jobs are created. Keys: irmsd, lrmsd, ilrmsd, fnat,
        dockq, global_rmsd.
    """
    params: dict[str, bool] = {p: False for p in _METRIC_TO_PARAM.values()}
    for metric in filter_by:
        params[_METRIC_TO_PARAM[metric]] = True
        if metric == "dockq":
            for dep in _DOCKQ_DEPS:
                params[dep] = True
    return params


def collect_metrics(capri_objects: list) -> dict[PDBFile, dict[str, float]]:
    """Return {model: {metric: value}} for all CAPRI metrics.

    Parameters
    ----------
    capri_objects : list[CAPRI]
        List of CAPRI objects, one per model
        (if multiple references given - one best reference is already selected).

    Returns
    -------
    metrics_data : dict[PDBFile, dict[str, float]]
        Dictionary mapping each model's PDBFile to a dict of metric values.
    """
    return {
        capri.model: {m: getattr(capri, m) for m in METRIC_NAMES}
        for capri in capri_objects
    }


def filter_models(
    metrics_data: dict[PDBFile, dict[str, float]],
    filter_specs: dict[str, tuple[float, str]],
) -> tuple[list[PDBFile], list[PDBFile]]:
    """Split models into kept and filtered_out based on multiple metric filters.

    All filters are applied simultaneously with AND logic. A model must pass
    every active filter to be kept. Models with NaN for any filtered metric are
    always removed.

    Parameters
    ----------
    metrics_data : dict[PDBFile, dict[str, float]]
        Metric values per model, as returned by collect_metrics.
    filter_specs : dict[str, tuple[float, str]]
        {metric: (cutoff, filter_out)} where filter_out is 'above' or 'below'.
        'above': filter out models with value > cutoff (keep value <= cutoff).
        'below': filter out models with value < cutoff (keep value >= cutoff).

    Returns
    -------
    kept, filtered_out : tuple[list[PDBFile], list[PDBFile]]
        'filtered_out' is not used in the caprifilter module, it is kept for tests.
    """
    kept: list[PDBFile] = []
    filtered_out: list[PDBFile] = []

    for model, metric_vals in metrics_data.items():
        passes = True
        for metric, (cutoff, filter_out) in filter_specs.items():
            value = metric_vals.get(metric, float("nan"))
            if isnan(value):
                passes = False
                break
            if filter_out == "above" and value > cutoff:
                passes = False
                break
            if filter_out == "below" and value < cutoff:
                passes = False
                break
        (kept if passes else filtered_out).append(model)

    return kept, filtered_out


def write_caprifilter(
    kept: list[PDBFile],
    capri_objects: list,
    filter_specs: dict[str, tuple[float, str]],
    fname: str = "caprifilter.tsv",
) -> None:
    """Write TSV with kept models and only the user-requested metric columns.

    Parameters
    ----------
    kept : list[PDBFile]
        Models that passed all filters.
    capri_objects : list[CAPRI]
        CAPRI objects (one per model, best reference already selected).
    filter_specs : dict[str, tuple[float, str]]
        Active filter specifications {metric: (cutoff, filter_out)}.
    fname : str
        Output file name.
    """
    capri_by_model: dict[int, object] = {
        id(capri.model): capri for capri in capri_objects
    }
    metrics = list(filter_specs.keys())

    filter_parts = []
    for metric, (cutoff, filter_out) in filter_specs.items():
        op = ">" if filter_out == "above" else "<"
        filter_parts.append(f"{metric}{op}{cutoff:.3f}")
    filter_summary = ", ".join(filter_parts) if filter_parts else "none"

    def fmt(v):
        return (
            "nan" if (v is None or (isinstance(v, float) and isnan(v))) else f"{v:.3f}"
        )

    with open(fname, "w") as fh:
        fh.write(
            f"# caprifilter: filter=[{filter_summary}]; "
            f"{len(kept)} model(s) kept."
            f"{os.linesep}"
        )
        fh.write("\t".join(["model", "score"] + metrics) + os.linesep)

        for model in kept:
            capri = capri_by_model.get(id(model))
            score_str = fmt(capri.score if capri is not None else float("nan"))
            metric_strs = [
                fmt(getattr(capri, m) if capri is not None else float("nan"))
                for m in metrics
            ]
            fh.write(
                "\t".join([str(model.rel_path), score_str] + metric_strs) + os.linesep
            )


def write_caprifilter_full(
    capri_objects: list,
    metrics_data: dict[PDBFile, dict[str, float]],
    kept: list[PDBFile],
    filter_specs: dict[str, tuple[float, str]],
    fname: str = "caprifilter_all_models.tsv",
) -> None:
    """Write TSV with all models, user-requested metric columns, and a status column.

    Parameters
    ----------
    capri_objects : list[CAPRI]
        CAPRI objects (one per model, best reference already selected).
    metrics_data : dict[PDBFile, dict[str, float]]
        Metric values per model, as returned by collect_metrics.
    kept : list[PDBFile]
        Models that passed all filters.
    filter_specs : dict[str, tuple[float, str]]
        Active filter specifications {metric: (cutoff, filter_out)}.
    fname : str
        Output file name.
    """
    kept_set = set(id(m) for m in kept)
    metrics = list(filter_specs.keys())

    # index CAPRI objects by model identity for score lookup
    capri_by_model: dict[int, object] = {
        id(capri.model): capri for capri in capri_objects
    }

    n_kept = len(kept)
    n_total = len(metrics_data)
    pct_filtered = (1 - n_kept / n_total) * 100 if n_total else 0.0

    filter_parts = []
    for metric, (cutoff, filter_out) in filter_specs.items():
        op = ">" if filter_out == "above" else "<"
        filter_parts.append(f"{metric}{op}{cutoff:.3f}")
    filter_summary = ", ".join(filter_parts) if filter_parts else "none"

    def fmt(v):
        return (
            "nan" if (v is None or (isinstance(v, float) and isnan(v))) else f"{v:.3f}"
        )

    with open(fname, "w") as fh:
        fh.write(
            f"# caprifilter: filter=[{filter_summary}]; "
            f"{n_kept} model(s) kept; {pct_filtered:.2f}% filtered out."
            f"{os.linesep}"
        )
        fh.write("\t".join(["model", "score"] + metrics + ["status"]) + os.linesep)

        for model, metric_vals in metrics_data.items():
            capri = capri_by_model.get(id(model))
            score_str = fmt(capri.score if capri is not None else float("nan"))
            status = "kept" if id(model) in kept_set else "filtered"
            metric_strs = [fmt(metric_vals[m]) for m in metrics]
            fh.write(
                "\t".join([str(model.rel_path), score_str] + metric_strs + [status])
                + os.linesep
            )
