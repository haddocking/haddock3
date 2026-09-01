"""The molecular systems and base workflows the corpus is built from.

Base runs are the main cost lever: a *small* set, each serving many test
cases.  Systems are chosen small enough that a ``topoaa``/``rigidbody`` job
lands well under its 10 s ceiling, which is a prerequisite for Gate 2 rather
than a detail.

The corpus needs at least two sizes.  A **tiny** run absorbs the perturbations
that invalidate *everything* -- a toppar file changed, ``iniseed`` changed --
which on a large run would produce hundreds of misses and be unaffordable.  A
**larger** run reaching a clustering module with enough models for clustering
to be meaningful is what Axis 5 requires.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from .config import Config, Step

REPO_ROOT = Path(__file__).resolve().parents[3]
EXAMPLES = REPO_ROOT / "examples"


@dataclass(frozen=True)
class System:
    """Input files copied into the corpus so perturbations may edit them."""

    name: str
    source: Path
    files: tuple[str, ...]


SYSTEMS: dict[str, System] = {
    system.name: system
    for system in (
        System(
            "glycan",
            EXAMPLES / "docking-protein-glycan" / "data",
            ("1LMQ_r_u.pdb", "1LMQ_l_u.pdb", "ambig.tbl", "target.pdb"),
        ),
        System(
            "pp",
            EXAMPLES / "docking-protein-protein" / "data",
            (
                "e2aP_1F3G.pdb",
                "hpr_ensemble.pdb",
                "e2a-hpr_air.tbl",
                "e2a-hpr_1GGR.pdb",
            ),
        ),
        System(
            "scoring",
            EXAMPLES / "scoring" / "data",
            ("protein-protein_1w.pdb", "protein-protein_2w.pdb"),
        ),
        System(
            "cg",
            EXAMPLES / "docking-protein-protein-shape" / "data",
            (
                "2r15_A.pdb",
                "2r15_B.pdb",
                "shape.pdb",
                "ambig-shape.tbl",
                "2r15_reference.pdb",
            ),
        ),
    )
}


@dataclass(frozen=True)
class BaseRunSpec:
    """One base run: a system, a workflow, and why it exists."""

    name: str
    system: str
    purpose: str
    #: Rough expected wall time, for the generator's progress report only.
    cost: str = "fast"
    tags: tuple[str, ...] = field(default_factory=tuple)


def _top(run_dir: Path, molecules: list[Path], ncores: int = 8) -> dict:
    return {
        "run_dir": str(run_dir),
        "mode": "local",
        "ncores": ncores,
        # The main corpus keeps all three off: they rewrite or compress
        # outputs and thereby destroy Gate 1.  Axis 3 gives them dedicated
        # cases asserted on Gate 2 and decompressed content instead.
        "clean": False,
        "postprocess": False,
        "gen_archive": False,
        "molecules": [str(molecule) for molecule in molecules],
    }


def build_config(name: str, run_dir: Path, data: dict[str, Path]) -> Config:
    """Return the configuration for base run ``name``."""
    builder = _BUILDERS[name]
    return builder(run_dir, data)


def _glycan_molecules(data: dict[str, Path]) -> list[Path]:
    return [data["1LMQ_r_u.pdb"], data["1LMQ_l_u.pdb"]]


def _tiny(run_dir: Path, data: dict[str, Path]) -> Config:
    return Config(
        top=_top(run_dir, _glycan_molecules(data), ncores=4),
        steps=[
            Step("topoaa"),
            Step(
                "rigidbody",
                {
                    "tolerance": 20,
                    "ambig_fname": str(data["ambig.tbl"]),
                    "sampling": 4,
                    "w_vdw": 1,
                },
            ),
        ],
    )


def _tiny_wide(run_dir: Path, data: dict[str, Path]) -> Config:
    """The tiny workflow with many more sampling jobs and nothing else changed.

    Exists solely so that ``t_hit`` can be measured rather than guessed: the
    only difference from ``tiny`` is the number of CNS jobs, so the marginal
    wall time per job between the two all-hit reruns is the per-job cost of
    cache resolution, with startup and analysis overhead cancelling out.
    """
    config = _tiny(run_dir, data)
    config.find("rigidbody").params["sampling"] = 40
    return config


def _tiny_topoaa(run_dir: Path, data: dict[str, Path]) -> Config:
    return Config(
        top=_top(run_dir, _glycan_molecules(data), ncores=4),
        steps=[Step("topoaa")],
    )


def _refine(run_dir: Path, data: dict[str, Path]) -> Config:
    restraints = str(data["ambig.tbl"])
    return Config(
        top=_top(run_dir, _glycan_molecules(data), ncores=4),
        steps=[
            Step("topoaa"),
            Step(
                "rigidbody",
                {
                    "tolerance": 20,
                    "ambig_fname": restraints,
                    "sampling": 2,
                    "w_vdw": 1,
                },
            ),
            Step("caprieval", {"reference_fname": str(data["target.pdb"])}),
            # An analysis module that only observes the models. Its presence
            # is what makes Axis 2.4 -- reorder two independent modules --
            # expressible at all.
            Step("rmsdmatrix"),
            Step("flexref", {"tolerance": 20, "ambig_fname": restraints}),
            Step("emref", {"tolerance": 20, "ambig_fname": restraints}),
            Step("mdref", {"tolerance": 20, "ambig_fname": restraints}),
        ],
    )


def _pp_topoaa_subsections() -> dict:
    return {
        "mol1": {"nhisd": 0, "nhise": 1, "hise_1": 75},
        "mol2": {"nhisd": 1, "hisd_1": 76, "nhise": 1, "hise_1": 15},
    }


def _pp(run_dir: Path, data: dict[str, Path]) -> Config:
    restraints = str(data["e2a-hpr_air.tbl"])
    reference = str(data["e2a-hpr_1GGR.pdb"])
    molecules = [data["e2aP_1F3G.pdb"], data["hpr_ensemble.pdb"]]
    return Config(
        top=_top(run_dir, molecules, ncores=12),
        steps=[
            Step("topoaa", {"autohis": False}, _pp_topoaa_subsections()),
            Step(
                "rigidbody",
                {"tolerance": 20, "ambig_fname": restraints, "sampling": 40},
            ),
            Step("caprieval", {"reference_fname": reference}),
            Step("seletop", {"select": 5}),
            Step("flexref", {"tolerance": 20, "ambig_fname": restraints}),
            Step("caprieval", {"reference_fname": reference}),
            Step("emref", {"tolerance": 20, "ambig_fname": restraints}),
            Step("clustfcc", {"min_population": 1}),
            Step("seletopclusts", {"top_models": 4}),
            Step("caprieval", {"reference_fname": reference}),
        ],
    )


def _pp_cluster(run_dir: Path, data: dict[str, Path]) -> Config:
    restraints = str(data["e2a-hpr_air.tbl"])
    reference = str(data["e2a-hpr_1GGR.pdb"])
    molecules = [data["e2aP_1F3G.pdb"], data["hpr_ensemble.pdb"]]
    return Config(
        top=_top(run_dir, molecules, ncores=12),
        steps=[
            Step("topoaa", {"autohis": False}, _pp_topoaa_subsections()),
            Step(
                "rigidbody",
                {"tolerance": 20, "ambig_fname": restraints, "sampling": 40},
            ),
            Step("caprieval", {"reference_fname": reference}),
            Step("clustfcc", {"min_population": 2}),
            Step("seletopclusts", {"top_clusters": 3, "top_models": 2}),
            Step("flexref", {"tolerance": 20, "ambig_fname": restraints}),
            Step("emref", {"tolerance": 20, "ambig_fname": restraints}),
            Step("caprieval", {"reference_fname": reference}),
        ],
    )


def _scoring(run_dir: Path, data: dict[str, Path]) -> Config:
    molecules = [data["protein-protein_1w.pdb"], data["protein-protein_2w.pdb"]]
    return Config(
        top=_top(run_dir, molecules, ncores=4),
        steps=[
            Step("topoaa", {"autohis": True}),
            Step("emscoring", {"tolerance": 20}),
            Step("mdscoring", {"tolerance": 20}),
        ],
    )


def _cg(run_dir: Path, data: dict[str, Path]) -> Config:
    restraints = str(data["ambig-shape.tbl"])
    reference = str(data["2r15_reference.pdb"])
    molecules = [data["2r15_A.pdb"], data["2r15_B.pdb"], data["shape.pdb"]]
    shape = {"mol_shape_3": True, "mol_fix_origin_3": True}
    return Config(
        top=_top(run_dir, molecules, ncores=8),
        steps=[
            Step("topoaa"),
            Step("topocg"),
            Step(
                "rigidbody",
                {"tolerance": 20, "ambig_fname": restraints, "sampling": 4, **shape},
            ),
            Step("seletop", {"select": 2}),
            Step("flexref", {"tolerance": 20, "ambig_fname": restraints, **shape}),
            Step("cgtoaa", {"mol_shape_3": True}),
            Step(
                "emref",
                {"tolerance": 20, "ambig_fname": restraints, **shape},
            ),
            Step(
                "caprieval",
                {"reference_fname": reference, "fnat_cutoff": 7.0},
            ),
        ],
    )


_BUILDERS = {
    "tiny": _tiny,
    "tiny-wide": _tiny_wide,
    "tiny-topoaa": _tiny_topoaa,
    "refine": _refine,
    "pp": _pp,
    "pp-cluster": _pp_cluster,
    "scoring": _scoring,
    "cg": _cg,
}


BASE_RUNS: tuple[BaseRunSpec, ...] = (
    BaseRunSpec(
        "tiny",
        "glycan",
        "Global-invalidation probes (toppar, iniseed, canonicalization) and the "
        "cheap Axis 1/2/3/6/7/8/11/12 material. Small enough that a "
        "perturbation invalidating everything costs a handful of misses.",
        cost="fast",
        tags=("topoaa", "rigidbody"),
    ),
    BaseRunSpec(
        "tiny-wide",
        "glycan",
        "The tiny workflow with 40 sampling jobs instead of 4. Its only "
        "purpose is to make t_hit measurable: subtracting the two all-hit "
        "reruns cancels startup and analysis overhead and leaves the per-job "
        "cost of cache resolution.",
        cost="medium",
        tags=("topoaa", "rigidbody"),
    ),
    BaseRunSpec(
        "tiny-topoaa",
        "glycan",
        "A topology-only run. Half of the Axis 11.2 disjoint-coverage pair, and "
        "an Axis 11.5 superset-workflow source.",
        cost="fast",
        tags=("topoaa",),
    ),
    BaseRunSpec(
        "refine",
        "glycan",
        "The flexref/emref/mdref job shapes on the tiny system, so Axis 13 does "
        "not have to be paid for at protein-protein scale.",
        cost="medium",
        tags=("topoaa", "rigidbody", "flexref", "emref", "mdref"),
    ),
    BaseRunSpec(
        "pp",
        "pp",
        "The seletop selection route, and a ten-step workflow -- one step below "
        "the zero-fill boundary, so Axis 2.5 is a single added module away.",
        cost="medium",
        tags=("topoaa", "rigidbody", "flexref", "emref"),
    ),
    BaseRunSpec(
        "pp-cluster",
        "pp",
        "The Axis 5 workhorse: clustering with enough models for the cluster "
        "membership and rank of the refined set to actually move.",
        cost="slow",
        tags=("topoaa", "rigidbody", "flexref", "emref"),
    ),
    BaseRunSpec(
        "scoring",
        "scoring",
        "The emscoring/mdscoring job shape -- scoring rather than sampling.",
        cost="medium",
        tags=("topoaa", "emscoring", "mdscoring"),
    ),
    BaseRunSpec(
        "cg",
        "cg",
        "The topocg and cgtoaa job shapes: coarse-grained topology and the "
        "coarse-grained to all-atom conversion.",
        cost="medium",
        tags=("topoaa", "topocg", "rigidbody", "flexref", "cgtoaa", "emref"),
    ),
)
