"""Load and execute Phase 2 cases.

A case is data, not code::

    - case: axis5.3-rank-changed
      taxonomy: "5.3"
      title: Same structure, different rank, same filename
      why: >
        A refinement job takes one input model, so that model binds to the same
        pin whatever its rank. Rank is a pure locator and the job must hit.
      base: pp-cluster
      sources: [pp-cluster]
      mode: complete
      perturb:
        - {op: module, module: clustfcc, set: {clust_cutoff: 0.65}}
      expect: {default: auto}

The harness runs ``haddock3 <config> --cache <s1> --cache <s2> ...`` as a
subprocess with the given environment, measures wall time, then ``stat()``s
each declared pair.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from . import thresholds
from .corpus import Manifest
from .harness import RunResult, run_haddock3
from .perturb import Context, apply_all, cosmetic_text
from .systems import BASE_RUNS, SYSTEMS, build_config

CASES_DIR = Path(__file__).resolve().parents[1] / "cases"


@dataclass
class Case:
    """One declarative test case."""

    case: str
    axis: str = ""
    taxonomy: str = ""
    title: str = ""
    why: str = ""
    base: str = "tiny"
    sources: list[str] = field(default_factory=list)
    mode: str = "complete"
    perturb: list[dict] = field(default_factory=list)
    expect: dict = field(default_factory=dict)
    env: dict = field(default_factory=dict)
    misses: dict | None = None
    overhead: float | None = None
    floor: float | None = None
    calibrate: str | None = None
    witness: str | None = None
    #: An explicitly recorded scope boundary.  Never a silent skip.
    skip: str = ""
    #: Requires a second filesystem.
    needs_other_filesystem: bool = False
    #: Gate 1 is blind for this case; only Gate 2 and Gate 3 remain.
    degraded: str = ""
    #: A known defect *outside* caching that this case runs into. Recorded so
    #: a failure here is not read as a caching finding. It never causes a
    #: skip: the case still runs, and still reports.
    known_unrelated: str = ""
    expect_failure: bool = False
    #: Axis 5 anti-vacuity: at least this many hits must come from a source
    #: file with a *different* name, or the perturbation did not reorder
    #: anything and the case proves nothing.
    require_renamed_hits: int = 0
    #: An everyday case: something a user would plausibly do on purpose.
    #: These are the ones worth reading first, and the ones ``try_case.py``
    #: offers under ``--common``.
    common: bool = False
    source_file: str = ""

    @property
    def identifier(self) -> str:
        return self.case


def load_cases(directory: Path | None = None) -> list[Case]:
    """Read every ``cases/*.yaml`` file, in axis order."""
    directory = directory or CASES_DIR
    cases: list[Case] = []
    for path in sorted(directory.glob("*.yaml")):
        document = yaml.safe_load(path.read_text(encoding="utf-8")) or []
        for entry in document:
            entry = dict(entry)
            entry.setdefault("axis", path.stem.replace("axis", ""))
            entry["source_file"] = path.name
            cases.append(Case(**entry))
    names = [case.case for case in cases]
    duplicates = {name for name in names if names.count(name) > 1}
    if duplicates:
        raise ValueError(f"duplicate case names: {sorted(duplicates)}")
    return cases


#: Perturbation operations that write to the input directory.  Only these
#: justify giving B a private copy of the inputs.
_INPUT_OPS = frozenset(
    {"edit_input", "copy_input", "move_input", "restraints_archive"}
)


def _touches_inputs(case: "Case") -> bool:
    return any(operation.get("op") in _INPUT_OPS for operation in case.perturb)


@dataclass
class Prepared:
    """A case, materialized into a config, an environment and a source list."""

    case: Case
    config_path: Path
    run_dir: Path
    #: Where B's input files live -- the corpus's own directory unless a
    #: perturbation edits an input, in which case a private copy.
    data_dir: Path
    sources: dict[str, Path]
    env: dict[str, str]
    install: Path | None
    cwd: Path
    notes: list[str]
    overhead: float
    #: How long the base run took to build from scratch, from the manifest.
    #: A case is a variation on that workflow, so this bounds how long it can
    #: legitimately take even when every single job misses.
    base_duration: float = 0.0


def prepare(case: Case, corpus: Manifest, corpus_dir: Path, tmp: Path) -> Prepared:
    """Build B: its configuration, its inputs, and its environment."""
    spec = next(entry for entry in BASE_RUNS if entry.name == case.base)
    system = SYSTEMS[spec.system]

    # B must differ from A by exactly one controlled perturbation. Copying the
    # inputs to a private directory would be a *second*, uncontrolled one --
    # it moves every input file, which is itself Axis 1.6 -- so the corpus's
    # own input directory is used unless a perturbation actually edits an
    # input, and only then is a private copy made.
    corpus_data = corpus_dir / "systems" / system.name
    if _touches_inputs(case):
        data_dir = tmp / "data"
        data_dir.mkdir(parents=True, exist_ok=True)
        data = {}
        for name in system.files:
            destination = data_dir / name
            shutil.copy2(corpus_data / name, destination)
            data[name] = destination
    else:
        data_dir = corpus_data
        data = {name: corpus_data / name for name in system.files}

    run_dir = tmp / "run"
    config = build_config(case.base, run_dir, data)
    context = Context(
        config=config,
        data=data,
        data_dir=data_dir,
        tmp=tmp,
        env=dict(case.env),
        cwd=tmp,
    )
    apply_all(case.perturb, context)
    config_path = context.config.write(tmp / f"{case.case}.cfg")
    if any(operation.get("op") == "cosmetic" for operation in case.perturb):
        config_path.write_text(
            cosmetic_text(config_path.read_text(encoding="utf-8")), encoding="utf-8"
        )

    sources = {}
    for name in case.sources:
        fixture = corpus.require(name)
        if not fixture.path:
            raise FileNotFoundError(
                f"corpus fixture {name!r} is unusable: {fixture.notes}"
            )
        sources[name] = fixture.run_dir(corpus_dir)

    base_fixture = corpus.fixtures.get(case.base)
    overhead = case.overhead
    if overhead is None:
        overhead = (
            base_fixture.allhit
            if base_fixture and base_fixture.allhit
            else thresholds.DEFAULT_OVERHEAD
        )

    env = {"HADDOCK_CACHE_HARDLINK": "1"}
    env.update(context.env)
    return Prepared(
        case=case,
        config_path=config_path,
        run_dir=Path(str(context.config.top["run_dir"])),
        data_dir=context.data_dir,
        sources=sources,
        env=env,
        install=context.install,
        cwd=context.cwd or tmp,
        notes=context.notes,
        overhead=float(overhead),
        base_duration=float(base_fixture.duration if base_fixture else 0.0),
    )


def execute(prepared: Prepared, timeout: float | None) -> RunResult:
    """Invoke ``haddock3`` for a prepared case."""
    return run_haddock3(
        prepared.config_path,
        cache=list(prepared.sources.values()),
        env=prepared.env,
        cwd=prepared.cwd,
        timeout=timeout,
        install=prepared.install,
    )
