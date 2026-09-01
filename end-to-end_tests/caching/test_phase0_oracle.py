"""Phase 0 -- does the oracle identify a source entry by content?

Part of the instrument, not of the feature: nothing here runs ``haddock3``.

Every case in this suite is judged against an expected mapping, and that
mapping is computed by looking for the job that produced each output in each
source run. If that search identifies a job by its *name* -- same module, same
occurrence, same filename -- then the suite asks the question the feature is
forbidden to ask, and the consequences are not symmetric:

* a correct implementation fails cases it should pass, wherever a job moved to
  a different name;
* an incorrect one passes them, and one of those cases can only be passed by
  serving a result from the wrong entry -- the single failure this suite
  exists to catch.

So the search itself is checked here, against small run directories written by
hand, before it is used to judge anything.
"""

import json
from pathlib import Path

import pytest

from cachesuite.expectations import ExpectationError, _find_sources
from cachesuite.harness import Artifact

pytestmark = pytest.mark.phase0

PDBFILE = "haddock.libs.libontology.PDBFile"
TOPOLOGYFILE = "haddock.libs.libontology.TopologyFile"


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _topology(folder: Path, stem: str) -> dict:
    return {
        "py/object": TOPOLOGYFILE,
        "file_name": f"{stem}.psf",
        "path": str(folder),
    }


def _model(folder: Path, name: str, **fields) -> dict:
    return {
        "py/object": PDBFILE,
        "file_name": name,
        "path": str(folder),
        "score": 0.0,
        **fields,
    }


def _run(
    root: Path,
    *,
    member: str,
    topology_name: str,
    sampling_name: str,
    seed: int,
) -> Path:
    """A two-step run: one topology job, one sampling job that docks it."""
    topoaa = root / "0_topoaa"
    rigidbody = root / "1_rigidbody"

    # The structure the topology job was built from, and its outputs.
    _write(topoaa / "member.pdb", member)
    _write(topoaa / f"{topology_name}.pdb", "ATOM\nEND\n")
    _write(topoaa / f"{topology_name}.psf", "PSF\n")
    topology_output = _model(
        topoaa,
        f"{topology_name}.pdb",
        ori_name="member",
        topology=_topology(topoaa, topology_name),
    )
    _write(topoaa / "io.json", json.dumps({"input": [], "output": [topology_output]}))

    # The sampling job that docked it.
    _write(rigidbody / f"{sampling_name}.pdb", "ATOM\nEND\n")
    sampling_output = _model(
        rigidbody,
        f"{sampling_name}.pdb",
        seed=seed,
        topology=[_topology(topoaa, topology_name)],
    )
    _write(
        rigidbody / "io.json",
        json.dumps({"input": [topology_output], "output": [sampling_output]}),
    )
    return root


def _artifact(module: str, folder: str, name: str, kind: str = "pdb") -> Artifact:
    return Artifact(f"{folder}/{name}", module, 0, 0, kind)


def test_a_topology_job_is_found_under_another_name(tmp_path):
    """The output moved to a different filename; the job did not move."""
    source = _run(
        tmp_path / "a",
        member="ATOM  the same structure\n",
        topology_name="member_1_haddock",
        sampling_name="rigidbody_1",
        seed=11,
    )
    target = _run(
        tmp_path / "b",
        member="ATOM  the same structure\n",
        topology_name="member_7_haddock",
        sampling_name="rigidbody_1",
        seed=11,
    )

    found = _find_sources(
        _artifact("topoaa", "0_topoaa", "member_7_haddock.pdb"),
        target,
        {"a": source},
    )

    assert found == (source / "0_topoaa" / "member_1_haddock.pdb",)


def test_a_topology_job_reading_other_bytes_is_not_found(tmp_path):
    """Same name, different structure: nothing in the source holds this job."""
    source = _run(
        tmp_path / "a",
        member="ATOM  one structure\n",
        topology_name="member_1_haddock",
        sampling_name="rigidbody_1",
        seed=11,
    )
    target = _run(
        tmp_path / "b",
        member="ATOM  a different structure\n",
        topology_name="member_1_haddock",
        sampling_name="rigidbody_1",
        seed=11,
    )

    assert (
        _find_sources(
            _artifact("topoaa", "0_topoaa", "member_1_haddock.pdb"),
            target,
            {"a": source},
        )
        == ()
    )


def test_a_topology_hit_delivers_the_matching_topology_file(tmp_path):
    """Both outputs of one job come from that one job, not from one name."""
    source = _run(
        tmp_path / "a",
        member="ATOM  the same structure\n",
        topology_name="member_1_haddock",
        sampling_name="rigidbody_1",
        seed=11,
    )
    target = _run(
        tmp_path / "b",
        member="ATOM  the same structure\n",
        topology_name="member_7_haddock",
        sampling_name="rigidbody_1",
        seed=11,
    )

    found = _find_sources(
        _artifact("topoaa", "0_topoaa", "member_7_haddock.psf", kind="psf"),
        target,
        {"a": source},
    )

    assert found == (source / "0_topoaa" / "member_1_haddock.psf",)


def test_a_sampling_job_is_found_by_its_combination_and_seed(tmp_path):
    """A docking job renumbered by a longer schedule is the same job."""
    source = _run(
        tmp_path / "a",
        member="ATOM  the same structure\n",
        topology_name="member_1_haddock",
        sampling_name="rigidbody_1",
        seed=11,
    )
    target = _run(
        tmp_path / "b",
        member="ATOM  the same structure\n",
        topology_name="member_1_haddock",
        sampling_name="rigidbody_9",
        seed=11,
    )

    found = _find_sources(
        _artifact("rigidbody", "1_rigidbody", "rigidbody_9.pdb"),
        target,
        {"a": source},
    )

    assert found == (source / "1_rigidbody" / "rigidbody_1.pdb",)


def test_a_sampling_job_with_another_seed_is_a_different_job(tmp_path):
    """Same inputs, different seed: a repeat, and repeats are new work."""
    source = _run(
        tmp_path / "a",
        member="ATOM  the same structure\n",
        topology_name="member_1_haddock",
        sampling_name="rigidbody_1",
        seed=11,
    )
    target = _run(
        tmp_path / "b",
        member="ATOM  the same structure\n",
        topology_name="member_1_haddock",
        sampling_name="rigidbody_1",
        seed=12,
    )

    assert (
        _find_sources(
            _artifact("rigidbody", "1_rigidbody", "rigidbody_1.pdb"),
            target,
            {"a": source},
        )
        == ()
    )


def test_two_duplicate_jobs_resolve_to_the_one_entry_that_holds_them(tmp_path):
    """A duplicated ensemble member is not new work, and must say so."""
    source = _run(
        tmp_path / "a",
        member="ATOM  the same structure\n",
        topology_name="member_1_haddock",
        sampling_name="rigidbody_1",
        seed=11,
    )
    target = _run(
        tmp_path / "b",
        member="ATOM  the same structure\n",
        topology_name="member_1_haddock",
        sampling_name="rigidbody_1",
        seed=11,
    )
    # A second member of the target's ensemble, byte-identical to the first.
    topoaa = target / "0_topoaa"
    _write(topoaa / "member_copy.pdb", "ATOM  the same structure\n")
    _write(topoaa / "member_11_haddock.pdb", "ATOM\nEND\n")
    document = json.loads((topoaa / "io.json").read_text(encoding="utf-8"))
    document["output"].append(
        _model(
            topoaa,
            "member_11_haddock.pdb",
            ori_name="member_copy",
            topology=_topology(topoaa, "member_11_haddock"),
        )
    )
    _write(topoaa / "io.json", json.dumps(document))

    original = _find_sources(
        _artifact("topoaa", "0_topoaa", "member_1_haddock.pdb"), target, {"a": source}
    )
    duplicate = _find_sources(
        _artifact("topoaa", "0_topoaa", "member_11_haddock.pdb"), target, {"a": source}
    )

    assert duplicate == original == (source / "0_topoaa" / "member_1_haddock.pdb",)


def test_a_source_without_the_step_answers_for_nothing(tmp_path):
    """A partial source is a source that holds fewer jobs, not a broken one."""
    source = _run(
        tmp_path / "a",
        member="ATOM  the same structure\n",
        topology_name="member_1_haddock",
        sampling_name="rigidbody_1",
        seed=11,
    )
    (source / "1_rigidbody" / "io.json").unlink()
    (source / "1_rigidbody" / "rigidbody_1.pdb").unlink()
    target = _run(
        tmp_path / "b",
        member="ATOM  the same structure\n",
        topology_name="member_1_haddock",
        sampling_name="rigidbody_1",
        seed=11,
    )

    assert (
        _find_sources(
            _artifact("rigidbody", "1_rigidbody", "rigidbody_1.pdb"),
            target,
            {"a": source},
        )
        == ()
    )


def test_a_sampling_job_is_found_through_jsonpickle_back_references(tmp_path):
    """Real ``io.json`` writes a docked complex's topologies as pointers.

    A docking job carries the very same topology objects its inputs carry, so
    jsonpickle writes them out once and refers back to them afterwards. The
    field that says which models were combined is therefore a set of pointers
    in every real run, and a reader that skipped them would be unable to say
    what any docking job docked -- while still answering confidently for every
    other shape.
    """
    source = _run(
        tmp_path / "a",
        member="ATOM  the same structure\n",
        topology_name="member_1_haddock",
        sampling_name="rigidbody_1",
        seed=11,
    )
    target = _run(
        tmp_path / "b",
        member="ATOM  the same structure\n",
        topology_name="member_1_haddock",
        sampling_name="rigidbody_9",
        seed=11,
    )
    for run in (source, target):
        path = run / "1_rigidbody" / "io.json"
        document = json.loads(path.read_text(encoding="utf-8"))
        # The input keeps the topology written out in full; the output refers
        # back to it, exactly as jsonpickle does.
        document["output"][0]["topology"] = [{"py/id": 6}]
        _write(path, json.dumps(document))

    found = _find_sources(
        _artifact("rigidbody", "1_rigidbody", "rigidbody_9.pdb"),
        target,
        {"a": source},
    )

    assert found == (source / "1_rigidbody" / "rigidbody_1.pdb",)


def test_an_unresolvable_reference_is_reported_rather_than_guessed(tmp_path):
    """More pointers than topologies written out: the reader must not guess."""
    source = _run(
        tmp_path / "a",
        member="ATOM  the same structure\n",
        topology_name="member_1_haddock",
        sampling_name="rigidbody_1",
        seed=11,
    )
    target = _run(
        tmp_path / "b",
        member="ATOM  the same structure\n",
        topology_name="member_1_haddock",
        sampling_name="rigidbody_1",
        seed=11,
    )
    path = target / "1_rigidbody" / "io.json"
    document = json.loads(path.read_text(encoding="utf-8"))
    document["output"][0]["topology"] = [{"py/id": 6}, {"py/id": 9}]
    _write(path, json.dumps(document))

    with pytest.raises(ExpectationError):
        _find_sources(
            _artifact("rigidbody", "1_rigidbody", "rigidbody_1.pdb"),
            target,
            {"a": source},
        )
