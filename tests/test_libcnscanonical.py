"""Tests for CNS canonical input representations."""

import gzip
import re
from pathlib import Path

import pytest

from haddock.libs.libcnscanonical import (
    build_canonical_mapping,
    compression_transparent_checksum,
)
from haddock.libs.libsubprocess import CNSJob


def test_canonical_mapping_is_independent_of_run_step_and_install_names(tmp_path):
    first, _ = _mapping(tmp_path, "one", "2_rigidbody", "install-a")
    second, _ = _mapping(tmp_path, "two", "02_rigidbody", "install-b")

    assert first.canonical_script == second.canonical_script
    assert first.checksums == second.checksums
    assert first.invariant_dependencies == (
        "canonical-cns",
        "module/protocol.cns",
        "toppar/protein.top",
    )


def test_canonical_mapping_is_independent_of_the_cns_executable(tmp_path):
    """The executable evaluates the computation; it is not part of it.

    Two installations hold different binaries, always -- nobody compiles or
    downloads the same bytes twice -- so an identity that read the
    executable's content could never be shared between them.  The pin is
    still there and still occupies its position; it is bound to a policy
    constant rather than to bytes.
    """
    first, _ = _mapping(tmp_path, "one", "2_rigidbody", "install-a")
    second, _ = _mapping(tmp_path, "two", "2_rigidbody", "install-b")

    assert compression_transparent_checksum(
        first.cns_exec
    ) != compression_transparent_checksum(second.cns_exec)
    assert first.checksums["canonical-cns"] == second.checksums["canonical-cns"]
    assert "canonical-cns" in first.invariant_dependencies


def test_canonical_mapping_is_independent_of_input_and_output_filenames(tmp_path):
    first = _named_model_mapping(tmp_path, "run-a", "rank_1.pdb", "flexref_1.pdb")
    second = _named_model_mapping(tmp_path, "run-b", "rank_9.pdb", "flexref_9.pdb")

    assert first.canonical_script == second.canonical_script
    assert first.checksums == second.checksums


def test_canonical_mapping_keeps_model_content_in_identity(tmp_path):
    first = _named_model_mapping(tmp_path, "run-a", "rank_1.pdb", "flexref_1.pdb")
    second = _named_model_mapping(
        tmp_path,
        "run-b",
        "rank_9.pdb",
        "flexref_9.pdb",
        model_content="ATOM changed\n",
    )

    assert first.canonical_script == second.canonical_script
    assert first.checksums["canonical-input-1.pdb"] != (
        second.checksums["canonical-input-1.pdb"]
    )


def test_canonical_mapping_erases_count_but_keeps_seed(tmp_path):
    first = _indexed_model_mapping(tmp_path, "run-a", index=1, seed=1001)
    second = _indexed_model_mapping(tmp_path, "run-b", index=9, seed=1001)
    different_seed = _indexed_model_mapping(tmp_path, "run-c", index=1, seed=1002)

    assert first.canonical_script == second.canonical_script
    assert "canonical-count" in first.canonical_script
    assert "1001" in first.canonical_script
    assert first.canonical_script != different_seed.canonical_script


def test_canonical_mapping_resolves_count_suffix_restraint_file(tmp_path):
    mapping = _counted_restraint_mapping(tmp_path, count=3, suffixed_exists=True)

    assert mapping.dependency_paths[
        (tmp_path / "run" / "03_flexref" / "ambig.tbl_3").resolve()
    ] == "canonical-ambig.tbl"


def test_canonical_mapping_erases_count_suffix_base_alias(tmp_path):
    first = _renamed_counted_restraint_mapping(
        tmp_path,
        "run-a",
        base_name="z_ambig.tbl",
        count=3,
    )
    second = _renamed_counted_restraint_mapping(
        tmp_path,
        "run-b",
        base_name="a_ambig.tbl",
        count=9,
    )

    assert first.canonical_script == second.canonical_script
    assert first.checksums == second.checksums
    assert "z_ambig.tbl" not in first.canonical_script
    assert "a_ambig.tbl" not in second.canonical_script


def test_canonical_mapping_falls_back_to_base_restraint_file(tmp_path):
    mapping = _counted_restraint_mapping(tmp_path, count=3, suffixed_exists=False)

    assert mapping.dependency_paths[
        (tmp_path / "run" / "03_flexref" / "ambig.tbl").resolve()
    ] == "canonical-ambig.tbl"


def test_canonical_mapping_replaces_relative_sibling_dependency_path(tmp_path):
    work_dir = tmp_path / "run" / "01_rigidbody"
    input_pdb = tmp_path / "run" / "data" / "00_topoaa" / "structure_1.pdb"
    module, toppar, cns = _install(tmp_path, "install")
    work_dir.mkdir(parents=True)
    input_pdb.parent.mkdir(parents=True)
    input_pdb.write_text("ATOM\n", encoding="utf-8")

    mapping = build_canonical_mapping(
        'evaluate ($input_pdb = "../data/00_topoaa/structure_1.pdb")\n'
        "coor @@$input_pdb\n",
        envvars={"MODULE": str(module), "TOPPAR": str(toppar)},
        cns_exec=cns,
        output_files=[work_dir / "result.pdb"],
        work_dir=work_dir,
    )

    assert "canonical-input-1.pdb" in mapping.canonical_script
    assert "00_topoaa" not in mapping.canonical_script


def test_canonical_mapping_accepts_dependency_named_after_its_role(tmp_path):
    work_dir = tmp_path / "run" / "01_rigidbody"
    module, toppar, cns = _install(tmp_path, "install")
    restraints = tmp_path / "run" / "data" / "ambig.tbl"
    work_dir.mkdir(parents=True)
    restraints.parent.mkdir(parents=True, exist_ok=True)
    restraints.write_text("assign\n", encoding="utf-8")
    script = work_dir / "job.inp"
    script.write_text(
        f'evaluate ($ambig_fname = "{restraints}")\n'
        "noe @@$ambig_fname end\n"
        "inline @@MODULE:protocol.cns\n"
        "inline @@TOPPAR:protein.top\n",
        encoding="utf-8",
    )

    mapping = build_canonical_mapping(
        script,
        envvars={"MODULE": str(module), "TOPPAR": str(toppar)},
        cns_exec=cns,
        output_files=[work_dir / "result.pdb"],
        work_dir=work_dir,
    )

    assert "canonical-ambig.tbl" in mapping.canonical_script
    assert str(restraints) not in mapping.canonical_script


def test_canonical_dependency_names_follow_first_reference_order(tmp_path):
    first = _two_input_mapping(
        tmp_path,
        "run-a",
        first_name="z_model.pdb",
        second_name="a_model.pdb",
    )
    second = _two_input_mapping(
        tmp_path,
        "run-b",
        first_name="a_model.pdb",
        second_name="z_model.pdb",
    )

    assert first.canonical_script == second.canonical_script
    assert first.checksums["canonical-input-1.pdb"] == (
        second.checksums["canonical-input-1.pdb"]
    )
    assert first.checksums["canonical-input-2.pdb"] == (
        second.checksums["canonical-input-2.pdb"]
    )


def test_canonical_mapping_resolves_module_absolute_suffix_and_toppar_slash(tmp_path):
    work_dir = tmp_path / "run" / "03_emref"
    module, toppar, cns = _install(tmp_path, "install")
    work_dir.mkdir(parents=True)
    (module / "protein-ss-restraints-all.cns").write_text(
        "{ ss restraints }\n",
        encoding="utf-8",
    )
    (toppar / "dmso.pdb").write_text("DMSO\n", encoding="utf-8")

    mapping = build_canonical_mapping(
        "inline @@MODULE:/protein-ss-restraints-all.cns\n"
        "coor @@TOPPAR/dmso.pdb\n",
        envvars={"MODULE": str(module), "TOPPAR": str(toppar)},
        cns_exec=cns,
        output_files=[work_dir / "result.pdb"],
        work_dir=work_dir,
    )

    assert mapping.dependency_paths[
        (module / "protein-ss-restraints-all.cns").resolve()
    ] == "module/protein-ss-restraints-all.cns"
    assert mapping.dependency_paths[
        (toppar / "dmso.pdb").resolve()
    ] == "toppar/dmso.pdb"


def test_canonical_mapping_rejects_unresolved_reads(tmp_path):
    work_dir = tmp_path / "run" / "01_rigidbody"
    module, toppar, cns = _install(tmp_path, "install")
    work_dir.mkdir(parents=True)

    with pytest.raises(ValueError, match="unresolved reads"):
        build_canonical_mapping(
            "coor @@$missing\n",
            envvars={"MODULE": str(module), "TOPPAR": str(toppar)},
            cns_exec=cns,
            output_files=[work_dir / "result.pdb"],
            work_dir=work_dir,
        )


def test_cnsjob_exposes_canonical_mapping(monkeypatch, tmp_path):
    mapping, step = _mapping(tmp_path, "run", "1_rigidbody", "install")
    monkeypatch.chdir(step)
    job = CNSJob(
        step / "job.inp",
        envvars={
            "MODULE": str(tmp_path / "install" / "module"),
            "TOPPAR": str(tmp_path / "install" / "toppar"),
        },
        cns_exec=mapping.cns_exec,
        output_files=[step / "result.pdb"],
    )

    assert job.canonical_mapping().checksums == mapping.checksums


def test_compression_transparent_checksum(tmp_path):
    plain = tmp_path / "input.pdb"
    compressed = tmp_path / "input.pdb.gz"
    plain.write_bytes(b"ATOM\n")
    with gzip.open(compressed, "wb") as handle:
        handle.write(plain.read_bytes())
    assert compression_transparent_checksum(plain) == (
        compression_transparent_checksum(compressed)
    )


def test_canonical_mapping_keeps_install_reference_spelling(tmp_path):
    """``MODULE:``/``TOPPAR:`` references resolve through the environment already."""
    mapping, _ = _mapping(tmp_path, "run", "1_rigidbody", "install")

    assert "@@MODULE:protocol.cns" in mapping.canonical_script
    assert "@@TOPPAR:protein.top" in mapping.canonical_script
    assert "MODULE:module/" not in mapping.canonical_script
    assert "TOPPAR:toppar/" not in mapping.canonical_script
    # the pin names are unchanged; only the script text keeps its own spelling
    assert "module/protocol.cns" in mapping.checksums
    assert "toppar/protein.top" in mapping.checksums


def test_canonical_mapping_rewrites_unread_absolute_install_paths(tmp_path):
    """An install path reaches the key even when the recipe never reads it."""
    root = tmp_path / "run" / "1_rigidbody"
    module, toppar, cns = _install(tmp_path, "install")
    root.mkdir(parents=True)
    (toppar / "boxtyp20.pdb").write_text("BOX\n", encoding="utf-8")

    mapping = build_canonical_mapping(
        f'evaluate ($boxtyp20 = "{toppar / "boxtyp20.pdb"}")\n'
        "inline @@MODULE:protocol.cns\n",
        envvars={"MODULE": str(module), "TOPPAR": str(toppar)},
        cns_exec=cns,
        output_files=[root / "result.pdb"],
        work_dir=root,
    )

    assert 'evaluate ($boxtyp20 = "TOPPAR:boxtyp20.pdb")' in mapping.canonical_script
    assert str(toppar) not in mapping.canonical_script


def test_canonical_mapping_preserves_cns_variable_names(tmp_path):
    """Rewriting a dependency must not rename the symbol that refers to it."""
    root = tmp_path / "run" / "0_topoaa"
    module, toppar, cns = _install(tmp_path, "install")
    root.mkdir(parents=True)
    (root / "molA.pdb").write_text("ATOM\n", encoding="utf-8")

    mapping = build_canonical_mapping(
        'eval ($file = "molA.pdb")\n'
        "evaluate ($coor_infile = $file)\n"
        "coor @@$coor_infile\n",
        envvars={"MODULE": str(module), "TOPPAR": str(toppar)},
        cns_exec=cns,
        output_files=[root / "result.pdb"],
        work_dir=root,
    )

    assert 'eval ($file = "canonical-input-1.pdb")' in mapping.canonical_script
    assert "evaluate ($coor_infile = $file)" in mapping.canonical_script
    assert "eval (canonical-input-1.pdb" not in mapping.canonical_script


def test_canonical_mapping_output_survives_input_basename_collision(tmp_path):
    """An input and an output may share a basename in different directories."""
    root = tmp_path / "run" / "1_topocg"
    upstream = tmp_path / "run" / "0_topoaa"
    module, toppar, cns = _install(tmp_path, "install")
    root.mkdir(parents=True)
    upstream.mkdir(parents=True)
    (upstream / "shape_haddock.pdb").write_text("ATOM\n", encoding="utf-8")

    mapping = build_canonical_mapping(
        'evaluate ($input_pdb = "../0_topoaa/shape_haddock.pdb")\n'
        'evaluate ($output_pdb_filename = "shape_haddock.pdb")\n'
        "coor @@$input_pdb\n",
        envvars={"MODULE": str(module), "TOPPAR": str(toppar)},
        cns_exec=cns,
        output_files=[root / "shape_haddock.pdb"],
        work_dir=root,
    )

    assert (
        'evaluate ($output_pdb_filename = "canonical-output.pdb")'
        in mapping.canonical_script
    )
    assert "canonical-input-1.pdb" in mapping.canonical_script


def test_canonical_mapping_rejects_output_bound_to_an_input_pin(tmp_path):
    """The guard refuses a script that writes to something other than its output."""
    root = tmp_path / "run" / "1_rigidbody"
    module, toppar, cns = _install(tmp_path, "install")
    root.mkdir(parents=True)
    (root / "model.pdb").write_text("ATOM\n", encoding="utf-8")

    with pytest.raises(ValueError, match="output_pdb_filename"):
        build_canonical_mapping(
            'evaluate ($input_pdb = "model.pdb")\n'
            'evaluate ($output_pdb_filename = "elsewhere.pdb")\n'
            "coor @@$input_pdb\n",
            envvars={"MODULE": str(module), "TOPPAR": str(toppar)},
            cns_exec=cns,
            output_files=[root / "result.pdb"],
            work_dir=root,
        )


def test_canonical_mapping_normalizes_logging_only_parameter(tmp_path):
    root = tmp_path / "run" / "1_rigidbody"
    module, toppar, cns = _install(tmp_path, "install")
    root.mkdir(parents=True)
    (root / "input.pdb").write_text("ATOM\n", encoding="utf-8")

    mapping = build_canonical_mapping(
        'eval ($log_level="verbose")\n'
        'eval ($input_pdb="input.pdb")\n'
        'eval ($output_pdb_filename="result.pdb")\n'
        "coor @@$input_pdb\n",
        envvars={"MODULE": str(module), "TOPPAR": str(toppar)},
        cns_exec=cns,
        output_files=[root / "result.pdb"],
        work_dir=root,
    )

    assert 'eval ($log_level="canonical-log-level")' in mapping.canonical_script


def test_canonical_mapping_resolves_cgtoaa_indexed_variable_references(tmp_path):
    root = tmp_path / "run" / "2_cgtoaa"
    module, toppar, cns = _install(tmp_path, "install")
    root.mkdir(parents=True)
    for index in (1, 2):
        (root / f"input_{index}.psf").write_text("PSF\n", encoding="utf-8")

    mapping = build_canonical_mapping(
        "eval ($input_aa_psf_filename_1=\"input_1.psf\")\n"
        "eval ($input_aa_psf_filename_2=\"input_2.psf\")\n"
        "while ($nchain < 2) loop nloop1\n"
        "  structure @@$input_aa_psf_filename_$nchain end\n"
        "end loop nloop1\n",
        envvars={"MODULE": str(module), "TOPPAR": str(toppar)},
        cns_exec=cns,
        output_files=[root / "result.pdb"],
        work_dir=root,
    )

    assert [dependency.original_path.name for dependency in mapping.dependencies] == [
        "input_1.psf",
        "input_2.psf",
    ]


def test_logging_and_count_canonicalizations_are_limited_to_known_cns_uses():
    module_root = Path(__file__).parents[1] / "src" / "haddock" / "modules"
    recipes = module_root.glob("*/*/cns/*.cns")
    log_level_lines = []
    count_lines = []
    for recipe in recipes:
        for line in recipe.read_text(encoding="utf-8").splitlines():
            if "$log_level" in line:
                log_level_lines.append(line.strip())
            if re.search(r"\$count(?![A-Za-z0-9_])", line):
                count_lines.append(line.strip())

    assert log_level_lines
    assert all(
        line.startswith(("if ( $log_level", "elseif ( $log_level"))
        for line in log_level_lines
    )
    assert count_lines
    assert all(
        "display STRUCTURE NUMBER $count" in line
        or "$ambig_fname + \"_\" + encode($count)" in line
        for line in count_lines
    )


def _mapping(tmp_path: Path, run_name: str, step_name: str, install_name: str):
    root = tmp_path / run_name / step_name
    module, toppar, cns = _install(tmp_path, install_name)
    root.mkdir(parents=True)
    (root / "renamed-model.pdb").write_text("ATOM\n", encoding="utf-8")
    script = root / "job.inp"
    script.write_text(
        'evaluate ($input_pdb = "renamed-model.pdb")\n'
        "coor @@$input_pdb\n"
        "inline @@MODULE:protocol.cns\n"
        "inline @@TOPPAR:protein.top\n",
        encoding="utf-8",
    )
    return (
        build_canonical_mapping(
            script,
            envvars={"MODULE": str(module), "TOPPAR": str(toppar)},
            cns_exec=cns,
            output_files=[root / "result.pdb"],
            work_dir=root,
        ),
        root,
    )


def _named_model_mapping(
    tmp_path: Path,
    run_name: str,
    input_name: str,
    output_name: str,
    model_content: str = "ATOM\n",
):
    root = tmp_path / run_name / "03_flexref"
    module, toppar, cns = _install(tmp_path, f"{run_name}-install")
    root.mkdir(parents=True)
    (root / input_name).write_text(model_content, encoding="utf-8")
    return build_canonical_mapping(
        f'evaluate ($input_pdb = "{input_name}")\n'
        f'evaluate ($output_pdb_filename = "{output_name}")\n'
        "coor @@$input_pdb\n",
        envvars={"MODULE": str(module), "TOPPAR": str(toppar)},
        cns_exec=cns,
        output_files=[root / output_name],
        work_dir=root,
    )


def _indexed_model_mapping(
    tmp_path: Path,
    run_name: str,
    index: int,
    seed: int,
):
    root = tmp_path / run_name / "03_flexref"
    module, toppar, cns = _install(tmp_path, f"{run_name}-install")
    root.mkdir(parents=True)
    input_name = f"rank_{index}.pdb"
    output_name = f"flexref_{index}.pdb"
    (root / input_name).write_text("ATOM\n", encoding="utf-8")
    return build_canonical_mapping(
        f'evaluate ($input_pdb = "{input_name}")\n'
        f'evaluate ($output_pdb_filename = "{output_name}")\n'
        f"evaluate ($count = {index})\n"
        f"evaluate ($seed = {seed})\n"
        "coor @@$input_pdb\n",
        envvars={"MODULE": str(module), "TOPPAR": str(toppar)},
        cns_exec=cns,
        output_files=[root / output_name],
        work_dir=root,
    )


def _counted_restraint_mapping(
    tmp_path: Path,
    count: int,
    suffixed_exists: bool,
):
    root = tmp_path / "run" / "03_flexref"
    module, toppar, cns = _install(tmp_path, "install")
    root.mkdir(parents=True)
    (root / "model.pdb").write_text("ATOM\n", encoding="utf-8")
    (root / "ambig.tbl").write_text("assign base\n", encoding="utf-8")
    if suffixed_exists:
        (root / f"ambig.tbl_{count}").write_text("assign counted\n", encoding="utf-8")
    return build_canonical_mapping(
        'evaluate ($input_pdb = "model.pdb")\n'
        'evaluate ($ambig_fname = "ambig.tbl")\n'
        f"evaluate ($count = {count})\n"
        'evaluate ($filenam0 = $ambig_fname + "_" + encode($count))\n'
        "noe class ambi @@$filenam0 end\n"
        "coor @@$input_pdb\n",
        envvars={"MODULE": str(module), "TOPPAR": str(toppar)},
        cns_exec=cns,
        output_files=[root / "result.pdb"],
        work_dir=root,
    )


def _renamed_counted_restraint_mapping(
    tmp_path: Path,
    run_name: str,
    base_name: str,
    count: int,
):
    root = tmp_path / run_name / "03_flexref"
    module, toppar, cns = _install(tmp_path, f"{run_name}-install")
    root.mkdir(parents=True)
    (root / "model.pdb").write_text("ATOM\n", encoding="utf-8")
    (root / f"{base_name}_{count}").write_text("assign same\n", encoding="utf-8")
    return build_canonical_mapping(
        'evaluate ($input_pdb = "model.pdb")\n'
        f'evaluate ($ambig_fname = "{base_name}")\n'
        f"evaluate ($count = {count})\n"
        'evaluate ($filenam0 = $ambig_fname + "_" + encode($count))\n'
        "noe class ambi @@$filenam0 end\n"
        "coor @@$input_pdb\n",
        envvars={"MODULE": str(module), "TOPPAR": str(toppar)},
        cns_exec=cns,
        output_files=[root / "result.pdb"],
        work_dir=root,
    )


def _two_input_mapping(
    tmp_path: Path,
    run_name: str,
    first_name: str,
    second_name: str,
):
    root = tmp_path / run_name / "03_flexref"
    module, toppar, cns = _install(tmp_path, f"{run_name}-install")
    root.mkdir(parents=True)
    (root / first_name).write_text("ATOM first\n", encoding="utf-8")
    (root / second_name).write_text("ATOM second\n", encoding="utf-8")
    return build_canonical_mapping(
        f'evaluate ($input_pdb_1 = "{first_name}")\n'
        f'evaluate ($input_pdb_2 = "{second_name}")\n'
        "coor @@$input_pdb_1\n"
        "coor @@$input_pdb_2\n",
        envvars={"MODULE": str(module), "TOPPAR": str(toppar)},
        cns_exec=cns,
        output_files=[root / "result.pdb"],
        work_dir=root,
    )


def _install(tmp_path: Path, install_name: str) -> tuple[Path, Path, Path]:
    module = tmp_path / install_name / "module"
    toppar = tmp_path / install_name / "toppar"
    module.mkdir(parents=True, exist_ok=True)
    toppar.mkdir(parents=True, exist_ok=True)
    (module / "protocol.cns").write_text("{ module }\n", encoding="utf-8")
    (toppar / "protein.top").write_text("{ toppar }\n", encoding="utf-8")
    cns = tmp_path / install_name / "cns"
    # Two installations never hold the same executable bytes.  Writing the
    # install name into it keeps that true here, so that every test comparing
    # two installations is comparing what a real pair of them would look like.
    cns.write_text(f"#!/bin/sh\n# {install_name}\n", encoding="utf-8")
    cns.chmod(0o755)
    return module, toppar, cns
