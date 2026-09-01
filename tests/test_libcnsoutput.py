import gzip
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from haddock.core.exceptions import CNSRunningError
from haddock.libs.libcnsoutput import (
    is_normalized_cns_pdb,
    normalize_cns_pdb,
)
from haddock.libs.libsubprocess import CNSJob


def test_normalize_cns_pdb_removes_run_volatile_remarks(tmp_path):
    pdb = tmp_path / "model.pdb"
    pdb.write_text(
        "\n".join(
            [
                "REMARK FILENAME= /tmp/run_1/model.pdb",
                "REMARK initial structure 1 - ../run_1/input.pdb",
                "REMARK DATE: 27-Jun-2026 12:34:56",
                "REMARK HADDOCK stats for rigidbody_1.pdb",
                "REMARK score: 1.23",
                "ATOM      1  CA  ALA A   1       0.000   0.000   0.000",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    changed = normalize_cns_pdb(pdb)

    assert changed is True
    assert pdb.read_text(encoding="utf-8") == (
        "REMARK score: 1.23\n"
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n"
        "END\n"
    )


def test_normalize_cns_pdb_leaves_stable_file_unchanged(tmp_path):
    pdb = tmp_path / "model.pdb"
    content = "REMARK score: 1.23\nEND\n"
    pdb.write_text(content, encoding="utf-8")

    changed = normalize_cns_pdb(pdb)

    assert changed is False
    assert pdb.read_text(encoding="utf-8") == content


def test_normalize_cns_pdb_handles_gzip_artifacts(tmp_path):
    pdb = tmp_path / "model.pdb.gz"
    pdb.write_bytes(gzip.compress(b"REMARK DATE: volatile\nATOM\n", mtime=0))

    assert normalize_cns_pdb(pdb) is True
    assert gzip.decompress(pdb.read_bytes()) == b"ATOM\n"


def test_normalize_cns_pdb_preserves_non_utf8_stable_bytes(tmp_path):
    pdb = tmp_path / "model.pdb"
    pdb.write_bytes(
        b"REMARK DATE: volatile\n"
        b"ATOM      1  X   BIN A   1       0.000   0.000   0.000\xff\n"
    )

    assert normalize_cns_pdb(pdb) is True
    assert pdb.read_bytes() == (
        b"ATOM      1  X   BIN A   1       0.000   0.000   0.000\xff\n"
    )


def test_normalize_cns_pdb_does_not_mutate_hardlink_source(tmp_path):
    source = tmp_path / "source.pdb"
    destination = tmp_path / "destination.pdb"
    source.write_text("REMARK DATE: volatile\nATOM\n", encoding="utf-8")
    destination.hardlink_to(source)

    assert normalize_cns_pdb(destination) is True
    assert source.read_text(encoding="utf-8") == "REMARK DATE: volatile\nATOM\n"
    assert destination.read_text(encoding="utf-8") == "ATOM\n"


def test_is_normalized_cns_pdb(tmp_path):
    pdb = tmp_path / "model.pdb"
    pdb.write_text("REMARK DATE: volatile\nATOM\n", encoding="utf-8")

    assert is_normalized_cns_pdb(pdb) is False
    normalize_cns_pdb(pdb)
    assert is_normalized_cns_pdb(pdb) is True


def test_cnsjob_run_normalizes_output_pdb(monkeypatch, tmp_path):
    output_pdb = tmp_path / "model.pdb"
    output_psf = tmp_path / "model.psf"
    output_pdb.write_text(
        "REMARK DATE: volatile\n"
        "REMARK initial structure 1 - ../run_1/input.pdb\n"
        "REMARK score: 1.23\n",
        encoding="utf-8",
    )
    output_psf.write_text(PSF_WITH_DATE, encoding="utf-8")
    job = CNSJob(
        input_file="cns input stream",
        cns_exec=_executable(tmp_path),
        output_files=[output_pdb, output_psf],
    )
    _mock_popen(monkeypatch)

    job.run()

    assert output_pdb.read_text(encoding="utf-8") == "REMARK score: 1.23\n"
    assert "created by user" not in output_psf.read_text(encoding="utf-8")


def test_cnsjob_normalizes_relative_outputs_from_creation_dir(tmp_path, monkeypatch):
    work_dir = tmp_path / "run" / "1_rigidbody"
    other_dir = tmp_path / "other"
    work_dir.mkdir(parents=True)
    other_dir.mkdir()
    output_pdb = work_dir / "model.pdb"
    output_pdb.write_text("REMARK DATE: volatile\nATOM\n", encoding="utf-8")

    monkeypatch.chdir(work_dir)
    job = CNSJob(
        input_file="cns input stream",
        cns_exec=_executable(tmp_path),
        output_files=["model.pdb"],
    )
    monkeypatch.chdir(other_dir)

    job.normalize_outputs()

    assert output_pdb.read_text(encoding="utf-8") == "ATOM\n"


def test_cnsjob_rejects_mismatched_declared_output(tmp_path):
    with pytest.raises(ValueError, match="does not match"):
        CNSJob(
            input_file='eval ($output_pdb_filename="actual.pdb")\n',
            cns_exec=_executable(tmp_path),
            output_files=[tmp_path / "declared.pdb"],
        )


def test_cnsjob_publishes_complete_normalized_staged_outputs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    public = tmp_path / "model.pdb"
    job = CNSJob(
        input_file='eval ($output_pdb_filename="model.pdb")\n',
        cns_exec=_executable(tmp_path),
        output_files=[public],
    )
    script, partial_outputs = job._partial_output_script()
    partial = partial_outputs[public]
    partial.write_text("REMARK DATE: volatile\nATOM\n", encoding="utf-8")

    assert "model.partial.pdb" in script
    job._publish_partial_outputs(partial_outputs)

    assert public.read_text(encoding="utf-8") == "ATOM\n"
    assert not partial.exists()


def test_cnsjob_retains_the_executed_input(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    input_file = tmp_path / "job.inp"
    input_file.write_text('eval ($output_pdb_filename="model.pdb")\n')
    job = CNSJob(
        input_file=input_file,
        cns_exec=_executable(tmp_path),
        output_files=[tmp_path / "model.pdb"],
    )

    script = job.prepare_execution_input()

    assert input_file.read_text() == script
    assert "model.partial.pdb" in script


def test_cnsjob_discards_incomplete_staged_outputs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    job = CNSJob(
        input_file='eval ($output_pdb_filename="model.pdb")\n',
        cns_exec=_executable(tmp_path),
        output_files=[tmp_path / "model.pdb"],
    )
    job.prepare_execution_input()
    partial = next(iter(job._pending_partial_outputs.values()))
    partial.write_text("incomplete\n")

    job._discard_partial_outputs()

    assert not partial.exists()


def test_cnsjob_ignores_only_the_benign_denormal_stderr(tmp_path, monkeypatch):
    job = CNSJob(input_file="stop\n", cns_exec=_executable(tmp_path))
    popen = _mock_popen(monkeypatch)
    popen.return_value.communicate.return_value = (
        b"normal output",
        b"Note: The following floating-point exceptions are signalling: "
        b"IEEE_DENORMAL\n",
    )
    popen.return_value.returncode = 0

    assert job.run() == b"normal output"


@pytest.mark.parametrize(
    ("stdout", "stderr", "returncode"),
    [
        (b"^^^^^\n", b"", 0),
        (b"normal output", b"real stderr\n", 0),
        (b"normal output", b"", 2),
    ],
)
def test_cnsjob_preserves_hard_failure_detection(
    tmp_path, monkeypatch, stdout, stderr, returncode
):
    job = CNSJob(input_file="stop\n", cns_exec=_executable(tmp_path))
    popen = _mock_popen(monkeypatch)
    popen.return_value.communicate.return_value = (stdout, stderr)
    popen.return_value.returncode = returncode

    with pytest.raises(CNSRunningError):
        job.run()


def _executable(tmp_path: Path) -> Path:
    path = tmp_path / "cns"
    path.write_text("#!/bin/sh\n", encoding="utf-8")
    path.chmod(0o755)
    return path


def _mock_popen(monkeypatch):
    popen = MagicMock()
    popen.return_value.communicate.return_value = (b"output", None)
    monkeypatch.setattr("subprocess.Popen", popen)
    return popen


PSF_WITH_DATE = """data_cns_mtf

_cns_mtf.title
; FILENAME="model_haddock.psf"
  disulphide added: from A    6    to A    127
  DATE:31-Aug-2026  01:17:08       created by user: unknown
  VERSION:1.3U
;

_cns_mtf.id   1
"""


def test_psf_normalization_removes_only_the_date_stamp(tmp_path):
    from haddock.libs.libcnsoutput import (
        is_normalized_cns_psf,
        normalize_cns_psf,
    )

    path = tmp_path / "model_haddock.psf"
    path.write_text(PSF_WITH_DATE, encoding="utf-8")
    assert not is_normalized_cns_psf(path)

    assert normalize_cns_psf(path) is True
    text = path.read_text(encoding="utf-8")

    assert "created by user" not in text
    assert "model_haddock.psf" not in text
    assert text.splitlines() == [
        "data_cns_mtf",
        "",
        "_cns_mtf.title",
        "; HADDOCK3 normalized topology",
        "  disulphide added: from A    6    to A    127",
        "  VERSION:1.3U",
        ";",
        "",
        "_cns_mtf.id   1",
    ]
    assert is_normalized_cns_psf(path)
    assert normalize_cns_psf(path) is False


def test_psf_normalization_makes_two_runs_agree(tmp_path):
    from haddock.libs.libcnsoutput import normalize_cns_psf_bytes

    monday = PSF_WITH_DATE.encode("utf-8")
    tuesday = PSF_WITH_DATE.replace("01:17:08", "09:42:55").encode("utf-8")

    assert monday != tuesday
    assert normalize_cns_psf_bytes(monday) == normalize_cns_psf_bytes(tuesday)


def test_psf_normalization_does_not_touch_structural_data(tmp_path):
    from haddock.libs.libcnsoutput import normalize_cns_psf_bytes

    body = "  DATE: 1 2 3 4\n  1 A    6    CYS  SG   SG     0.000000\n"
    assert normalize_cns_psf_bytes(body.encode()) == body.encode()
