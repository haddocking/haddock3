import pytest
import tempfile
from haddock.libs.libsubprocess import CNSJob
from pathlib import Path
from typing import Generator

from . import golden_data, CNS_EXEC, has_cns


@pytest.fixture
def cns_output_pdb_filename() -> Generator[str, None, None]:
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as output_f:
        yield output_f.name


@pytest.fixture
def cns_output_filename() -> Generator[str, None, None]:
    with tempfile.NamedTemporaryFile(suffix=".out", delete=False) as out_f:
        yield out_f.name


@pytest.fixture
def cnsjob(cns_input_filename, cns_output_filename):
    return CNSJob(
        input_file=Path(cns_input_filename),
        output_file=Path(cns_output_filename),
        cns_exec=CNS_EXEC,
    )


@pytest.fixture
def cnsjob_no_files(cns_inp_str):

    yield CNSJob(
        input_file=cns_inp_str,
        cns_exec=CNS_EXEC,
    )


@pytest.fixture
def cns_seed_filename(cns_output_filename) -> Generator[str, None, None]:
    seed_filename = Path(Path(cns_output_filename).stem).with_suffix(".seed")
    yield str(seed_filename)
    seed_filename.unlink(missing_ok=True)


@pytest.fixture
def cns_inp_str(cns_seed_filename, cns_output_pdb_filename):
    yield f"""
structure
    @@{golden_data}/prot.psf
end
coor @@{golden_data}/prot.pdb

write coordinates format=pdbo output={cns_output_pdb_filename} end

set display={cns_seed_filename} end
display evaluate($seed=42)
close {cns_seed_filename} end

stop"""


@pytest.fixture
def cns_input_filename(cns_inp_str) -> Generator[str, None, None]:
    with tempfile.NamedTemporaryFile(suffix=".inp", delete=False) as inp_f:
        inp_f.write(cns_inp_str.encode())
        inp_f.seek(0)

        yield inp_f.name


def test_cnsjob_run(cnsjob, cns_output_pdb_filename):

    cnsjob.run(
        compress_inp=False,
        compress_out=False,
        compress_seed=False,
    )

    assert Path(cns_output_pdb_filename).exists()
    assert Path(cns_output_pdb_filename).stat().st_size > 0


def test_cnsjob_run_compress_inp(cnsjob, cns_input_filename, cns_output_pdb_filename):
    cnsjob.run(
        compress_inp=True,
        compress_out=False,
        compress_seed=False,
    )

    assert Path(f"{cns_input_filename}.gz").exists()
    assert Path(f"{cns_input_filename}.gz").stat().st_size > 0

    assert Path(cns_output_pdb_filename).exists()
    assert Path(cns_output_pdb_filename).stat().st_size > 0


def test_cnsjob_run_compress_out(cnsjob, cns_output_filename, cns_output_pdb_filename):
    cnsjob.run(
        compress_inp=False,
        compress_out=True,
        compress_seed=False,
    )

    assert Path(f"{cns_output_filename}.gz").exists()
    assert Path(f"{cns_output_filename}.gz").stat().st_size > 0

    assert Path(cns_output_pdb_filename).exists()
    assert Path(cns_output_pdb_filename).stat().st_size > 0


def test_cnsjob_compress_seed(cnsjob, cns_output_pdb_filename, cns_seed_filename):
    cnsjob.run(
        compress_inp=False,
        compress_out=False,
        compress_seed=True,
    )

    assert Path(f"{cns_seed_filename}.gz").exists()
    assert Path(f"{cns_seed_filename}.gz").stat().st_size > 0
    Path(f"{cns_seed_filename}.gz").unlink()

    assert Path(cns_output_pdb_filename).exists()
    assert Path(cns_output_pdb_filename).stat().st_size > 0


def test_cnsjob_nofiles(cnsjob_no_files, cns_output_pdb_filename):

    cnsjob_no_files.run(
        compress_inp=False,
        compress_out=False,
        compress_seed=False,
    )

    assert Path(cns_output_pdb_filename).exists()
    assert Path(cns_output_pdb_filename).stat().st_size > 0
