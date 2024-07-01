import pytest
import tempfile
from haddock.libs.libsubprocess import CNSJob
from pathlib import Path
from typing import Generator

from . import golden_data, CNS_EXEC, has_cns


@pytest.fixture
def cns_output_pdb() -> Generator[str, None, None]:
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as output_f:
        yield output_f.name


@pytest.fixture
def cns_output() -> Generator[str, None, None]:
    with tempfile.NamedTemporaryFile(suffix=".out", delete=False) as out_f:
        yield out_f.name


@pytest.fixture
def cnsjob(cns_input, cns_output):
    return CNSJob(
        input_file=Path(cns_input),
        output_file=Path(cns_output),
        cns_exec=CNS_EXEC,
    )


@pytest.fixture
def cns_seed_file(cns_output) -> Generator[str, None, None]:
    yield str(Path(Path(cns_output).stem).with_suffix(".seed"))


@pytest.fixture
def cns_input(cns_output_pdb: str, cns_seed_file) -> Generator[str, None, None]:
    with tempfile.NamedTemporaryFile(suffix=".inp", delete=False) as inp_f:
        inp_str = f"""
structure
    @@{golden_data}/prot.psf
end
coor @@{golden_data}/prot.pdb

write coordinates format=pdbo output={cns_output_pdb} end

set display={cns_seed_file} end
display evaluate($seed=42)
close {cns_seed_file} end

stop"""
        inp_f.write(inp_str.encode())
        inp_f.seek(0)

        yield inp_f.name


def test_cnsjob_run(cnsjob, cns_output_pdb):

    cnsjob.run(
        compress_inp=False,
        compress_out=False,
        compress_seed=False,
    )

    assert Path(cns_output_pdb).exists()
    assert Path(cns_output_pdb).stat().st_size > 0


def test_cnsjob_run_compress_inp(cnsjob, cns_input, cns_output_pdb):
    cnsjob.run(
        compress_inp=True,
        compress_out=False,
        compress_seed=False,
    )

    assert Path(f"{cns_input}.gz").exists()
    assert Path(f"{cns_input}.gz").stat().st_size > 0

    assert Path(cns_output_pdb).exists()
    assert Path(cns_output_pdb).stat().st_size > 0


def test_cnsjob_run_compress_out(cnsjob, cns_output, cns_output_pdb):
    cnsjob.run(
        compress_inp=False,
        compress_out=True,
        compress_seed=False,
    )

    assert Path(f"{cns_output}.gz").exists()
    assert Path(f"{cns_output}.gz").stat().st_size > 0

    assert Path(cns_output_pdb).exists()
    assert Path(cns_output_pdb).stat().st_size > 0


def test_cnsjob_compress_seed(cnsjob, cns_output_pdb, cns_seed_file):
    cnsjob.run(
        compress_inp=False,
        compress_out=False,
        compress_seed=True,
    )

    assert Path(f"{cns_seed_file}.gz").exists()
    assert Path(f"{cns_seed_file}.gz").stat().st_size > 0

    assert Path(cns_output_pdb).exists()
    assert Path(cns_output_pdb).stat().st_size > 0
