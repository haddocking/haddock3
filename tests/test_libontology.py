"""Test functions and methods in haddock.libs.libontology."""
import json
import math
import tempfile
import shutil
from pathlib import Path

import pytest

from haddock.core.typing import Generator
from haddock.libs.libontology import (
    Format,
    ModuleIO,
    Molecule,
    PDBFile,
    Persistent,
    RMSDFile,
    TopologyFile,
)

from . import golden_data


@pytest.fixture
def input_pdbfile() -> Generator[PDBFile, None, None]:
    with tempfile.NamedTemporaryFile() as f:
        yield PDBFile(file_name=f.name)


@pytest.fixture
def output_pdbfile() -> Generator[PDBFile, None, None]:
    with tempfile.NamedTemporaryFile() as f:
        yield PDBFile(file_name=f.name)


@pytest.fixture
def moduleio_with_pdbfile_list(input_pdbfile, output_pdbfile):
    m = ModuleIO()
    m.input = [input_pdbfile]
    m.output = [output_pdbfile, output_pdbfile]
    return m


@pytest.fixture
def moduleio_with_pdbfile_dict(output_pdbfile):
    m = ModuleIO()
    m.input = []
    m.output = [
        {0: output_pdbfile, 1: output_pdbfile},
        {0: output_pdbfile, 1: output_pdbfile},
    ]
    return m


@pytest.fixture
def module_io_with_persistent():
    m = ModuleIO()

    # Create 10 random PDBFile instances and add them to the output list
    file_list = [tempfile.NamedTemporaryFile(delete=False) for _ in range(10)]

    for f in file_list:
        p = Persistent(file_name=f.name, file_type=Format.PDB, path=Path(f.name).parent)
        m.output.append(p)

    yield m

    for f in file_list:
        if Path(f.name).exists():
            Path(f.name).unlink()


@pytest.fixture
def io_data() -> dict:
    return {"input": ["input"], "output": ["output"]}


@pytest.fixture
def io_json_file(io_data) -> Generator[Path, None, None]:
    with tempfile.NamedTemporaryFile(mode="w+") as f:

        json.dump(io_data, f)

        f.seek(0)

        yield Path(f.name)


@pytest.fixture
def molecule():
    return Molecule(None)


@pytest.fixture
def protein():
    return Path(golden_data, "protein.pdb")


@pytest.fixture
def ensemble_header_w_md5():
    return Path(golden_data, "ens_header.pdb")


def test_get_md5(
        molecule: Molecule,
        ensemble_header_w_md5: Path,
        protein: Path,
        ):
    """Test get_md5 method."""
    observed_md5_dic = molecule.get_md5(ensemble_header_w_md5)
    expected_md5_dic = {
        1: '71098743056e0b95fbfafff690703761',
        2: 'f7ab0b7c751adf44de0f25f53cfee50b',
        3: '41e028d8d28b8d97148dc5e548672142',
        4: '761cb5da81d83971c2aae2f0b857ca1e',
        5: '6c438f941cec7c6dc092c8e48e5b1c10',
        }

    assert observed_md5_dic == expected_md5_dic
    observed_md5_dic = molecule.get_md5(protein)
    assert observed_md5_dic == {}


def test_get_ensemble_origin(
        molecule: Molecule,
        ensemble_header_w_md5: Path,
        protein: Path,
        ):
    """Test get_ensemble_origin method."""
    expected_origin_dic = {
        1: 'T161-hybrid-fit-C2-NCS_complex_100w',
        2: 'T161-hybrid-fit-C2-NCS_complex_101w',
        3: 'T161-hybrid-fit-C2-NCS_complex_102w',
        4: 'T161-hybrid-fit-C2-NCS_complex_103w',
        5: 'T161-hybrid-fit-C2-NCS_complex_104w',
        6: '73b07fb2ab6b3245_t264_1',
        }
    observed_origin = molecule.get_ensemble_origin(ensemble_header_w_md5)
    assert observed_origin == expected_origin_dic
    observed_origin = molecule.get_ensemble_origin(protein)
    assert observed_origin == {}


def test_load_single_pdb(molecule: Molecule, protein: Path):
    """Test casting into PDBFile."""
    with tempfile.TemporaryDirectory('.') as tempdir:
        tmp_protein = Path(tempdir, protein.name)
        shutil.copyfile(protein, tmp_protein)
        # Re-initialize with a actual protein
        molecule.__init__(tmp_protein)
        assert isinstance(molecule.pdb_files, dict)
        for pdbfile in molecule.pdb_files.values():
            assert isinstance(pdbfile, PDBFile)
        assert len(molecule) == 1


def test_persistent_init():
    target_path = Path(".")
    file_name = Path("/some/path/test")

    p = Persistent(file_name=file_name, file_type=Format.PDB, path=target_path)

    assert p.created is not None
    assert p.file_name == file_name.name
    assert p.file_type == Format.PDB
    assert p.path == str(target_path.resolve())
    assert p.full_name == str(target_path / file_name.name)
    assert p.rel_path == Path("..", Path(target_path).name, file_name)
    assert p.md5 is None
    assert p.restr_fname is None


def test_persistent_is_present():
    with tempfile.NamedTemporaryFile() as f:

        p = Persistent(file_name=f.name, file_type=Format.PDB, path=Path(f.name).parent)

        is_present = p.is_present()

        assert is_present

    p = Persistent(file_name="non_existent_file", file_type=Format.PDB)

    assert not p.is_present()


def test_pdbfile_init_empty():

    pdbfile = PDBFile(file_name="test.pdb")

    assert pdbfile.file_name == "test.pdb"
    assert pdbfile.file_type == Format.PDB
    assert pdbfile.path == str(Path.cwd())
    assert pdbfile.md5 is None
    assert pdbfile.restr_fname is None
    assert pdbfile.topology is None
    assert math.isnan(pdbfile.score)
    assert pdbfile.ori_name is None
    assert pdbfile.clt_id is None
    assert pdbfile.clt_rank is None
    assert pdbfile.clt_model_rank is None
    assert math.isnan(pdbfile.len)
    assert pdbfile.unw_energies is None
    assert pdbfile.seed is None


def test_pdbfile_init():
    target_path = Path(".")
    file_name = Path("/some/path/test")
    topology = "some topology"
    score = 10.0
    md5 = "some md5"
    restr_fname = "some restraint file"

    pdbfile = PDBFile(
        file_name=file_name,
        topology=topology,
        path=target_path,
        score=score,
        md5=md5,
        restr_fname=restr_fname,
    )

    assert pdbfile.file_name == file_name.name
    assert pdbfile.file_type == Format.PDB
    assert pdbfile.path == str(target_path.resolve())
    assert pdbfile.full_name == str(target_path / file_name.name)
    assert pdbfile.rel_path == Path("..", Path(target_path).name, file_name)
    assert pdbfile.md5 == md5
    assert pdbfile.restr_fname == restr_fname
    assert pdbfile.topology == topology
    assert pdbfile.score == score
    assert pdbfile.ori_name is None
    assert pdbfile.clt_id is None
    assert pdbfile.clt_rank is None
    assert pdbfile.clt_model_rank is None
    assert pdbfile.len == score
    assert pdbfile.unw_energies is None


def test_rmsdfile_init():

    npairs = 42
    file_name = Path("/some/path/test.rmsd")
    rmsdfile = RMSDFile(file_name=file_name, npairs=npairs)

    assert rmsdfile.file_name == file_name.name
    assert rmsdfile.file_type == Format.MATRIX
    assert rmsdfile.path == str(Path.cwd())
    assert rmsdfile.npairs == npairs


def test_topologyfile_init():

    file_name = Path("/some/path/test.top")
    topologyfile = TopologyFile(file_name=file_name)

    assert topologyfile.file_name == file_name.name
    assert topologyfile.file_type == Format.TOPOLOGY
    assert topologyfile.path == str(Path.cwd())


def test_moduleio_init():

    moduleio = ModuleIO()

    assert isinstance(moduleio.input, list)
    assert isinstance(moduleio.output, list)


def test_moduleio_add_anything():

    input = "literally anything"

    moduleio = ModuleIO()

    moduleio.add(input, mode="i")

    assert moduleio.input == [input]
    assert moduleio.output == []

    moduleio = ModuleIO()

    moduleio.add(input, mode="anything not i")

    assert moduleio.input == []
    assert moduleio.output == [input]


def test_moduleio_add_list():

    input = ["literally", "anything"]

    moduleio = ModuleIO()

    moduleio.add(input, mode="i")

    assert moduleio.input == ["literally", "anything"]
    assert moduleio.output == []

    moduleio = ModuleIO()

    moduleio.add(input, mode="anything not i")

    assert moduleio.input == []
    assert moduleio.output == ["literally", "anything"]


def test_moduleio_save(mocker, moduleio_with_pdbfile_list):

    with tempfile.NamedTemporaryFile() as temp_module_io_f:
        mocker.patch("haddock.core.defaults", temp_module_io_f.name)

        observed_io_filename = moduleio_with_pdbfile_list.save(
            filename=temp_module_io_f.name
        )

        assert observed_io_filename.exists()
        assert observed_io_filename.stat().st_size > 0

        expected_io_filename = Path.cwd() / temp_module_io_f.name

        assert observed_io_filename == expected_io_filename

        with open(observed_io_filename, "r") as f:
            observed_data = json.load(f)

        assert isinstance(observed_data, dict)


def test_moduleio_load(io_json_file, io_data):

    moduleio = ModuleIO()
    moduleio.load(filename=io_json_file)

    assert moduleio.input == io_data["input"]
    assert moduleio.output == io_data["output"]


def test_moduleio_retrieve_models_list(moduleio_with_pdbfile_list):

    result = moduleio_with_pdbfile_list.retrieve_models()

    assert isinstance(result, list)
    assert isinstance(result[0], PDBFile)
    assert isinstance(result[1], PDBFile)


def test_moduleio_retrieve_models_dict(moduleio_with_pdbfile_dict):

    result = moduleio_with_pdbfile_dict.retrieve_models(
        crossdock=True, individualize=True
    )

    assert isinstance(result, list)
    assert isinstance(result[0], PDBFile)

    result = moduleio_with_pdbfile_dict.retrieve_models(
        crossdock=True, individualize=False
    )

    assert isinstance(result, list)
    assert isinstance(result[0], tuple)
    assert isinstance(result[0][0], PDBFile)

    result = moduleio_with_pdbfile_dict.retrieve_models(
        crossdock=False, individualize=True
    )

    assert isinstance(result, list)
    assert isinstance(result[0], PDBFile)

    result = moduleio_with_pdbfile_dict.retrieve_models(
        crossdock=False, individualize=False
    )

    assert isinstance(result, list)
    assert isinstance(result[0], tuple)
    assert isinstance(result[0][0], PDBFile)


def test_moduleio_check_faulty(mocker, module_io_with_persistent):

    mocker.patch.object(module_io_with_persistent, "remove_missing", return_value=None)

    result = module_io_with_persistent.check_faulty()

    assert isinstance(result, float)
    assert result == pytest.approx(0.0)

    # Remove 10% of the files
    total = len(module_io_with_persistent.output)
    to_remove = int(total * 0.1)

    for i in range(to_remove):
        module_io_with_persistent.output[i].rel_path.unlink()

    result = module_io_with_persistent.check_faulty()

    assert result == pytest.approx(10.0)


def test_moduleio_remove_missing(module_io_with_persistent):

    # Remove the first file
    first_file = module_io_with_persistent.output[0].rel_path
    first_file.unlink()

    module_io_with_persistent.remove_missing()

    assert len(module_io_with_persistent.output) == 9

    # Make sure the first file is not in the list anymore
    assert first_file not in [p.rel_path for p in module_io_with_persistent.output]
