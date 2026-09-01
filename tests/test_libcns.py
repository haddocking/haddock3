"""Test libcns."""

import os
import random
import tempfile
from pathlib import Path

import pytest

import haddock.modules
from haddock import EmptyPath
from haddock.libs import libcns
from haddock.libs.libcns import (
    _add_cg_backmapping_arguments,
    prepare_cns_input,
    prepare_expected_pdb,
    prepare_multiple_input,
)
from haddock.libs.libontology import Format, PDBFile, Persistent


@pytest.mark.parametrize(
    "value",
    [
        {"dic": 1},
        tuple([1, 2]),
        [1, 2],
        set([1, 2]),
    ],
)
def test_empty_vars_error(value):
    """Test empty vars of types that are not supported."""
    with pytest.raises(TypeError):
        libcns.filter_empty_vars(value)


@pytest.mark.parametrize(
    "value",
    [
        True,
        False,
        1,
        12,
        3453.543,
        "str",
        Path("path"),
        EmptyPath(),  # empty paths need to be written as ""
    ],
)
def test_empty_vars_True(value):
    """Test empty vars of types that are not supported."""
    result = libcns.filter_empty_vars(value)
    assert result is True
    assert type(result) is bool


@pytest.mark.parametrize(
    "value",
    [
        None,
        "",
        float("nan"),
    ],
)
def test_empty_vars_False(value):
    """Test empty vars of types that are not supported."""
    result = libcns.filter_empty_vars(value)
    assert result is False
    assert type(result) is bool


@pytest.fixture
def pdbfile():
    """Create a temporary file."""
    with tempfile.NamedTemporaryFile() as pdb_f, tempfile.NamedTemporaryFile() as top_f:
        pdb = PDBFile(
            file_name=pdb_f.name,
            path=Path(pdb_f.name).parent,
            topology=Persistent(
                file_name=top_f.name,
                path=Path(top_f.name).parent,
                file_type=Format.TOPOLOGY,
            ),
        )
        pdb.seed = 42  # type: ignore
        yield pdb


def test_load_workflow_params():
    """Test workflow params."""
    params = {
        "var1": 1,
        "var2": "some string",
        "var3": True,
        "var4": Path("some/path"),
        "var5": 5.5,
        "var6": "",
        "var7": None,
    }

    result = libcns.load_workflow_params(**params)

    expected = (
        f"{os.linesep}"
        f"! Parameters{os.linesep}"
        f"eval ($var1=1){os.linesep}"
        f'eval ($var2="some string"){os.linesep}'
        f"eval ($var3=true){os.linesep}"
        f'eval ($var4="some/path"){os.linesep}'
        f"eval ($var5=5.5){os.linesep}"
    )

    assert result == expected


def test_prepare_cns_input(pdbfile):

    # TODO: Improve this test to increase coverage AND refactor `prepare_cns_input`

    model_number = 1
    observed_cns_input = prepare_cns_input(
        model_number=model_number,
        input_element=pdbfile,
        step_path=Path("."),
        recipe_str="",
        defaults={},
        identifier="",
        ambig_fname="",
        native_segid=False,
        default_params_path=None,
        debug=False,
        seed=pdbfile.seed,
    )

    expected_cns_input = f"""
! Parameters
eval ($ambig_fname="")

! Input structure
structure
  @@{pdbfile.topology.rel_path}
end
coor @@{pdbfile.rel_path}
eval ($input_pdb_filename_1="{pdbfile.rel_path}")
eval ($ncomponents=0)
eval ($seed={pdbfile.seed})

! Output structure
eval ($output_pdb_filename="_{model_number}.pdb")
eval ($count=1)
"""

    assert observed_cns_input == expected_cns_input


def test_prepare_cns_input_omits_seed_when_one_is_not_declared(pdbfile):
    observed_cns_input = prepare_cns_input(
        model_number=1,
        input_element=pdbfile,
        step_path=Path("."),
        recipe_str="",
        defaults={},
        identifier="model",
        seed=None,
    )

    assert "$seed" not in observed_cns_input


def test_prepare_cns_input_replaces_archive_restraint_assignment(pdbfile):
    observed_cns_input = prepare_cns_input(
        model_number=1,
        input_element=pdbfile,
        step_path=Path("."),
        recipe_str="",
        defaults={"ambig_fname": Path("ambig.tbl.tgz")},
        identifier="model",
        ambig_fname=Path("ambig_1.tbl"),
    )

    assert "ambig.tbl.tgz" not in observed_cns_input
    assert observed_cns_input.count("$ambig_fname") == 1
    assert 'eval ($ambig_fname="ambig_1.tbl")' in observed_cns_input


def test_prepare_cns_input_expands_molecule_defaults_for_direct_call(pdbfile, mocker):
    """Direct callers receive one molecule-family value per component."""
    mocker.patch(
        "haddock.libs.libpdb.identify_chainseg",
        return_value=(["A"], []),
    )

    observed = prepare_cns_input(
        model_number=1,
        input_element=[pdbfile, pdbfile],
        step_path=Path("."),
        recipe_str="",
        defaults={"mol_shape_1": False},
        identifier="model",
    )

    assert "eval ($mol_shape_1=false)" in observed
    assert "eval ($mol_shape_2=false)" in observed


def test_cg_backmapping_keeps_psf_tbl_pairs_aligned_with_shapes(tmp_path):
    """A leading shape component must not shift CG-to-AA file pairings."""
    model = PDBFile(file_name="complex.pdb", path=tmp_path)
    model.aa_topology = [
        Persistent("shape.psf", Format.TOPOLOGY, tmp_path),
        Persistent("protein.psf", Format.TOPOLOGY, tmp_path),
    ]
    model.cgtoaa_tbl = [None, Path("protein.tbl")]
    model.shape = [True, False]

    observed = _add_cg_backmapping_arguments(model)

    assert "shape.psf" not in observed
    assert f'input_aa_psf_filename_1="{model.aa_topology[1].rel_path}"' in observed
    assert 'input_cgtbl_filename_1="protein.tbl"' in observed


def test_prepare_multiple_input(mocker):

    mocker.patch("haddock.libs.libpdb.identify_chainseg", return_value="A")

    observed_input = prepare_multiple_input(
        pdb_input_list=["pdb_1", "pdb_2"],
        psf_input_list=["psf_1"],
    )

    expected_input = """
! Input structure
structure
  @@psf_1
end
coor @@pdb_1
eval ($input_pdb_filename_1="pdb_1")
coor @@pdb_2
eval ($input_pdb_filename_2="pdb_2")
eval ($ncomponents=1)
"""

    assert observed_input == expected_input


def test_prepare_expected_pdb(pdbfile):

    model_number = random.randint(1, 100)
    output_pdb = prepare_expected_pdb(
        model_obj=pdbfile,
        model_nb=model_number,
        path=Path.cwd(),
        identifier="",
    )

    assert isinstance(output_pdb, PDBFile)
    assert output_pdb.file_name == f"_{model_number}.pdb"
    assert output_pdb.topology is not None
    assert output_pdb.topology.file_name == pdbfile.topology.file_name


def test_prepare_cns_input_with_ligand_files():
    """Test CNS input generation with per-model ligand files."""
    with tempfile.NamedTemporaryFile() as pdb_f, tempfile.NamedTemporaryFile() as top_f:
        # Create PDBFile with ligand files
        pdb = PDBFile(
            file_name=pdb_f.name,
            path=Path(pdb_f.name).parent,
            topology=Persistent(
                file_name=top_f.name,
                path=Path(top_f.name).parent,
                file_type=Format.TOPOLOGY,
            ),
            ligand_top_fname="ligand.top",
            ligand_param_fname="ligand.param",
        )
        pdb.seed = 42

        model_number = 1
        cns_input = prepare_cns_input(
            model_number=model_number,
            input_element=pdb,
            step_path=Path("."),
            recipe_str="",
            defaults={},
            identifier="",
            ambig_fname="",
            native_segid=False,
            default_params_path=None,
            debug=False,
            seed=pdb.seed,
        )

        # The per-model ligand *parameter* file is included: every recipe
        # reached through `prepare_cns_input` reads it.
        assert "ligand.param" in cns_input
        assert 'eval ($ligand_param_fname="ligand.param")' in cns_input

        # The ligand *topology* file is not.  Only the two topology recipes
        # read it, and neither of them comes through here, so emitting it
        # would leave a step-folder path in the generated input that CNS never
        # opens.  The topology it describes is already in the PSF by now.
        assert "ligand.top" not in cns_input
        assert "$ligand_top_fname" not in cns_input


def _seeded_model(directory: Path, name: str, content: str) -> PDBFile:
    """A model whose bytes are on disk, so its seed can be derived."""
    path = Path(directory, name)
    path.write_text(content, encoding="utf-8")
    return PDBFile(file_name=name, path=directory)


def test_derive_seed_ignores_the_job_position_and_name(tmp_path, monkeypatch):
    """The same structure yields the same seed wherever it is scheduled."""
    monkeypatch.chdir(tmp_path)
    first = _seeded_model(tmp_path, "rigidbody_1.pdb", "ATOM  one\n")
    renamed = _seeded_model(tmp_path, "rank_37.pdb", "ATOM  one\n")

    assert libcns.derive_seed(917, first) == libcns.derive_seed(917, renamed)


def test_derive_seed_follows_the_content_it_reads(tmp_path, monkeypatch):
    """A changed structure is a different computation and gets a new seed."""
    monkeypatch.chdir(tmp_path)
    original = _seeded_model(tmp_path, "model.pdb", "ATOM  one\n")
    changed = _seeded_model(tmp_path, "other.pdb", "ATOM  two\n")

    assert libcns.derive_seed(917, original) != libcns.derive_seed(917, changed)


def test_derive_seed_separates_repeats_of_one_job(tmp_path, monkeypatch):
    """Repeats of one job are additional sampling, not duplicates of it."""
    monkeypatch.chdir(tmp_path)
    model = _seeded_model(tmp_path, "model.pdb", "ATOM  one\n")

    seeds = [libcns.derive_seed(917, model, repeat) for repeat in range(8)]

    assert len(set(seeds)) == len(seeds)


def test_derive_seed_still_follows_iniseed(tmp_path, monkeypatch):
    """``iniseed`` keeps its meaning: changing it changes every seed."""
    monkeypatch.chdir(tmp_path)
    model = _seeded_model(tmp_path, "model.pdb", "ATOM  one\n")

    assert libcns.derive_seed(917, model) != libcns.derive_seed(4242, model)


def test_derive_seed_binds_the_inputs_in_order(tmp_path, monkeypatch):
    """Two molecules in the other order is a different docking job."""
    monkeypatch.chdir(tmp_path)
    receptor = _seeded_model(tmp_path, "receptor.pdb", "ATOM  receptor\n")
    ligand = _seeded_model(tmp_path, "ligand.pdb", "ATOM  ligand\n")

    assert libcns.derive_seed(917, [receptor, ligand]) != libcns.derive_seed(
        917, [ligand, receptor]
    )


def test_derive_seed_stays_in_the_range_cns_represents_exactly(tmp_path, monkeypatch):
    """CNS holds numbers as doubles; a derived seed must be an exact integer."""
    monkeypatch.chdir(tmp_path)

    for index in range(200):
        model = _seeded_model(tmp_path, f"model_{index}.pdb", f"ATOM  {index}\n")
        seed = libcns.derive_seed(9999999999999999, model, index)
        assert 1 <= seed < libcns.SEED_CEILING
        assert float(seed) == seed


def test_derive_seed_reads_a_compressed_model(tmp_path, monkeypatch):
    """A cleaned step holds ``model.pdb.gz``; that is the same content."""
    import gzip

    monkeypatch.chdir(tmp_path)
    plain = _seeded_model(tmp_path, "model.pdb", "ATOM  one\n")
    expected = libcns.derive_seed(917, plain)

    Path(tmp_path, "model.pdb.gz").write_bytes(
        gzip.compress(b"ATOM  one\n", mtime=0)
    )
    Path(tmp_path, "model.pdb").unlink()

    assert libcns.derive_seed(917, plain) == expected


def test_refinement_schedule_emits_replicas_in_rounds():
    """Every model gets its first replica before any model gets its second."""
    schedule = libcns.refinement_schedule(3, 2)

    assert schedule == [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)]


def test_refinement_schedule_is_prefix_stable_in_sampling_factor():
    """Raising ``sampling_factor`` appends jobs and renumbers none of them."""
    for factor in range(1, 6):
        schedule = libcns.refinement_schedule(4, factor)
        assert len(schedule) == 4 * factor
        for smaller in range(1, factor + 1):
            assert schedule[: 4 * smaller] == libcns.refinement_schedule(4, smaller)


def test_no_prepare_cns_input_module_reads_the_ligand_topology_variable():
    """`prepare_cns_input` omits `$ligand_top_fname`; this pins why that is safe.

    The per-model assignment is a step-folder path.  Being unread it is not a
    declared dependency, so nothing rewrites it and it would leak a locator
    into the job's canonical form.  Only the two topology recipes read it, and
    neither reaches CNS through `prepare_cns_input`.  If that ever changes,
    this fails rather than the omission silently dropping a real input.
    """
    modules_dir = Path(haddock.modules.__file__).parent
    offenders = sorted(
        module_dir.name
        for module_dir in modules_dir.glob("*/*")
        if (module_dir / "cns").is_dir()
        and "prepare_cns_input"
        in (module_dir / "__init__.py").read_text(encoding="utf-8")
        and any(
            "$ligand_top_fname" in recipe.read_text(encoding="utf-8", errors="ignore")
            for recipe in (module_dir / "cns").rglob("*.cns")
        )
    )
    assert offenders == []
