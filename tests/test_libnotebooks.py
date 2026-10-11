import pytest
import gzip
from unittest.mock import patch, MagicMock
from haddock.libs import libnotebooks
from . import has_notebook


# Reference structure: two chains, non-collinear CA atoms so that the
# superposition is well defined
REFERENCE_PDB = """ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  CA  ALA A   2       3.800   0.000   0.000  1.00  0.00           C
ATOM      3  CA  ALA A   3       3.800   3.800   0.000  1.00  0.00           C
ATOM      4  CA  ALA B   1       0.000   0.000  10.000  1.00  0.00           C
ATOM      5  CA  ALA B   2       3.800   0.000  10.000  1.00  0.00           C
ATOM      6  CA  ALA B   3       3.800   3.800  10.000  1.00  0.00           C
END
"""

# Ensemble of two models, each one a rigid translation of the reference
# (+5 Å along x for model 1, +7 Å along y for model 2) so that the RMSD
# after superposition is zero
ENSEMBLE_PDB = """MODEL        1
ATOM      1  CA  ALA A   1       5.000   0.000   0.000  1.00  0.00           C
ATOM      2  CA  ALA A   2       8.800   0.000   0.000  1.00  0.00           C
ATOM      3  CA  ALA A   3       8.800   3.800   0.000  1.00  0.00           C
ATOM      4  CA  ALA B   1       5.000   0.000  10.000  1.00  0.00           C
ATOM      5  CA  ALA B   2       8.800   0.000  10.000  1.00  0.00           C
ATOM      6  CA  ALA B   3       8.800   3.800  10.000  1.00  0.00           C
ENDMDL
MODEL        2
ATOM      1  CA  ALA A   1       0.000   7.000   0.000  1.00  0.00           C
ATOM      2  CA  ALA A   2       3.800   7.000   0.000  1.00  0.00           C
ATOM      3  CA  ALA A   3       3.800  10.800   0.000  1.00  0.00           C
ATOM      4  CA  ALA B   1       0.000   7.000  10.000  1.00  0.00           C
ATOM      5  CA  ALA B   2       3.800   7.000  10.000  1.00  0.00           C
ATOM      6  CA  ALA B   3       3.800  10.800  10.000  1.00  0.00           C
ENDMDL
END
"""


# Sample minimal PDB string for tests
@pytest.fixture
def sample_pdb():
    pdb_string = """ATOM      1  CA  MET A   1      11.104  13.207   2.527  1.00  0.00           CA 
ATOM      2  CA  MET A   2      16.104  13.207   2.527  1.00  0.00           CA 
ATOM      3  CA  MET A   3      21.104  13.207   2.527  1.00  0.00           CA 
END
"""
    return pdb_string


@pytest.fixture
def sample_pdb_as_file(tmp_path, sample_pdb):
    fname = tmp_path / "test.pdb"
    with open(fname, "w") as f:
        f.write(sample_pdb)

    yield fname


@pytest.fixture
def sample_pdb_as_gz(tmp_path, sample_pdb):
    fname = tmp_path / "test.pdb.gz"
    with gzip.open(fname, "wt") as f:
        f.write(sample_pdb)

    yield fname


def test_load_pdb_file_regular(sample_pdb_as_file, sample_pdb):
    result = libnotebooks.load_pdb_file(str(sample_pdb_as_file))
    assert result == sample_pdb


def test_load_pdb_file_gz(sample_pdb_as_gz, sample_pdb):
    result = libnotebooks.load_pdb_file(str(sample_pdb_as_gz))
    assert result == sample_pdb


def test_load_pdb_file_not_found():
    result = libnotebooks.load_pdb_file("non_existent_file.pdb")
    assert result is None


def test_pdb_string_to_structure_and_structure_to_pdb_string(sample_pdb):
    # Convert PDB string to structure and back
    structure = libnotebooks.pdb_string_to_structure(sample_pdb, "test")
    assert structure.id == "test"

    out_pdb = libnotebooks.structure_to_pdb_string(structure)
    assert "ATOM" in out_pdb and "END" in out_pdb


@has_notebook
def test_align_full_success(tmp_path, sample_pdb):
    # Patch py3Dmol and create two similar minimal pdb files
    pdb1 = tmp_path / "model1.pdb"
    pdb2 = tmp_path / "model2.pdb"
    pdb1.write_text(sample_pdb)
    pdb2.write_text(sample_pdb)

    view = libnotebooks.align_full(
        str(pdb1), str(pdb2), chains=["A"], atom_types=["CA"], show_labels=False
    )
    # Check this is a view
    # NOTE: Yes this import needs to be here and not on top
    import py3Dmol

    assert isinstance(view, py3Dmol.view)


@has_notebook
def test_align_full_file_not_found():

    result = libnotebooks.align_full("nofile1.pdb", "nofile2.pdb")
    assert result[1] is None

    # NOTE: Yes this import needs to be here and not on top
    import py3Dmol

    assert isinstance(result[0], py3Dmol.view)


@pytest.fixture
def reference_pdb_file(tmp_path):
    fname = tmp_path / "reference.pdb"
    fname.write_text(REFERENCE_PDB)
    yield str(fname)


@pytest.fixture
def ensemble_pdb_file(tmp_path):
    fname = tmp_path / "ensemble.pdb"
    fname.write_text(ENSEMBLE_PDB)
    yield str(fname)


@pytest.fixture
def mock_view():
    """Replace `py3Dmol` by a mock and yield the mocked viewer."""
    view = MagicMock()
    mocked_py3dmol = MagicMock()
    mocked_py3dmol.view.return_value = view
    with patch.object(libnotebooks, "py3Dmol", mocked_py3dmol):
        yield view


@has_notebook
def test_align_full_ens_returns_view(reference_pdb_file, ensemble_pdb_file):
    view = libnotebooks.align_full_ens(
        reference_pdb_file, ensemble_pdb_file, chains=["A", "B"], atom_types=["CA"]
    )

    # NOTE: Yes this import needs to be here and not on top
    import py3Dmol

    assert isinstance(view, py3Dmol.view)


def test_align_full_ens_adds_reference_and_each_model(
    mock_view, reference_pdb_file, ensemble_pdb_file
):
    libnotebooks.align_full_ens(
        reference_pdb_file, ensemble_pdb_file, chains=["A", "B"], atom_types=["CA"]
    )

    # one model for the reference plus one per ensemble conformer
    assert mock_view.addModel.call_count == 3
    mock_view.addModelsAsFrames.assert_not_called()
    mock_view.animate.assert_not_called()

    # the reference is the first model added, the ensemble models follow
    ref_string = mock_view.addModel.call_args_list[0].args[0]
    assert "ALA A" in ref_string and "ALA B" in ref_string

    # each ensemble model is added on its own, aligned onto the reference
    for call in mock_view.addModel.call_args_list[1:]:
        assert call.args[0].count("CA") == 6

    # models 1 and 2 are styled as the ensemble
    styled_models = [
        call.args[0]["model"] for call in mock_view.setStyle.call_args_list
    ]
    assert styled_models == [0, 1, 2]


def test_align_full_ens_reports_rmsd_per_model(
    mock_view, reference_pdb_file, ensemble_pdb_file, capsys
):
    libnotebooks.align_full_ens(
        reference_pdb_file,
        ensemble_pdb_file,
        chains=["A", "B"],
        atom_types=["CA"],
        show_per_model_rmsd=True,
    )

    captured = capsys.readouterr().out
    # both ensemble models are rigid translations of the reference, so the
    # RMSD after superposition must be zero
    assert "Model 0: RMSD 0.000 Å (6 atom pairs)" in captured
    assert "Model 1: RMSD 0.000 Å (6 atom pairs)" in captured
    assert "Aligned 2 ensemble model(s) onto reference" in captured


def test_align_full_ens_hides_rmsd(
    mock_view, reference_pdb_file, ensemble_pdb_file, capsys
):
    libnotebooks.align_full_ens(
        reference_pdb_file,
        ensemble_pdb_file,
        chains=["A", "B"],
        atom_types=["CA"],
        show_per_model_rmsd=False,
    )

    assert "RMSD" not in capsys.readouterr().out


def test_align_full_ens_aligns_coordinates(
    mock_view, reference_pdb_file, ensemble_pdb_file
):
    libnotebooks.align_full_ens(
        reference_pdb_file, ensemble_pdb_file, chains=["A", "B"], atom_types=["CA"]
    )

    ref_struct = libnotebooks.pdb_string_to_structure(REFERENCE_PDB, "reference")
    ref_coords = [a.coord for a in ref_struct.get_atoms()]

    # the models handed over to the viewer must superimpose on the reference
    for call in mock_view.addModel.call_args_list[1:]:
        aligned = libnotebooks.pdb_string_to_structure(call.args[0], "aligned")
        for atom, ref_coord in zip(aligned.get_atoms(), ref_coords):
            assert atom.coord == pytest.approx(ref_coord, abs=1e-3)


def test_align_full_ens_animate(mock_view, reference_pdb_file, ensemble_pdb_file):
    libnotebooks.align_full_ens(
        reference_pdb_file,
        ensemble_pdb_file,
        chains=["A", "B"],
        atom_types=["CA"],
        animate=True,
        interval=500,
    )

    # only the reference is added as a model, the ensemble becomes frames
    assert mock_view.addModel.call_count == 1
    mock_view.addModelsAsFrames.assert_called_once()

    frames = mock_view.addModelsAsFrames.call_args.args[0]
    assert frames.count("CA") == 12

    mock_view.animate.assert_called_once_with(
        {"loop": "forward", "reps": 0, "interval": 500}
    )

    # the frames all live in a single viewer model
    styled_models = [
        call.args[0]["model"] for call in mock_view.setStyle.call_args_list
    ]
    assert styled_models == [0, 1]


def test_align_full_ens_animate_show_model_number(
    mock_view, reference_pdb_file, ensemble_pdb_file
):
    libnotebooks.align_full_ens(
        reference_pdb_file,
        ensemble_pdb_file,
        chains=["A", "B"],
        atom_types=["CA"],
        animate=True,
        show_model_number=True,
    )

    labels = [
        (call.args[0], call.args[1]["frame"])
        for call in mock_view.addLabel.call_args_list
        if "frame" in call.args[1]
    ]
    assert labels == [("Model 1", 0), ("Model 2", 1)]


def test_align_full_ens_no_labels_without_model_number(
    mock_view, reference_pdb_file, ensemble_pdb_file
):
    libnotebooks.align_full_ens(
        reference_pdb_file,
        ensemble_pdb_file,
        chains=["A", "B"],
        atom_types=["CA"],
        animate=True,
        show_model_number=False,
    )

    mock_view.addLabel.assert_not_called()


def test_align_full_ens_show_labels(mock_view, reference_pdb_file, ensemble_pdb_file):
    libnotebooks.align_full_ens(
        reference_pdb_file,
        ensemble_pdb_file,
        chains=["A", "B"],
        atom_types=["CA"],
        show_labels=True,
    )

    labels = [call.args[0] for call in mock_view.addLabel.call_args_list]
    assert labels == ["Reference", "Ensemble (2 models)", "Mean RMSD: 0.000 Å"]


def test_align_full_ens_file_not_found(mock_view, reference_pdb_file, capsys):
    view = libnotebooks.align_full_ens(reference_pdb_file, "nofile.pdb")

    assert view is mock_view
    mock_view.addModel.assert_not_called()
    assert "Failed to load" in capsys.readouterr().out


def test_align_full_ens_unknown_chain(
    mock_view, reference_pdb_file, ensemble_pdb_file, capsys
):
    libnotebooks.align_full_ens(
        reference_pdb_file, ensemble_pdb_file, chains=["Z"], atom_types=["CA"]
    )

    captured = capsys.readouterr().out
    assert "Could not find any reference atoms for alignment" in captured

    # both files are added unaligned
    assert mock_view.addModel.call_count == 2
    assert mock_view.addModel.call_args_list[0].args[0] == REFERENCE_PDB
    assert mock_view.addModel.call_args_list[1].args[0] == ENSEMBLE_PDB


def test_align_full_ens_skips_model_without_atoms(
    mock_view, reference_pdb_file, tmp_path, capsys
):
    # the second model of the ensemble has no chain B at all
    ensemble = ENSEMBLE_PDB.replace(
        "ATOM      4  CA  ALA B   1       0.000   7.000  10.000  1.00  0.00           C\n"
        "ATOM      5  CA  ALA B   2       3.800   7.000  10.000  1.00  0.00           C\n"
        "ATOM      6  CA  ALA B   3       3.800  10.800  10.000  1.00  0.00           C\n",
        "",
    )
    ensemble_path = tmp_path / "partial_ensemble.pdb"
    ensemble_path.write_text(ensemble)

    libnotebooks.align_full_ens(
        reference_pdb_file, str(ensemble_path), chains=["B"], atom_types=["CA"]
    )

    captured = capsys.readouterr().out
    assert "Model 1: no atoms for alignment, skipping" in captured
    assert "Aligned 1 ensemble model(s) onto reference" in captured
    # reference plus the single model that could be aligned
    assert mock_view.addModel.call_count == 2


def test_align_full_ens_alignment_failure(
    mock_view, reference_pdb_file, ensemble_pdb_file, capsys
):
    with patch.object(
        libnotebooks, "pdb_string_to_structure", side_effect=ValueError("boom")
    ):
        libnotebooks.align_full_ens(
            reference_pdb_file, ensemble_pdb_file, chains=["A"], atom_types=["CA"]
        )

    assert "Alignment failed: boom" in capsys.readouterr().out
    # falls back to displaying both files unaligned
    assert mock_view.addModel.call_count == 2
    assert mock_view.addModel.call_args_list[0].args[0] == REFERENCE_PDB
    assert mock_view.addModel.call_args_list[1].args[0] == ENSEMBLE_PDB


def test_model_to_pdb_string_single_model():
    structure = libnotebooks.pdb_string_to_structure(ENSEMBLE_PDB, "ensemble")
    first, second = structure.get_models()

    pdb_string = libnotebooks._model_to_pdb_string(second)

    assert pdb_string.startswith("ATOM")
    assert pdb_string.rstrip().endswith("END")
    # only the requested model is serialised
    assert pdb_string.count("CA") == 6
    assert "MODEL" not in pdb_string
    # coordinates of the first model are not present
    assert "5.000   0.000   0.000" not in pdb_string
    assert "0.000   7.000   0.000" in pdb_string


def test_model_to_pdb_string_roundtrip():
    structure = libnotebooks.pdb_string_to_structure(ENSEMBLE_PDB, "ensemble")
    model = next(structure.get_models())

    reparsed = libnotebooks.pdb_string_to_structure(
        libnotebooks._model_to_pdb_string(model), "reparsed"
    )

    original_atoms = list(model.get_atoms())
    reparsed_atoms = list(reparsed.get_atoms())
    assert len(reparsed_atoms) == len(original_atoms)
    for original, copy in zip(original_atoms, reparsed_atoms):
        assert copy.get_name() == original.get_name()
        assert (
            copy.get_parent().get_parent().id == original.get_parent().get_parent().id
        )
        assert copy.coord == pytest.approx(original.coord, abs=1e-3)


def test_model_to_pdb_string_keeps_chains():
    structure = libnotebooks.pdb_string_to_structure(ENSEMBLE_PDB, "ensemble")
    model = next(structure.get_models())

    pdb_string = libnotebooks._model_to_pdb_string(model)

    assert "ALA A" in pdb_string
    assert "ALA B" in pdb_string
