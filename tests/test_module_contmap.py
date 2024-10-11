"""Test the CONTact MAP module."""

import os
import tempfile
from pathlib import Path
from typing import Callable

import numpy as np
import pytest
from scipy.spatial.distance import pdist, squareform

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.contactmap import DEFAULT_CONFIG
from haddock.modules.analysis.contactmap import \
    HaddockModule as ContactMapModule
from haddock.modules.analysis.contactmap.contmap import (
    PI,
    ClusteredContactMap,
    ContactsMap,
    check_square_matrix,
    compute_distance_matrix,
    control_pts,
    ctrl_rib_chords,
    datakey_to_colorscale,
    extract_pdb_coords,
    extract_submatrix,
    gen_contact_dt,
    invPerm,
    make_chordchart,
    make_ideogram_arc,
    make_q_bezier,
    make_ribbon_arc,
    min_dist,
    moduloAB,
    topX_models,
    within_2PI,
    write_res_contacts,
    )

from . import golden_data


##########################
# Define pytest fixtures #
##########################
@pytest.fixture(name="contactmap_output_ext")
def fixture_contactmap_output_ext():
    """List of generated files suffixes."""
    return (
        "contacts.tsv",
        "heatmap.html",
        "chordchart.html",
        "interchain_contacts.tsv",
        "heavyatoms_interchain_contacts.tsv",
    )


@pytest.fixture(name="pdbline")
def fixture_pdbline():
    """???"""
    return "ATOM     28  O   GLU A  21      11.097   3.208   5.136  1.00"


@pytest.fixture(name="atoms_coordinates")
def fixture_atoms_coordinates() -> list[list]:
    """List of 3 atom coordinates."""
    return [
        [1, 2, 3],
        [1, 2, 2],
        [0, 2, 3],
    ]


@pytest.fixture(name="ref_one_line_dist_half_matrix")
def fixture_ref_one_line_dist_half_matrix():
    """???"""
    return [1, 1, np.sqrt(2)]


@pytest.fixture(name="ref_dist_matrix")
def fixture_ref_dist_matrix():
    """???"""
    return np.array(
        [
            [0, 1, 1],
            [1, 0, np.sqrt(2)],
            [1, np.sqrt(2), 0],
        ]
    )


@pytest.fixture(name="contactmap")
def fixture_contactmap():
    """Return contmap module."""
    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)
        yield ContactMapModule(
            order=1,
            path=Path("."),
            initial_params=DEFAULT_CONFIG,
        )


# @pytest.fixture
# def protprot_input_list() -> list:
#     """Prot-prot input."""
#     return [
#         PDBFile(Path(golden_data, "protprot_complex_1.pdb"), path=golden_data),
#         PDBFile(Path(golden_data, "protprot_complex_2.pdb"), path=golden_data),
#     ]


@pytest.fixture(name="cluster_input_iter")
def fixture_cluster_input_iter(protprot_input_list) -> Callable:
    """Generate an iterable of files path."""
    return iter(protprot_input_list)


@pytest.fixture(name="params")
def fixture_params() -> dict:
    """Set of parameters."""
    return {
        "ca_ca_dist_threshold": 9.0,
        "shortest_dist_threshold": 7.5,
        "color_ramp": "Greys",
        "single_model_analysis": False,
        "topX": 10,
        "generate_heatmap": True,
        "cluster_heatmap_datatype": "shortest-cont-probability",
        "generate_chordchart": True,
        "chordchart_datatype": "shortest-dist",
        "offline": False,
    }


@pytest.fixture(name="protprot_contactmap")
def fixture_protprot_contactmap(protprot_input_list, params):
    """???"""
    params["single_model_analysis"] = True
    with (
        tempfile.TemporaryDirectory() as tempdir,
        tempfile.NamedTemporaryFile() as temp_f,
    ):
        os.chdir(tempdir)
        return ContactsMap(
            model=Path(protprot_input_list[0].rel_path),
            output=Path(temp_f.name),
            params=params,
        )


@pytest.fixture(name="clustercontactmap")
def fixture_clustercontactmap(protprot_input_list, params):
    """???"""
    with (
        tempfile.TemporaryDirectory() as tempdir,
        tempfile.NamedTemporaryFile() as temp_f,
    ):
        os.chdir(tempdir)
        return ClusteredContactMap(
            models=[Path(m.rel_path) for m in protprot_input_list],
            output=Path(temp_f.name),
            params=params,
        )


@pytest.fixture(name="res_res_contacts")
def fixture_res_res_contacts():
    """???"""
    return [
        {"res1": "A-1-MET", "res2": "A-2-ALA", "ca-ca-dist": 3.0},
        {"res1": "A-1-MET", "res2": "A-3-VAL", "ca-ca-dist": 4.0},
    ]


@pytest.fixture(name="contact_matrix")
def fixture_contact_matrix():
    """???"""
    return np.array(
        [
            [1, 0, 0, 1, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [1, 0, 0, 1, 1, 0],
            [0, 0, 0, 1, 1, 1],
            [0, 0, 0, 0, 1, 1],
        ]
    )


@pytest.fixture(name="dist_matrix")
def fixture_dist_matrix():
    """???"""
    return np.array(
        [
            [0.0, 18.4, 16.0, 5.1, 11.4, 14.7],
            [18.4, 0.0, 18.0, 15.2, 11.7, 20.1],
            [16.0, 18.0, 0.0, 9.2, 12.8, 10.0],
            [5.1, 15.2, 9.2, 0.0, 6.2, 9.3],
            [11.4, 11.7, 12.8, 6.2, 0.0, 6.9],
            [14.7, 20.1, 10.0, 9.3, 6.9, 0.0],
        ]
    )


@pytest.fixture(name="intertype_matrix")
def fixture_intertype_matrix():
    """???"""
    return np.array(
        [
            [
                "self-self",
                "polar-positive",
                "polar-apolar",
                "polar-apolar",
                "polar-apolar",
                "polar-negative",
            ],
            [
                "polar-positive",
                "self-self",
                "positive-apolar",
                "positive-apolar",
                "positive-apolar",
                "positive-negative",
            ],
            [
                "polar-apolar",
                "positive-apolar",
                "self-self",
                "apolar-apolar",
                "apolar-apolar",
                "apolar-negative",
            ],
            [
                "polar-apolar",
                "positive-apolar",
                "apolar-apolar",
                "self-self",
                "apolar-apolar",
                "apolar-negative",
            ],
            [
                "polar-apolar",
                "positive-apolar",
                "apolar-apolar",
                "apolar-apolar",
                "self-self",
                "apolar-negative",
            ],
            [
                "polar-negative",
                "positive-negative",
                "apolar-negative",
                "apolar-negative",
                "apolar-negative",
                "self-self",
            ],
        ]
    )


@pytest.fixture(name="protein_labels")
def fixture_protein_labels():
    """???"""
    return [
        "A-69-THR",
        "A-255-LYS",
        "A-296-ILE",
        "A-323-LEU",
        "B-44-GLY",
        "B-70-GLU",
    ]


########################
# Test Scipy functions #
########################
def test_pdist(atoms_coordinates, ref_one_line_dist_half_matrix):
    """Test computation of distance matrix."""
    one_line_dist_half_matrix = pdist(atoms_coordinates)
    assert one_line_dist_half_matrix.tolist() == ref_one_line_dist_half_matrix


def test_squareform(ref_one_line_dist_half_matrix, ref_dist_matrix):
    """Test scipy on line half matrix to squareform behavior."""
    dist_matrix = squareform(ref_one_line_dist_half_matrix)
    assert np.array_equal(dist_matrix, ref_dist_matrix)


def test_numpy_set_diagonal(ref_dist_matrix):
    """Test numpy fil_diagonal function."""
    np.fill_diagonal(ref_dist_matrix, 1)
    assert ref_dist_matrix[0, 0] == 1
    assert ref_dist_matrix[1, 1] == 1
    assert ref_dist_matrix[2, 2] == 1


###################################
# Testing of class in __init__.py #
###################################
def test_confirm_installation(contactmap):
    """Test confirm install."""
    assert contactmap.confirm_installation() is None


def test_init(contactmap):
    """Test __init__ function."""
    contactmap.__init__(
        order=42,
        path=Path("0_anything"),
        initial_params=DEFAULT_CONFIG,
    )

    # Once a module is initialized, it should have the following attributes
    assert contactmap.path == Path("0_anything")
    assert contactmap._origignal_config_file == DEFAULT_CONFIG
    assert type(contactmap.params) == dict
    assert len(contactmap.params.keys()) != 0


class MockPreviousIO:
    """A mocking class holding specific methods."""

    # In the mocked method, add the arguments that are called by the original
    #  method that is being tested
    def retrieve_models(self, individualize: bool = False):
        """Provide a set of models."""
        models = [
            PDBFile(Path(golden_data, "protprot_complex_1.pdb"), path=golden_data),
            PDBFile(Path(golden_data, "protprot_complex_2.pdb"), path=golden_data),
        ]
        models[0].clt_id = None
        models[1].clt_id = 1
        return models


def test_contactmap_run(contactmap, mocker):
    """Test content of _run() function from __init__.py HaddockModule class."""
    # Mock some functions
    contactmap.previous_io = MockPreviousIO()
    mocker.patch("haddock.libs.libparallel.Scheduler.run", return_value=None)
    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models",
        return_value=None,
    )
    # run main module _run() function
    module_sucess = contactmap.run()
    assert module_sucess is None


#########################################
# Testing of previous_io errors handles #
#########################################
def test_contactmap_run_errors(contactmap, protprot_input_list, mocker):
    """Test content of _run() function from __init__.py HaddockModule class."""
    contactmap.previous_io = protprot_input_list
    # run main module _run() function
    with pytest.raises(RuntimeError):
        module_sucess = contactmap.run()
        assert module_sucess is None


def test_contactmap_run_iter_errors(contactmap, cluster_input_iter, mocker):
    """Test content of _run() function from __init__.py HaddockModule class."""
    contactmap.previous_io = cluster_input_iter
    # run main module _run() function.
    with pytest.raises(RuntimeError):
        module_sucess = contactmap.run()
        assert module_sucess is None


#####################################################
# Testing of Classes and function within contmap.py #
#####################################################
def test_single_model(protprot_contactmap, contactmap_output_ext):
    """Test ContactsMap run function."""
    contacts_dt = protprot_contactmap.run()
    # check return variable
    assert type(contacts_dt) == tuple
    assert type(contacts_dt[0]) == list
    assert type(contacts_dt[1]) == list
    # check generated output files
    output_bp = protprot_contactmap.output
    for output_ext in contactmap_output_ext:
        fpath = f"{output_bp}_{output_ext}"
        assert os.path.exists(fpath) is True
        assert Path(fpath).stat().st_size != 0
        Path(fpath).unlink(missing_ok=False)


def test_clustercontactmap_run(clustercontactmap, contactmap_output_ext):
    """Test ClusteredContactMap run function."""
    # run object
    clustercontactmap.run()
    # check terminated flag
    assert clustercontactmap.terminated is True
    # check outputs
    output_bp = clustercontactmap.output
    for output_ext in contactmap_output_ext:
        fpath = f"{output_bp}_{output_ext}"
        assert os.path.exists(fpath) is True
        assert Path(fpath).stat().st_size != 0
        Path(fpath).unlink(missing_ok=False)


def test_write_res_contacts(res_res_contacts):
    """Test list of dict to tsv generation."""

    with (
        tempfile.TemporaryDirectory() as tempdir,
        tempfile.NamedTemporaryFile(suffix=".tsv") as temp_f,
    ):
        os.chdir(tempdir)
        res_contact_output_f = write_res_contacts(
            res_res_contacts=res_res_contacts,
            header=["res1", "res2", "ca-ca-dist"],
            path=Path(temp_f.name),
        )
        assert os.path.exists(res_contact_output_f) is True

        with open(res_contact_output_f, "r", encoding="utf-8") as fh:
            lines = [_ for _ in fh.readlines() if not _.startswith("#")]

        assert lines[0].strip().split("\t") == ["res1", "res2", "ca-ca-dist"]
        assert lines[1].strip().split("\t") == ["A-1-MET", "A-2-ALA", "3.0"]
        assert lines[2].strip().split("\t") == ["A-1-MET", "A-3-VAL", "4.0"]


def test_topx_models(protprot_input_list):
    """Test sorting and X=1 topX function."""
    # set input PDBFiles scores
    protprot_input_list[0].score = -20.0
    protprot_input_list[1].score = -40.0  # `better` than -20.0
    topxmodels = topX_models(protprot_input_list, topX=1)
    assert len(topxmodels) == 1
    assert topxmodels[0].file_name == protprot_input_list[1].file_name


def test_topx_models_exception():
    """Test exception handling of non PDBFile type lists."""
    topxmodels = topX_models(["pdbpath1.pdb", "pdbpath2.pdb"], topX=1)
    assert len(topxmodels) == 1
    assert topxmodels == ["pdbpath1.pdb"]
    assert topxmodels[0] == "pdbpath1.pdb"


def test_compute_distance_matrix(atoms_coordinates, ref_dist_matrix):
    """Test compuation of all vs all distance matrix from coodrinates."""
    dist_matrix = compute_distance_matrix(atoms_coordinates)
    assert dist_matrix.tolist() == ref_dist_matrix.tolist()


def test_extract_submatrix_symetric(ref_dist_matrix):
    """Test extraction of submatrix with symetrical indices."""
    submat = extract_submatrix(ref_dist_matrix, [1, 2])
    assert np.array_equal(submat, np.array([[0, np.sqrt(2)], [np.sqrt(2), 0]]))


def test_extract_submatrix(ref_dist_matrix):
    """Test extraction of submatrix with non-symetrical indices."""
    submatrix = extract_submatrix(ref_dist_matrix, [1, 2], indices2=[0, 1])
    assert np.array_equal(submatrix, np.array([[1, 0], [1, np.sqrt(2)]]))


def test_extract_pdb_coords(pdbline):
    """Test to extract PDB atom coordinates."""
    atom_coods = extract_pdb_coords(pdbline)
    assert atom_coods == [11.097, 3.208, 5.136]


def test_min_dist(ref_dist_matrix):
    """Test extraction of lowest value in ND array."""
    assert min_dist(ref_dist_matrix) == 0


def test_no_reverse_coloscale():
    """Test non reversed color scale."""
    colorscale = "color"
    new_colorscale = datakey_to_colorscale(
        "probability",
        color_scale=colorscale,
    )
    assert new_colorscale == colorscale


def test_reverse_coloscale():
    """Test reversed color scale."""
    colorscale = "color"
    new_colorscale = datakey_to_colorscale("anything", color_scale=colorscale)
    assert new_colorscale == f"{colorscale}_r"
    assert "_r" in new_colorscale


def test_gen_contact_dt_ca_exception(ref_dist_matrix):
    """Test missing ca coodrinates for a residue."""
    cont_dt = gen_contact_dt(
        ref_dist_matrix,
        {
            "r1": {"atoms_indices": [0, 1], "CA": 0, "resname": "ALA"},
            "r2": {"atoms_indices": [2], "resname": "ALA"},
        },
        "r1",
        "r2",
    )
    assert cont_dt["ca-ca-dist"] == 9999


#################################
# Testing chord chart functions #
#################################
def test_moduloAB_errors():
    """Test argument error by moduloAB()."""
    with pytest.raises(ValueError):
        noreturn = moduloAB(1.0, 3.1, 2.6)
        assert noreturn is None


def test_moduloAB_execution():
    """Test proper functioning of moduloAB()."""
    moduloab = moduloAB(6.1035, 0, 2 * PI)
    assert np.isclose(moduloab, 6.1035, atol=0.001)
    moduloab2 = moduloAB(7, 0, 2 * PI)
    assert np.isclose(moduloab2, 0.7168, atol=0.001)
    moduloab3 = moduloAB(3, 0, 2)
    assert moduloab3 == 1


def test_within_2PI():
    """Test of within [0, 2 * Pi) range test function."""
    assert within_2PI(7) is False
    assert within_2PI(2 * PI) is False
    assert within_2PI(3)


def test_check_square_matrix_shape():
    """Test 2D array shape error."""
    with pytest.raises(ValueError):
        noreturn = check_square_matrix(np.array([[[1, 2]]]))
        assert noreturn is None


def test_check_square_matrix_error():
    """Test square matrix error."""
    with pytest.raises(ValueError):
        noreturn = check_square_matrix(np.array([[1, 2], [1, 2], [1, 2]]))
        assert noreturn is None


def test_check_square_matrix(contact_matrix):
    """Test square matrix."""
    length = check_square_matrix(contact_matrix)
    assert length == 6


def test_make_q_bezier_error():
    """Test error raising in make_q_bezier()."""
    with pytest.raises(ValueError):
        noreturn = make_q_bezier([1.0, 2.0])
        assert noreturn is None


def test_make_ribbon_arc_angle_error():
    """Test error raising in make_ribbon_arc()."""
    with pytest.raises(ValueError):
        noreturn = make_ribbon_arc(-1.0, 1.0)
        assert noreturn is None
        noreturn2 = make_ribbon_arc(1.0, -1.0)
        assert noreturn2 is None


def test_make_ribbon_arc_angle_baddef_error():
    """Test error raising in make_ribbon_arc()."""
    with pytest.raises(ValueError):
        noreturn = make_ribbon_arc(1.1, 1.2)
        assert noreturn is None


def test_make_chordchart(
    contact_matrix,
    dist_matrix,
    intertype_matrix,
    protein_labels,
):
    """Test main function."""
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        outputpath = "chord.html"
        graphpath = make_chordchart(
            contact_matrix,
            dist_matrix,
            intertype_matrix,
            protein_labels,
            output_fpath=outputpath,
            title="test",
        )
        assert graphpath == outputpath
        assert os.path.exists(outputpath)
        assert Path(outputpath).stat().st_size != 0
        Path(outputpath).unlink(missing_ok=False)


def test_invPerm():
    """Test permutation inversion."""
    inv = invPerm([1, 2, 3, 0])
    assert inv == [3, 0, 1, 2]


def test_ctrl_rib_chords_error():
    """Test error raising in ctrl_rib_chords()."""
    with pytest.raises(ValueError):
        noreturn = ctrl_rib_chords([1, 2, 3], [1, 2], 1.2)
        assert noreturn is None


def test_control_pts_error():
    """Test error raising in control_pts()."""
    with pytest.raises(ValueError):
        noreturn = control_pts([1, 0], 1.2)
        assert noreturn is None


def test_make_ideogram_arc_moduloAB():
    """Test usage of moduloAB while providing outranged angle values."""
    nb_points = 2
    arc_positions = make_ideogram_arc(1.1, (1, -1), nb_points=nb_points)
    excpected_output = np.array(
        [
            0.5943325364549538 + 0.9256180832886862j,
            0.5943325364549535 - 0.9256180832886863j,
        ]
    )
    assert arc_positions.shape == excpected_output.shape
    for i in range(nb_points):
        assert np.isclose(arc_positions[i], excpected_output[i], atol=0.0001)
