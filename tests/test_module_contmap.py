"""Test the CONTact MAP module."""
import os
import pytest
import pytest_mock  # noqa : F401
import tempfile
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist, squareform

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.contactmap import HaddockModule as ContactMapModule  # noqa : E501
from haddock.modules.analysis.contactmap import DEFAULT_CONFIG
from haddock.modules.analysis.contactmap.contmap import (
    ContactsMap,
    ClusteredContactMap,
    datakey_to_colorscale,
    write_res_contacts,
    min_dist,
    extract_submatrix,
    compute_distance_matrix,
    extract_pdb_coords,
    topX_models,
    gen_contact_dt,
    )

from . import golden_data


##########################
# Define pytest fixtures #
##########################
@pytest.fixture
def pdbline():
    return "ATOM     28  O   GLU A  21      11.097   3.208   5.136  1.00"


@pytest.fixture
def atoms_coordinates() -> list[list]:
    """List of 3 atom coordinates."""
    return [
        [1, 2, 3],
        [1, 2, 2],
        [0, 2, 3],
        ]


@pytest.fixture
def ref_one_line_dist_half_matrix():
    return [1, 1, np.sqrt(2)]


@pytest.fixture
def ref_dist_matrix():
    return np.array([
        [0, 1, 1],
        [1, 0, np.sqrt(2)],
        [1, np.sqrt(2), 0],
        ])


@pytest.fixture
def contactmap():
    """Return contmap module."""
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        yield ContactMapModule(
            order=1,
            path=Path(tmpdir),
            initial_params=DEFAULT_CONFIG,
            )


@pytest.fixture
def cluster_input_list() -> list:
    """Prot-prot input."""
    return [
        PDBFile(Path(golden_data, "protprot_complex_1.pdb"), path=golden_data),
        PDBFile(Path(golden_data, "protprot_complex_2.pdb"), path=golden_data),
        ]


@pytest.fixture
def cluster_input_iter(cluster_input_list) -> iter:
    """Generate an iterable of files path."""
    return iter(cluster_input_list)


@pytest.fixture
def params() -> dict:
    """Set of parameters."""
    return {
        "ca_ca_dist_threshold": 9.0,
        "shortest_dist_threshold": 7.5,
        "color_ramp": "Greys",
        "single_model_analysis": False,
        "generate_heatmap": True,
        "topX": 10,
        "cluster-heatmap-datatype": 'shortest-cont-ratio',
        }


@pytest.fixture
def protprot_contactmap(cluster_input_list, params):
    params["single_model_analysis"] = True
    return ContactsMap(
        Path(cluster_input_list[0].rel_path),
        Path('./contmap_test'),
        params,
        )


@pytest.fixture
def clustercontactmap(cluster_input_list, params):
    return ClusteredContactMap(
        [Path(m.rel_path) for m in cluster_input_list],
        Path('./clustcontmap_test'),
        params,
        )


@pytest.fixture
def res_res_contacts():
    return [
        {'res1': 'A-1-MET', 'res2': 'A-2-ALA', 'dist': 3.0},
        {'res1': 'A-1-MET', 'res2': 'A-3-VAL', 'dist': 4.0},
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
            PDBFile(Path(golden_data, "protprot_complex_1.pdb"),
                    path=golden_data),
            PDBFile(Path(golden_data, "protprot_complex_2.pdb"),
                    path=golden_data),
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
def test_contactmap_run_errors(contactmap, cluster_input_list, mocker):
    """Test content of _run() function from __init__.py HaddockModule class."""
    contactmap.previous_io = cluster_input_list
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


######################################################
# Testing of Classes and function withing contmap.py #
######################################################
def test_single_model(protprot_contactmap):
    """Test ContactsMap run function."""
    contacts_dt = protprot_contactmap.run()
    # check return variable
    assert type(contacts_dt) == list
    # check generated output files
    output_bp = protprot_contactmap.output
    assert os.path.exists(f'{output_bp}_contacts.tsv') is True
    assert Path(f'{output_bp}_contacts.tsv').stat().st_size != 0
    assert os.path.exists(f'{output_bp}_heatmap.html') is True
    assert Path(f'{output_bp}_heatmap.html').stat().st_size != 0
    Path(f'{output_bp}_contacts.tsv').unlink(missing_ok=False)
    Path(f'{output_bp}_heatmap.html').unlink(missing_ok=False)


def test_clustercontactmap_run(clustercontactmap):
    """Test ClusteredContactMap run function."""
    # run object
    clustercontactmap.run()
    # check terminated flag
    assert clustercontactmap.terminated is True
    # check outputs
    output_bp = clustercontactmap.output
    assert os.path.exists(f'{output_bp}_contacts.tsv') is True
    assert Path(f'{output_bp}_contacts.tsv').stat().st_size != 0
    assert os.path.exists(f'{output_bp}_heatmap.html') is True
    assert Path(f'{output_bp}_heatmap.html').stat().st_size != 0
    Path(f'{output_bp}_contacts.tsv').unlink(missing_ok=False)
    Path(f'{output_bp}_heatmap.html').unlink(missing_ok=False)


def test_write_res_contacts(res_res_contacts):
    """Test list of dict to tsv generation."""
    fpath = write_res_contacts(
        res_res_contacts,
        ['res1', 'res2', 'dist'],
        Path('./test-contacts.tsv'),
        )
    assert os.path.exists(fpath) is True
    with open(fpath, 'r') as filin:
        flines = filin.readlines()
    assert flines[0].strip().split('\t') == ['res1', 'res2', 'dist']
    assert flines[1].strip().split('\t') == ['A-1-MET', 'A-2-ALA', '3.0']
    assert flines[2].strip().split('\t') == ['A-1-MET', 'A-3-VAL', '4.0']
    fpath.unlink(missing_ok=True)


def test_topx_models(cluster_input_list):
    """Test sorting and X=1 topX function."""
    # set input PDBFiles scores
    cluster_input_list[0].score = -20.0
    cluster_input_list[1].score = -40.0  # `better` than -20.0
    topxmodels = topX_models(cluster_input_list, topX=1)
    assert len(topxmodels) == 1
    assert topxmodels[0].file_name == cluster_input_list[1].file_name


def test_topx_models_exception():
    """Test exception handling of non PDBFile type lists."""
    topxmodels = topX_models(['pdbpath1.pdb', 'pdbpath2.pdb'], topX=1)
    assert len(topxmodels) == 1
    assert topxmodels == ['pdbpath1.pdb']
    assert topxmodels[0] == 'pdbpath1.pdb'


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
    colorscale = 'color'
    new_colorscale = datakey_to_colorscale('ratio', color_scale=colorscale)
    assert new_colorscale == colorscale


def test_reverse_coloscale():
    """Test reversed color scale."""
    colorscale = 'color'
    new_colorscale = datakey_to_colorscale('anything', color_scale=colorscale)
    assert new_colorscale == f'{colorscale}_r'
    assert '_r' in new_colorscale


def test_gen_contact_dt_ca_exception(ref_dist_matrix):
    """Test missing ca coodrinates for a residue."""
    cont_dt = gen_contact_dt(
        ref_dist_matrix,
        {
            'r1': {'atoms_indices': [0, 1], 'CA': 0, 'resname': 'ALA'},
            'r2': {'atoms_indices': [2], 'resname': 'ALA'},
            },
        'r1',
        'r2',
        )
    assert cont_dt['ca-ca-dist'] == 9999
