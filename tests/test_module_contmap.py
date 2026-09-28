"""Test the CONTact MAP module."""

import os
import tempfile
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial.distance import pdist, squareform

from haddock.libs.libontology import PDBFile
from haddock.libs.libutil import get_available_memory
from haddock.modules.analysis.contactmap import DEFAULT_CONFIG
from haddock.modules.analysis.contactmap import HaddockModule as ContactMapModule
from haddock.modules.analysis.contactmap.contmap import (
    CONNECT_COLORS,
    OTHER_PAIR_COLOR,
    PI,
    ClusteredContactMap,
    ContactsMap,
    check_square_matrix,
    compute_residue_distances,
    control_pts,
    ctrl_rib_chords,
    datakey_to_colorscale,
    extract_heavyatom_contacts,
    extract_pdb_coords,
    extract_pdb_dt,
    extract_submatrix,
    gen_contacts_dt,
    get_all_ideograms_ends,
    get_cont_type,
    get_data_threshold,
    get_ordered_coords,
    get_pair_color,
    group_indices_by_chain,
    is_hydrogen,
    make_chordchart,
    make_contactmap_report,
    make_ideogram_arc,
    make_q_bezier,
    make_ribbon_arc,
    moduloAB,
    parse_reskey,
    split_labels_by_chains,
    to_color_weight,
    topX_models,
    tsv_to_chordchart,
    tsv_to_heatmap,
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
def fixture_contactmap(monkeypatch):
    """Return contmap module."""
    with tempfile.TemporaryDirectory() as tempdir:
        monkeypatch.chdir(tempdir)
        yield ContactMapModule(
            order=1,
            path=Path("."),
            initial_params=DEFAULT_CONFIG,
        )


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
def fixture_protprot_contactmap(protprot_input_list, params, monkeypatch):
    """???"""
    params["single_model_analysis"] = True
    with (
        tempfile.TemporaryDirectory() as tempdir,
        tempfile.NamedTemporaryFile() as temp_f,
    ):
        monkeypatch.chdir(tempdir)
        return ContactsMap(
            model=Path(protprot_input_list[0].rel_path),
            output=Path(temp_f.name),
            params=params,
        )


@pytest.fixture(name="clustercontactmap")
def fixture_clustercontactmap(protprot_input_list, params, monkeypatch):
    """???"""
    with (
        tempfile.TemporaryDirectory() as tempdir,
        tempfile.NamedTemporaryFile() as temp_f,
    ):
        monkeypatch.chdir(tempdir)
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


def test_write_res_contacts(res_res_contacts, monkeypatch):
    """Test list of dict to tsv generation."""

    with (
        tempfile.TemporaryDirectory() as tempdir,
        tempfile.NamedTemporaryFile(suffix=".tsv") as temp_f,
    ):
        monkeypatch.chdir(tempdir)
        res_contact_output_f = write_res_contacts(
            res_res_contacts=res_res_contacts,
            header=["res1", "res2", "ca-ca-dist"],
            path=Path(temp_f.name),
        )
        assert os.path.exists(res_contact_output_f) is True

        with open(res_contact_output_f, "r", encoding="utf-8") as fh:
            lines = [_ for _ in fh if not _.startswith("#")]

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


def test_topx_models_undefined_scores(protprot_input_list):
    """Test that undefined scores keep the input order."""
    protprot_input_list[0].score = None
    protprot_input_list[1].score = -40.0
    topxmodels = topX_models(protprot_input_list, topX=1)
    assert topxmodels == [protprot_input_list[0]]


def test_compute_residue_distances(protprot_input_list):
    """Test block computation of residue distances against brute force."""
    pdb_dt = extract_pdb_dt(Path(protprot_input_list[0].rel_path))
    coords, resid_keys, resid_dt = get_ordered_coords(pdb_dt)
    full = squareform(pdist(coords))
    # Small blocks to go through several iterations
    ref_dists, shortest_dists = compute_residue_distances(
        coords, resid_keys, resid_dt, max_block=len(coords) * 50
    )
    for i in (0, 5, len(resid_keys) - 1):
        for j in (1, 7, len(resid_keys) - 2):
            idx_i = resid_dt[resid_keys[i]]["atoms_indices"]
            idx_j = resid_dt[resid_keys[j]]["atoms_indices"]
            expected = full[np.ix_(idx_i, idx_j)].min()
            assert np.isclose(shortest_dists[i, j], expected)
            ref_i = resid_dt[resid_keys[i]]["ref"]
            ref_j = resid_dt[resid_keys[j]]["ref"]
            assert np.isclose(ref_dists[i, j], full[ref_i, ref_j])


def test_extract_heavyatom_contacts(protprot_input_list):
    """Test KD-tree interchain contacts against brute force."""
    pdb_dt = extract_pdb_dt(Path(protprot_input_list[0].rel_path))
    coords, resid_keys, resid_dt = get_ordered_coords(pdb_dt)
    contacts = extract_heavyatom_contacts(
        coords, resid_keys, resid_dt, contact_distance=4.5
    )
    full = squareform(pdist(coords))
    chains = np.array(
        [resid_dt[k]["chainid"] for k in resid_keys for _ in resid_dt[k]["atoms_order"]]
    )
    expected = (full <= 4.5) & (chains[:, None] != chains[None, :])
    assert len(contacts) == np.triu(expected, k=1).sum()
    assert len(contacts) > 0
    for cont in contacts:
        assert (
            parse_reskey(cont["atom1"].rsplit("-", 1)[0])[0]
            != (parse_reskey(cont["atom2"].rsplit("-", 1)[0])[0])
        )
        assert cont["dist"] <= 4.5


def test_extract_submatrix_empty_columns(ref_dist_matrix):
    """Test that an empty list of columns is not treated as None."""
    submat = extract_submatrix(ref_dist_matrix, [1, 2], indices2=[])
    assert submat.shape == (2, 0)


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


def _atom_line(
    atname: str,
    resname: str,
    chain: str,
    resseq: int,
    coords: tuple,
    record: str = "ATOM",
    altloc: str = " ",
    icode: str = " ",
    element: str = "",
) -> str:
    """Build a fixed width PDB ATOM/HETATM record."""
    name = f" {atname:<3s}" if len(atname) < 4 else atname
    x, y, z = coords
    return (
        f"{record:<6s}{1:5d} {name}{altloc}{resname:>3s} {chain}{resseq:4d}{icode}"
        f"   {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n"
    )


def test_extract_pdb_dt_insertion_altloc_model(tmp_path):
    """Test insertion codes, alternate locations, hydrogens and models."""
    lines = [
        "MODEL        1\n",
        _atom_line("N", "ALA", "A", 52, (0, 0, 0)),
        _atom_line("CA", "ALA", "A", 52, (1, 0, 0)),
        _atom_line("CA", "GLY", "A", 52, (5, 0, 0), icode="A"),
        _atom_line("CB", "SER", "A", 53, (9, 0, 0), altloc="A"),
        _atom_line("CB", "SER", "A", 53, (9, 9, 9), altloc="B"),
        _atom_line("1HB", "SER", "A", 53, (9, 1, 0)),
        _atom_line("HG", "HG", "B", 1, (20, 0, 0), record="HETATM"),
        _atom_line("CA", "CA", "B", 2, (30, 0, 0), record="HETATM"),
        _atom_line("CA", "ALA", "A", -5, (40, 0, 0)),
        "ENDMDL\n",
        "MODEL        2\n",
        _atom_line("CA", "ALA", "C", 1, (0, 0, 0)),
        "ENDMDL\n",
    ]
    pdb = tmp_path / "test.pdb"
    pdb.write_text("".join(lines))
    pdb_dt = extract_pdb_dt(pdb)
    # Second model ignored
    assert pdb_dt["chain_order"] == ["A", "B"]
    # Insertion code kept in residue identifier
    assert pdb_dt["A"]["order"] == ["52", "52A", "53", "-5"]
    # Only the first alternate location is kept, hydrogen skipped
    assert pdb_dt["A"]["53"]["atoms_order"] == ["CB"]
    assert pdb_dt["A"]["53"]["atoms"]["CB"] == [9.0, 0.0, 0.0]
    # Mercury is not a hydrogen
    assert pdb_dt["B"]["1"]["atoms_order"] == ["HG"]

    coords, resid_keys, resid_dt = get_ordered_coords(pdb_dt)
    assert resid_keys == [
        "A-52-ALA",
        "A-52A-GLY",
        "A-53-SER",
        "A--5-ALA",
        "B-1-HG",
        "B-2-CA",
    ]
    assert coords.shape == (7, 3)
    # CA of ALA 52 is a reference atom (N present), calcium is not
    assert resid_dt["A-52-ALA"]["ref"] == 1
    assert resid_dt["B-2-CA"]["ref"] is None

    ref_dists, _ = compute_residue_distances(coords, resid_keys, resid_dt)
    contacts = gen_contacts_dt(ref_dists, ref_dists, resid_keys, resid_dt)
    calcium_contacts = [c for c in contacts if c["res2"] == "B-2-CA"]
    assert all(np.isnan(c["ca-ca-dist"]) for c in calcium_contacts)


def test_is_hydrogen():
    """Test hydrogen detection."""
    assert is_hydrogen(_atom_line("HA", "ALA", "A", 1, (0, 0, 0)))
    assert is_hydrogen(_atom_line("2HG1", "VAL", "A", 1, (0, 0, 0)))
    assert not is_hydrogen(_atom_line("CA", "ALA", "A", 1, (0, 0, 0)))
    assert not is_hydrogen(_atom_line("HG", "HG", "A", 1, (0, 0, 0)))
    # Element column takes precedence over the atom name
    assert not is_hydrogen(_atom_line("HG", "XXX", "A", 1, (0, 0, 0), element="HG"))
    assert is_hydrogen(_atom_line("X1", "XXX", "A", 1, (0, 0, 0), element="H"))


def test_parse_reskey():
    """Test residue key parsing, including negative residue numbers."""
    assert parse_reskey("A-12-ALA") == ("A", "12", "ALA")
    assert parse_reskey("A--5-ALA") == ("A", "-5", "ALA")
    assert parse_reskey("B-52A-GLY") == ("B", "52A", "GLY")
    assert split_labels_by_chains(["A--5-ALA", "B-1-GLY"]) == {
        "A": ["A--5-ALA"],
        "B": ["B-1-GLY"],
    }


def test_residue_classes():
    """Test residue class pairs."""
    assert get_cont_type("LYS", "DA") == "positive-nucleotide"
    assert get_cont_type("SEP", "HIP") == "negative-positive"
    assert get_cont_type("ZN", "ALA") == "unknown-apolar"
    assert get_cont_type("NAG", "TRP") == "carbohydrate-polar"
    assert get_cont_type("SIA", "BGC") == "carbohydrate-carbohydrate"


def test_carbohydrate_reference_atom(tmp_path):
    """Test that C1 is the reference atom of carbohydrates only."""
    lines = [
        _atom_line("N", "ALA", "A", 1, (0, 0, 0)),
        _atom_line("CA", "ALA", "A", 1, (1, 0, 0)),
        _atom_line("C1", "BGC", "B", 1, (5, 0, 0)),
        _atom_line("C2", "BGC", "B", 1, (6, 0, 0)),
        # A ligand with a C1 atom has no reference atom
        _atom_line("C1", "LIG", "C", 1, (9, 0, 0), record="HETATM"),
    ]
    pdb = tmp_path / "glycan.pdb"
    pdb.write_text("".join(lines))
    coords, resid_keys, resid_dt = get_ordered_coords(extract_pdb_dt(pdb))
    assert resid_dt["B-1-BGC"]["ref"] == 2
    assert resid_dt["C-1-LIG"]["ref"] is None
    ref_dists, _ = compute_residue_distances(coords, resid_keys, resid_dt)
    assert np.isclose(ref_dists[0, 1], 4.0)
    assert np.isnan(ref_dists[0, 2])


def test_get_pair_color():
    """Test order independent pair colors."""
    assert get_pair_color("positive-negative") == CONNECT_COLORS["negative-positive"]
    assert get_pair_color("negative-positive") == CONNECT_COLORS["negative-positive"]
    assert get_pair_color("polar-apolar") == OTHER_PAIR_COLOR
    assert get_pair_color("self-self") == OTHER_PAIR_COLOR


def test_to_color_weight_bounds():
    """Test that color weights are clamped to their range."""
    assert to_color_weight(1.0, 9.0) == 0.9
    assert to_color_weight(9.0, 9.0) == 0.35
    assert to_color_weight(12.0, 9.0) == 0.35
    assert 0.35 < to_color_weight(5.0, 9.0) < 0.9


def test_get_data_threshold(params):
    """Test thresholds used for each data type."""
    assert get_data_threshold("shortest-cont-probability", params) == 1.0
    assert get_data_threshold("ca-ca-cont-probability", params) == 1.0
    assert get_data_threshold("ca-ca-dist", params) == 9.0
    assert get_data_threshold("shortest-dist", params) == 7.5


def test_ideograms_follow_labels_chain_order():
    """Test that residue arcs follow label order, not alphabetical order."""
    chains = split_labels_by_chains(["B-1-ALA", "B-2-ALA", "A-1-GLY"])
    assert list(chains) == ["B", "A"]
    ideo_ends, chain_ideo_ends = get_all_ideograms_ends(chains)
    # First two residues arcs lie within the first chain arc (chain B)
    for start, end in ideo_ends[:2]:
        assert chain_ideo_ends[0][0] <= start < end <= chain_ideo_ends[0][1] + 1e-9
    start, end = ideo_ends[2]
    assert chain_ideo_ends[1][0] <= start < end <= chain_ideo_ends[1][1] + 1e-9


def test_group_indices_by_chain():
    """Test that non contiguous chains are grouped."""
    labels = ["B-1-ALA", "A-1-GLY", "B-2-ALA"]
    assert group_indices_by_chain(labels) == [0, 2, 1]


def _write_contacts(path, contacts):
    return write_res_contacts(
        contacts,
        ["res1", "res2", "ca-ca-dist", "contact-type", "shortest-dist"],
        path,
    )


def test_tsv_to_chordchart_no_interchain(tmp_path):
    """Test that no chord chart is generated without interchain contacts."""
    tsv = _write_contacts(
        tmp_path / "contacts.tsv",
        [
            {
                "res1": "A-1-ALA",
                "res2": "A-2-ALA",
                "ca-ca-dist": 3.8,
                "contact-type": "apolar-apolar",
                "shortest-dist": 1.3,
            }
        ],
    )
    output = tmp_path / "chord.html"
    assert tsv_to_chordchart(tsv, output_fname=output) is None
    assert not output.exists()


def test_tsv_to_heatmap_distance_not_clipped_to_one(tmp_path, mocker):
    """Test that distances are bounded to the given threshold only."""
    tsv = _write_contacts(
        tmp_path / "contacts.tsv",
        [
            {
                "res1": "A-1-ALA",
                "res2": "B-1-ALA",
                "ca-ca-dist": 5.0,
                "contact-type": "apolar-apolar",
                "shortest-dist": 3.0,
            }
        ],
    )
    heatmap = mocker.patch(
        "haddock.modules.analysis.contactmap.contmap.heatmap_plotly",
        return_value="heatmap.html",
    )
    tsv_to_heatmap(tsv, data_key="shortest-dist", contact_threshold=4.5)
    matrix = heatmap.call_args[0][0]
    assert matrix[0, 1] == 3.0


def test_make_contactmap_report_escapes_and_filters(tmp_path, monkeypatch):
    """Test that the report lists only job files and escapes them."""
    monkeypatch.chdir(tmp_path)
    job = ClusteredContactMap([Path("model_1.pdb")], Path("cluster<1>"), {})
    Path("cluster<1>_contacts.tsv").write_text("")
    # Per-model file not reported without single model analysis
    Path("cluster<1>_model_1_contacts.tsv").write_text("")
    make_contactmap_report([job], "report.html")
    report = Path("report.html").read_text()
    assert "cluster&lt;1&gt;_contacts.tsv" in report
    assert "model_1" not in report
    assert "<1>" not in report


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


def test_compute_residue_distances_missing_reference(atoms_coordinates):
    """Test missing reference atom for a residue gives nan distances."""
    resid_dt = {
        "A-1-ALA": {"atoms_indices": [0, 1], "ref": 0, "resname": "ALA"},
        "B-1-ZN": {"atoms_indices": [2], "ref": None, "resname": "ZN"},
    }
    resid_keys = list(resid_dt)
    ref_dists, shortest_dists = compute_residue_distances(
        np.asarray(atoms_coordinates, dtype=float), resid_keys, resid_dt
    )
    contacts = gen_contacts_dt(ref_dists, shortest_dists, resid_keys, resid_dt)
    assert len(contacts) == 1
    assert np.isnan(contacts[0]["ca-ca-dist"])
    assert contacts[0]["shortest-dist"] == 1.0
    assert contacts[0]["contact-type"] == "apolar-unknown"


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
    monkeypatch,
):
    """Test main function."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
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


def test_get_available_memory():
    """Test get_available_memory function."""
    memory = get_available_memory()
    # Should return a positive float (system memory)
    assert isinstance(memory, float)
    assert memory > 0
