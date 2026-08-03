"""Test the rnascan module."""

import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from haddock.core.exceptions import ConfigurationError
from haddock.libs.libontology import PDBFile
from haddock.libs.libscan import add_zscores
from haddock.modules.analysis.rnascan import DEFAULT_CONFIG
from haddock.modules.analysis.rnascan import HaddockModule as RnascanModule
from haddock.modules.analysis.rnascan.rnascan import (
    AddDeltaBFactor,
    calc_score,
    ClusterOutputer,
    get_atoms_to_keep,
    group_scan_by_cluster,
    mutate,
    validate_scan_bases,
    write_scan_out,
    MutationResult,
    InterfaceScanner,
    ModelPointMutation,
    BACKBONE_ATOMS,
    PURINE_BASE_ATOMS,
    PYRIMIDINE_BASE_ATOMS,
    DEFAULT_SCAN_BASES,
)

from . import golden_data


# The protein-RNA complex 1XS7 has an RNA chain B with residues:
#   2=U, 3=G, 4=U, 5=U, 6=U


# fixtures
@pytest.fixture(name="rna_input_list")
def fixture_rna_input_list():
    """Prot-RNA input (two entries pointing at the same complex)."""
    return [
        PDBFile(Path(golden_data, "1XS7_clean.pdb"), path=golden_data),
        PDBFile(Path(golden_data, "1XS7_clean.pdb"), path=golden_data),
    ]


@pytest.fixture(name="complex_pdb")
def fixture_complex_pdb():
    """Return example prot-RNA complex pdb."""
    return Path(golden_data, "1XS7_clean.pdb")


@pytest.fixture(name="rna_model_list")
def fixture_rna_model_list(complex_pdb):
    """Prot-RNA input."""
    return [
        PDBFile(file_name=complex_pdb, path=str(golden_data)),
    ]


@pytest.fixture(name="params")
def fixture_params():
    """Parameter fixture to be used in the tests"""
    return {
        "int_cutoff": 3.0,
        "plot": False,
        "scan_bases": ["A", "C", "G", "U"],
        "resdic_B": [3, 4],
        "chains": [],
        "output_mutants": False,
    }


@pytest.fixture(name="rnascan")
def fixture_rnascan(monkeypatch):
    """Return rnascan module."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        yield RnascanModule(
            order=1,
            path=Path("."),
            initial_params=DEFAULT_CONFIG,
        )


@pytest.fixture
def example_df_scan_clt():
    """Return example rnascan clt DataFrame."""
    example_clt_data = [
        ["B", 3, "G", "A", "B-3-G-A", -2.0, -1.0, -0.4, -2.3, -0.5, -7.2, 1.0],
        ["B", 3, "G", "U", "B-3-G-U", -0.0, 1.0, -0.4, 0.8, 0.5, -7.2, 1.0],
    ]
    columns = [
        "chain",
        "resid",
        "resname",
        "target_resname",
        "full_resname",
        "score",
        "delta_score",
        "delta_vdw",
        "delta_elec",
        "delta_desolv",
        "delta_bsa",
        "frac_pres",
    ]

    example_df_scan_clt = pd.DataFrame(example_clt_data, columns=columns)
    yield example_df_scan_clt


@pytest.fixture(name="interface_scanner")
def fixture_interface_scanner(rna_model_list, params):
    """Interface scanner in default HADDOCK mode."""
    return InterfaceScanner(
        model=rna_model_list[0],
        scan_bases=["A", "C", "G", "U"],
        params=params,
    )


@pytest.fixture(name="successful_mutation_result")
def fixture_successful_mutation_result():
    """Successful mutation result for B3 G->A."""
    return MutationResult(
        model_id="1XS7_clean",
        chain="B",
        resid=3,
        ori_resname="G",
        target_resname="A",
        mutant_scores=(-108.540, -43.234, -275.271, -10.252, 1589.260),
        delta_scores=(5.401, 0.119, 28.482, -0.414, 10.470),
        success=True,
    )


@pytest.fixture(name="successful_mutation_result_2")
def fixture_successful_mutation_result_2():
    """Successful mutation result for B3 G->U."""
    return MutationResult(
        model_id="1XS7_clean",
        chain="B",
        resid=3,
        ori_resname="G",
        target_resname="U",
        mutant_scores=(-18.540, -41.234, -270.271, -11.252, 589.260),
        delta_scores=(4.401, 0.109, 27.482, -0.3, 0.470),
        success=True,
    )


@pytest.fixture(name="results_by_model")
def fixture_results_by_model(successful_mutation_result, successful_mutation_result_2):
    return {
        "1XS7_clean": [
            successful_mutation_result,
            successful_mutation_result_2,
        ],
    }


# tests start here
def test_init(rnascan):
    rnascan.__init__(
        order=42,
        path=Path("0_anything"),
        initial_params=DEFAULT_CONFIG,
    )
    assert rnascan.path == Path("0_anything")
    assert rnascan._origignal_config_file == DEFAULT_CONFIG
    assert isinstance(rnascan.params, dict)
    assert len(rnascan.params) != 0
    # scan_bases default must contain all four RNA bases
    assert set(rnascan.params["scan_bases"]) == {"A", "C", "G", "U"}


def test_confirm_installation(rnascan):
    assert rnascan.confirm_installation() is None


"""Test validate_scan_bases guard"""


def test_validate_scan_bases_canonical():
    """Canonical RNA names are accepted and returned unchanged."""
    assert validate_scan_bases(["A", "C", "G", "U"]) == ["A", "C", "G", "U"]


def test_validate_scan_bases_case_insensitive_and_dedup():
    """Lower-case entries are normalised and duplicates removed, order kept."""
    assert validate_scan_bases(["g", "A", "g"]) == ["G", "A"]


def test_validate_scan_bases_rejects_dna_names():
    """Two-letter (DNA) names must be rejected."""
    with pytest.raises(ConfigurationError):
        validate_scan_bases(["DA", "DC", "DG", "DT"])


def test_validate_scan_bases_rejects_thymine():
    """Thymine (T, a DNA-only base) must be rejected."""
    with pytest.raises(ConfigurationError):
        validate_scan_bases(["A", "T"])


def test_validate_scan_bases_rejects_unknown():
    """Unknown bases must be rejected."""
    with pytest.raises(ConfigurationError):
        validate_scan_bases(["A", "X"])


def test_validate_scan_bases_rejects_empty():
    """An empty list must be rejected."""
    with pytest.raises(ConfigurationError):
        validate_scan_bases([])


"""Test InterfaceScanner"""


def test_interface_scanner_init_default_bases():
    """Default scan_bases should be the four RNA bases."""
    scanner = InterfaceScanner(model=Path("test.pdb"))
    assert scanner.scan_bases == DEFAULT_SCAN_BASES
    assert scanner.params == {}
    assert scanner.filter_resdic == {}


def test_interface_scanner_run_four_mutations(mocker, interface_scanner):
    """Each interface nucleotide yields one job per non-WT base."""
    mocker.patch(
        "haddock.modules.analysis.rnascan.rnascan.calc_score",
        return_value=(-106.7, -29.6, -316.5, -13.8, 1494.7),
    )
    mocker.patch(
        "haddock.libs.libcapri.CAPRI.identify_interface",
        return_value={"B": [3, 4]},
    )
    mocker.patch("haddock.libs.libalign.get_atoms")
    mocker.patch(
        "haddock.libs.libalign.load_coords",
        return_value=(
            {
                ("B", 3, "C1'", "G"): [0, 0, 0],
                ("B", 4, "C1'", "U"): [1, 1, 1],
            }
        ),
    )
    mutation_jobs = interface_scanner.run()
    # 2 nucleotides x (4 bases - WT) = 6 jobs
    assert len(mutation_jobs) == 6
    assert all(isinstance(job, ModelPointMutation) for job in mutation_jobs)
    # WT base must never appear as a target
    for job in mutation_jobs:
        assert job.ori_resname != job.target_resname
    # the guanine (B3) should be mutated to A, C, U (not G)
    g_targets = {j.target_resname for j in mutation_jobs if j.resid == 3}
    assert g_targets == {"A", "C", "U"}


def test_interface_scanner_skips_non_rna(mocker, rna_model_list, params):
    """Protein residues at the interface must be ignored."""
    scanner = InterfaceScanner(
        model=rna_model_list[0],
        params={**params, "resdic_B": [3], "resdic_A": [40]},
    )
    mocker.patch(
        "haddock.modules.analysis.rnascan.rnascan.calc_score",
        return_value=(-106.7, -29.6, -316.5, -13.8, 1494.7),
    )
    mocker.patch(
        "haddock.libs.libcapri.CAPRI.identify_interface",
        return_value={"A": [40], "B": [3]},
    )
    mocker.patch("haddock.libs.libalign.get_atoms")
    mocker.patch(
        "haddock.libs.libalign.load_coords",
        return_value=(
            {
                ("A", 40, "CA", "THR"): [0, 0, 0],
                ("B", 3, "C1'", "G"): [1, 1, 1],
            }
        ),
    )
    mutation_jobs = scanner.run()
    # Both A40 (protein) and B3 (RNA) are selected by resdic, but only the RNA
    # nucleotide must be scanned; the protein residue is skipped.
    assert {job.chain for job in mutation_jobs} == {"B"}
    assert all(job.ori_resname == "G" for job in mutation_jobs)


def test_interface_scanner_restrict_scan_bases(mocker, rna_model_list, params):
    """Restricting scan_bases limits the number of mutations."""
    scanner = InterfaceScanner(
        model=rna_model_list[0],
        scan_bases=["G"],
        params=params,
    )
    mocker.patch(
        "haddock.modules.analysis.rnascan.rnascan.calc_score",
        return_value=(-106.7, -29.6, -316.5, -13.8, 1494.7),
    )
    mocker.patch(
        "haddock.libs.libcapri.CAPRI.identify_interface",
        return_value={"B": [3, 4]},
    )
    mocker.patch("haddock.libs.libalign.get_atoms")
    mocker.patch(
        "haddock.libs.libalign.load_coords",
        return_value=(
            {
                ("B", 3, "C1'", "G"): [0, 0, 0],
                ("B", 4, "C1'", "U"): [1, 1, 1],
            }
        ),
    )
    mutation_jobs = scanner.run()
    # G3 -> G skipped (WT), U4 -> G (1 job) => only 1 job
    assert len(mutation_jobs) == 1
    assert mutation_jobs[0].resid == 4
    assert mutation_jobs[0].target_resname == "G"


"""Test stand-alone functions"""


def test_get_atoms_to_keep_same_ring_type():
    """Same-type mutations keep the shared base ring atoms."""
    # purine -> purine
    keep = get_atoms_to_keep("G", "A")
    assert set(BACKBONE_ATOMS).issubset(keep)
    assert set(PURINE_BASE_ATOMS).issubset(keep)
    assert not set(PYRIMIDINE_BASE_ATOMS).issubset(keep)
    # pyrimidine -> pyrimidine
    keep = get_atoms_to_keep("U", "C")
    assert set(PYRIMIDINE_BASE_ATOMS).issubset(keep)
    # N9/C8 (purine-only) must not be kept for a pyrimidine mutation
    assert "N9" not in keep
    assert "C8" not in keep


def test_get_atoms_to_keep_cross_ring_type():
    """Cross-type mutations keep the backbone plus renamed anchor atoms."""
    # purine -> pyrimidine: keep glycosidic anchors, renamed N9->N1, C4->C2, C8->C6
    keep = get_atoms_to_keep("G", "C")
    assert set(BACKBONE_ATOMS).issubset(keep)
    assert keep["N9"] == "N1"
    assert keep["C4"] == "C2"
    assert keep["C8"] == "C6"
    # only the three anchor atoms are kept from the base
    assert set(keep) - set(BACKBONE_ATOMS) == {"N9", "C4", "C8"}
    # pyrimidine -> purine: keep glycosidic anchors, renamed N1->N9, C2->C4, C6->C8
    keep = get_atoms_to_keep("U", "A")
    assert set(BACKBONE_ATOMS).issubset(keep)
    assert keep["N1"] == "N9"
    assert keep["C2"] == "C4"
    assert keep["C6"] == "C8"
    assert set(keep) - set(BACKBONE_ATOMS) == {"N1", "C2", "C6"}


def test_backbone_includes_ribose_hydroxyl():
    """The RNA backbone must retain the ribose 2'-hydroxyl oxygen."""
    assert "O2'" in BACKBONE_ATOMS


def test_mutate_cross_type(rna_model_list, monkeypatch):
    """A purine->pyrimidine mutation keeps only backbone atoms."""
    with tempfile.TemporaryDirectory() as tmpdir:
        mut_fname = Path(golden_data, rna_model_list[0].file_name)
        monkeypatch.chdir(path=tmpdir)
        # residue B3 is a guanine (purine); U is a pyrimidine -> cross-type
        mut_pdb_fname = mutate(mut_fname, "B", 3, "U")

        assert mut_pdb_fname == Path("1XS7_clean-B_G3U.pdb")
        assert os.path.exists(mut_pdb_fname)

        kept_atoms = []
        with open(mut_pdb_fname) as fh:
            for line in fh:
                if (
                    line.startswith("ATOM")
                    and line[21] == "B"
                    and int(line[22:26]) == 3
                ):
                    assert line[17:20] == "  U"
                    kept_atoms.append(line[12:16].strip().replace("*", "'"))
        # purine -> pyrimidine: the glycosidic anchors are kept and renamed
        # (N9 -> N1, C4 -> C2, C8 -> C6); no other base atoms remain
        assert "N9" not in kept_atoms
        assert "C8" not in kept_atoms
        assert "N1" in kept_atoms
        assert "C2" in kept_atoms
        assert "C6" in kept_atoms
        # only the three renamed anchors survive from the base
        base_atoms = [a for a in kept_atoms if a not in BACKBONE_ATOMS]
        assert set(base_atoms) == {"N1", "C2", "C6"}
        assert "C1'" in kept_atoms
        assert "P" in kept_atoms
        # ribose 2'-hydroxyl must be preserved
        assert "O2'" in kept_atoms

        with pytest.raises(KeyError):
            mutate(mut_fname, "B", 3, "HOH")


def test_mutate_same_type_keeps_ring(rna_model_list, monkeypatch):
    """A purine->purine mutation keeps the shared purine ring atoms."""
    with tempfile.TemporaryDirectory() as tmpdir:
        mut_fname = Path(golden_data, rna_model_list[0].file_name)
        monkeypatch.chdir(path=tmpdir)
        # residue B3 is a guanine (purine); A is also a purine -> same-type
        mut_pdb_fname = mutate(mut_fname, "B", 3, "A")

        kept_atoms = []
        with open(mut_pdb_fname) as fh:
            for line in fh:
                if (
                    line.startswith("ATOM")
                    and line[21] == "B"
                    and int(line[22:26]) == 3
                ):
                    assert line[17:20] == "  A"
                    kept_atoms.append(line[12:16].strip().replace("*", "'"))
        # shared purine ring atoms must be preserved to keep base orientation
        for atom in PURINE_BASE_ATOMS:
            assert atom in kept_atoms
        # guanine-specific atoms O6/N2 must be dropped (rebuilt as N6 for A)
        assert "O6" not in kept_atoms
        assert "N2" not in kept_atoms
        assert "C1'" in kept_atoms


def test_write_delta_score_to_pdb(rna_model_list, results_by_model, monkeypatch):
    """Test the write_delta_score_to_pdb method with per-residue averaging."""
    with tempfile.TemporaryDirectory() as tmpdir:
        mut_fname = Path(golden_data, rna_model_list[0].file_name)
        monkeypatch.chdir(path=tmpdir)
        mut_pdb_fname = mutate(mut_fname, "B", 3, "A")
        rna_model_list[0].file_name = mut_pdb_fname
        rna_model_list[0].path = tmpdir
        add_delta_to_bfactor_obj = AddDeltaBFactor(
            rna_model_list[0],
            tmpdir,
            results_by_model["1XS7_clean"],
        )
        add_delta_to_bfactor_obj.reorder_results()
        # Both results are for residue B3; the stored score is their mean.
        assert np.isclose(
            add_delta_to_bfactor_obj.model_results["B-3"],
            np.mean([5.401, 4.401]),
        )


def test_add_zscores(example_df_scan_clt):
    """Test the add_zscores function."""
    obs_df_scan_clt = add_zscores(example_df_scan_clt)
    assert np.isclose(obs_df_scan_clt["z_score"].values[0], -1.0)
    assert np.isclose(obs_df_scan_clt["z_score"].values[1], 1.0)


def test_calc_score(mocker):
    """Test the calc_score function."""
    mocker.patch(
        "haddock.libs.libscan.get_score_string",
        return_value=[
            "> starting calculations...",
            "> HADDOCK-score = (1.0 * vdw) + (0.2 * elec) + (1.0 * desolv) + (0.0 * air) + (0.0 * bsa)",
            "> HADDOCK-score (emscoring) = -106.7376",
            "> vdw=-29.5808,elec=-316.542,desolv=-13.8484,air=0.0,bsa=1494.73",
        ],
    )
    scores = calc_score(Path(golden_data, "1XS7_clean.pdb"), run_dir="tmp")
    assert scores == (-106.7376, -29.5808, -316.542, -13.8484, 1494.73)


def test_rnascan_group_scan_by_cluster(rna_input_list, results_by_model, monkeypatch):
    """Group scan data per cluster, keyed by mutation (base included)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        clt_scan, clt_pops = group_scan_by_cluster(rna_input_list, results_by_model)
        assert "unclustered" in clt_scan.keys()
        assert clt_pops["unclustered"] == 2
        # The two mutations of nucleotide B3 must be reported separately
        idents = set(clt_scan["unclustered"].keys())
        assert "B-3-G-A" in idents
        assert "B-3-G-U" in idents


def test_rnascan_cluster_full_outputs(rna_input_list, results_by_model, monkeypatch):
    """Test ClusterOutputer writes tsv and plot."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        clt_scan, clt_pops = group_scan_by_cluster(rna_input_list, results_by_model)
        ClusterOutputer(
            clt_scan["unclustered"],
            "unclustered",
            clt_pops["unclustered"],
            scan_residue="RNA base",
            generate_plot=True,
            offline=False,
        ).run()
        assert Path("scan_clt_unclustered.tsv").exists()
        assert Path("scan_clt_unclustered.html").exists()


def test_write_scan_out_with_mutation_results(
    successful_mutation_result, successful_mutation_result_2, monkeypatch
):
    """Test write_scan_out."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        results = [successful_mutation_result, successful_mutation_result_2]
        write_scan_out(results, "1XS7_clean")
        output_file = Path("scan_1XS7_clean.tsv")
        assert output_file.exists()
        content = output_file.read_text()
        assert "# `rnascan` results for 1XS7_clean" in content
        assert "delta_score" in content
        assert "z_score" in content
        df = pd.read_csv(output_file, sep="\t", comment="#")
        assert len(df) == 2
        assert df["ori_resname"].tolist() == ["G", "G"]
        assert set(df["end_resname"].tolist()) == {"A", "U"}
