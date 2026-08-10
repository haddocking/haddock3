"""Test the dnascan module."""

import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from haddock.core.exceptions import ConfigurationError
from haddock.libs.libontology import PDBFile
from haddock.libs.libscan import add_zscores
from haddock.modules.analysis.dnascan import DEFAULT_CONFIG
from haddock.modules.analysis.dnascan import HaddockModule as DnascanModule
from haddock.modules.analysis.dnascan.dnascan import (
    AddDeltaBFactor,
    calc_score,
    ClusterOutputer,
    COMPLEMENT,
    find_base_pairs,
    get_atoms_to_keep,
    group_scan_by_cluster,
    is_cross_type,
    mutate,
    validate_scan_bases,
    wc_atom,
    write_scan_out,
    MutationResult,
    InterfaceScanner,
    ModelBasePairMutation,
    BACKBONE_ATOMS,
    PURINE_BASE_ATOMS,
    PYRIMIDINE_BASE_ATOMS,
    DEFAULT_SCAN_BASES,
)

from . import golden_data


# The protein-DNA complex protdna_complex_1 has a protein chain A and a
# double-stranded DNA chain B. A few relevant base pairs (chain B):
#   2=DG paired with 37=DC ; 7=DA paired with 32=DT ; 5=DC paired with 34=DG


# fixtures
@pytest.fixture(name="dna_input_list")
def fixture_dna_input_list():
    """Prot-DNA input (two entries pointing at the same complex)."""
    return [
        PDBFile(Path(golden_data, "protdna_complex_1.pdb"), path=golden_data),
        PDBFile(Path(golden_data, "protdna_complex_1.pdb"), path=golden_data),
    ]


@pytest.fixture(name="complex_pdb")
def fixture_complex_pdb():
    """Return example prot-DNA complex pdb."""
    return Path(golden_data, "protdna_complex_1.pdb")


@pytest.fixture(name="dna_model_list")
def fixture_dna_model_list(complex_pdb):
    """Prot-DNA input."""
    return [
        PDBFile(file_name=complex_pdb, path=str(golden_data)),
    ]


@pytest.fixture(name="params")
def fixture_params():
    """Parameter fixture to be used in the tests"""
    return {
        "int_cutoff": 5.0,
        "bp_cutoff": 3.5,
        "plot": False,
        "scan_bases": ["DA", "DC", "DG", "DT"],
        # resdic_ with the trailing underscore mimics the default (no user
        # residue filter), so that the whole interface is scanned.
        "resdic_": [],
        "chains": [],
        "output_mutants": False,
    }


@pytest.fixture(name="dnascan")
def fixture_dnascan(monkeypatch):
    """Return dnascan module."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        yield DnascanModule(
            order=1,
            path=Path("."),
            initial_params=DEFAULT_CONFIG,
        )


@pytest.fixture
def example_df_scan_clt():
    """Return example dnascan clt DataFrame."""
    example_clt_data = [
        ["B", 2, "DG", "DA", "B-2-DG-DA", -2.0, -1.0, -0.4, -2.3, -0.5, -7.2, 1.0],
        ["B", 2, "DG", "DT", "B-2-DG-DT", -0.0, 1.0, -0.4, 0.8, 0.5, -7.2, 1.0],
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
def fixture_interface_scanner(dna_model_list, params):
    """Interface scanner in default HADDOCK mode."""
    return InterfaceScanner(
        model=dna_model_list[0],
        scan_bases=["DA", "DC", "DG", "DT"],
        params=params,
    )


@pytest.fixture(name="successful_mutation_result")
def fixture_successful_mutation_result():
    """Successful base-pair mutation result for B2 DG->DA / B37 DC->DT."""
    return MutationResult(
        model_id="protdna_complex_1",
        chain="B",
        resid=2,
        ori_resname="DG",
        target_resname="DA",
        partner_chain="B",
        partner_resid=37,
        partner_ori_resname="DC",
        partner_target_resname="DT",
        mutant_scores=(-108.540, -43.234, -275.271, -10.252, 1589.260),
        delta_scores=(5.401, 0.119, 28.482, -0.414, 10.470),
        success=True,
    )


@pytest.fixture(name="successful_mutation_result_2")
def fixture_successful_mutation_result_2():
    """Successful base-pair mutation result for B2 DG->DC / B37 DC->DG."""
    return MutationResult(
        model_id="protdna_complex_1",
        chain="B",
        resid=2,
        ori_resname="DG",
        target_resname="DC",
        partner_chain="B",
        partner_resid=37,
        partner_ori_resname="DC",
        partner_target_resname="DG",
        mutant_scores=(-18.540, -41.234, -270.271, -11.252, 589.260),
        delta_scores=(4.401, 0.109, 27.482, -0.3, 0.470),
        success=True,
    )


@pytest.fixture(name="results_by_model")
def fixture_results_by_model(successful_mutation_result, successful_mutation_result_2):
    return {
        "protdna_complex_1": [
            successful_mutation_result,
            successful_mutation_result_2,
        ],
    }


# tests start here
def test_init(dnascan):
    dnascan.__init__(
        order=42,
        path=Path("0_anything"),
        initial_params=DEFAULT_CONFIG,
    )
    assert dnascan.path == Path("0_anything")
    assert dnascan._origignal_config_file == DEFAULT_CONFIG
    assert isinstance(dnascan.params, dict)
    assert len(dnascan.params) != 0
    # scan_bases default must contain all four DNA bases
    assert set(dnascan.params["scan_bases"]) == {"DA", "DC", "DG", "DT"}


def test_confirm_installation(dnascan):
    assert dnascan.confirm_installation() is None


"""Test validate_scan_bases guard"""


def test_validate_scan_bases_canonical():
    """Canonical DNA names are accepted and returned unchanged."""
    assert validate_scan_bases(["DA", "DC", "DG", "DT"]) == ["DA", "DC", "DG", "DT"]


def test_validate_scan_bases_case_insensitive_and_dedup():
    """Lower-case/mixed entries are normalised and duplicates removed."""
    assert validate_scan_bases(["dg", "DA", "DG"]) == ["DG", "DA"]


def test_validate_scan_bases_rejects_one_letter_rna_names():
    """One-letter names (RNA bases) must be rejected."""
    with pytest.raises(ConfigurationError):
        validate_scan_bases(["A", "C", "G", "T"])
    with pytest.raises(ConfigurationError):
        validate_scan_bases(["DA", "U"])


def test_validate_scan_bases_rejects_unknown():
    """Unknown bases must be rejected."""
    with pytest.raises(ConfigurationError):
        validate_scan_bases(["DA", "X"])


def test_validate_scan_bases_rejects_empty():
    """An empty list must be rejected."""
    with pytest.raises(ConfigurationError):
        validate_scan_bases([])


"""Test base-pair detection"""


def test_wc_atom():
    """Purines pair through N1, pyrimidines through N3."""
    assert wc_atom("DA") == "N1"
    assert wc_atom("DG") == "N1"
    assert wc_atom("DC") == "N3"
    assert wc_atom("DT") == "N3"


def test_find_base_pairs_real_duplex(complex_pdb):
    """Watson-Crick base pairs are detected in a real DNA duplex."""
    from haddock.libs.libalign import get_atoms, load_coords

    atoms = get_atoms(complex_pdb)
    coords, _ = load_coords(complex_pdb, atoms, add_resname=True)
    pairs = find_base_pairs(coords)
    # canonical pairs of the duplex
    assert pairs[("B", 2)] == ("B", 37, "DC")  # DG:DC
    assert pairs[("B", 7)] == ("B", 32, "DT")  # DA:DT
    # pairing must be symmetric
    assert pairs[("B", 37)] == ("B", 2, "DG")


def test_find_base_pairs_only_complementary_ring_types():
    """A purine can only be paired with a pyrimidine."""
    coords = {
        # two purines close together must NOT be paired to each other
        ("B", 1, "N1", "DA"): np.array([0.0, 0.0, 0.0]),
        ("B", 2, "N1", "DG"): np.array([0.0, 0.0, 2.0]),
    }
    assert find_base_pairs(coords) == {}


"""Test InterfaceScanner"""


def test_interface_scanner_init_default_bases():
    """Default scan_bases should be the four DNA bases."""
    scanner = InterfaceScanner(model=Path("test.pdb"))
    assert scanner.scan_bases == DEFAULT_SCAN_BASES
    assert scanner.params == {}
    assert scanner.filter_resdic == {}


def test_interface_scanner_double_mutations(mocker, interface_scanner):
    """Each interface base pair yields double-mutation jobs with a WC partner."""
    mocker.patch(
        "haddock.modules.analysis.dnascan.dnascan.calc_score",
        return_value=(-106.7, -29.6, -316.5, -13.8, 1494.7),
    )
    # A protein residue (A44) plus two DNA nucleotides that pair with each
    # other (B7:B32) and one DNA nucleotide whose partner (B37) is not itself
    # at the interface.
    mocker.patch(
        "haddock.libs.libcapri.CAPRI.identify_interface",
        return_value={"A": [44], "B": [2, 7, 32]},
    )
    mutation_jobs = interface_scanner.run()

    # all jobs are base-pair (double) mutations
    assert all(isinstance(job, ModelBasePairMutation) for job in mutation_jobs)
    # the protein residue must be ignored: only DNA chain B is mutated
    assert {job.chain for job in mutation_jobs} == {"B"}
    # WT base must never appear as a target and partner must be its complement
    for job in mutation_jobs:
        assert job.ori_resname != job.target_resname
        assert job.partner_target_resname == COMPLEMENT[job.target_resname]
    # B7 and B32 form a single base pair: it must be scanned only once, from
    # the first-encountered nucleotide (B7). B32 must never be a primary.
    assert 32 not in {job.resid for job in mutation_jobs}
    # base pair B2:B37 -> 3 jobs (DG to DA/DC/DT); base pair B7:B32 -> 3 jobs
    assert len(mutation_jobs) == 6
    # B2 (DG) mutated to A, C, T (not G); partner B37 always mutated
    b2_targets = {j.target_resname for j in mutation_jobs if j.resid == 2}
    assert b2_targets == {"DA", "DC", "DT"}
    assert all(j.partner_resid == 37 for j in mutation_jobs if j.resid == 2)


def test_interface_scanner_skips_unpaired(mocker, dna_model_list, params):
    """DNA nucleotides without a base-pairing partner are skipped."""
    scanner = InterfaceScanner(
        model=dna_model_list[0],
        params=params,
    )
    mocker.patch(
        "haddock.modules.analysis.dnascan.dnascan.calc_score",
        return_value=(-106.7, -29.6, -316.5, -13.8, 1494.7),
    )
    mocker.patch("haddock.modules.analysis.dnascan.dnascan.get_atoms")
    # a single isolated nucleotide with no complementary partner
    mocker.patch(
        "haddock.modules.analysis.dnascan.dnascan.load_coords",
        return_value=(
            {("B", 5, "N1", "DG"): np.array([0.0, 0.0, 0.0])},
            {},
        ),
    )
    mocker.patch(
        "haddock.libs.libcapri.CAPRI.identify_interface",
        return_value={"B": [5]},
    )
    mutation_jobs = scanner.run()
    assert mutation_jobs == []


def test_interface_scanner_restrict_scan_bases(mocker, dna_model_list, params):
    """Restricting scan_bases limits the number of mutations."""
    scanner = InterfaceScanner(
        model=dna_model_list[0],
        scan_bases=["DA"],
        params=params,
    )
    mocker.patch(
        "haddock.modules.analysis.dnascan.dnascan.calc_score",
        return_value=(-106.7, -29.6, -316.5, -13.8, 1494.7),
    )
    mocker.patch(
        "haddock.libs.libcapri.CAPRI.identify_interface",
        return_value={"B": [2]},
    )
    mutation_jobs = scanner.run()
    # B2 is DG -> only the DG->DA mutation remains (partner DC->DT)
    assert len(mutation_jobs) == 1
    assert mutation_jobs[0].resid == 2
    assert mutation_jobs[0].target_resname == "DA"
    assert mutation_jobs[0].partner_target_resname == "DT"


"""Test stand-alone functions"""


def test_get_atoms_to_keep_same_ring_type():
    """Same-type mutations keep the shared base ring atoms."""
    # purine -> purine
    keep = get_atoms_to_keep("DG", "DA")
    assert set(BACKBONE_ATOMS).issubset(keep)
    assert set(PURINE_BASE_ATOMS).issubset(keep)
    assert not set(PYRIMIDINE_BASE_ATOMS).issubset(keep)
    # pyrimidine -> pyrimidine
    keep = get_atoms_to_keep("DT", "DC")
    assert set(PYRIMIDINE_BASE_ATOMS).issubset(keep)
    # N9/C8 (purine-only) must not be kept for a pyrimidine mutation
    assert "N9" not in keep
    assert "C8" not in keep


def test_get_atoms_to_keep_cross_ring_type():
    """Cross-type mutations keep the backbone plus renamed anchor atoms."""
    # purine -> pyrimidine: keep glycosidic anchors, renamed N9->N1, C4->C2, C8->C6
    keep = get_atoms_to_keep("DG", "DC")
    assert set(BACKBONE_ATOMS).issubset(keep)
    assert keep["N9"] == "N1"
    assert keep["C4"] == "C2"
    assert keep["C8"] == "C6"
    assert set(keep) - set(BACKBONE_ATOMS) == {"N9", "C4", "C8"}
    # pyrimidine -> purine: keep glycosidic anchors, renamed N1->N9, C2->C4, C6->C8
    keep = get_atoms_to_keep("DT", "DA")
    assert set(BACKBONE_ATOMS).issubset(keep)
    assert keep["N1"] == "N9"
    assert keep["C2"] == "C4"
    assert keep["C6"] == "C8"
    assert set(keep) - set(BACKBONE_ATOMS) == {"N1", "C2", "C6"}


def test_backbone_excludes_ribose_hydroxyl():
    """The DNA (deoxyribose) backbone must NOT contain the 2'-hydroxyl O2'."""
    assert "O2'" not in BACKBONE_ATOMS
    assert "C2'" in BACKBONE_ATOMS


def test_mutate_double_cross_type(dna_model_list, monkeypatch):
    """A base-pair mutation mutates both nucleotides at once."""
    with tempfile.TemporaryDirectory() as tmpdir:
        mut_fname = Path(golden_data, dna_model_list[0].file_name)
        monkeypatch.chdir(path=tmpdir)
        # B2 is DG (purine); mutate to DC (pyrimidine, cross-type). Partner B37
        # is DC and must be mutated to the complement of DC, i.e. DG.
        mut_pdb_fname = mutate(mut_fname, "B", 2, "DC", "B", 37, "DG")

        assert mut_pdb_fname == Path("protdna_complex_1-B_G2C-C37G.pdb")
        assert os.path.exists(mut_pdb_fname)

        primary_atoms, partner_atoms = [], []
        with open(mut_pdb_fname) as fh:
            for line in fh:
                if not line.startswith("ATOM") or line[21] != "B":
                    continue
                resid = int(line[22:26])
                atom = line[12:16].strip().replace("*", "'")
                if resid == 2:
                    assert line[17:20] == " DC"  # DG -> DC
                    primary_atoms.append(atom)
                elif resid == 37:
                    assert line[17:20] == " DG"  # DC -> DG
                    partner_atoms.append(atom)
        # cross-type primary DG -> DC (purine -> pyrimidine): the glycosidic
        # anchors are kept and renamed N9 -> N1, C4 -> C2, C8 -> C6
        assert "N9" not in primary_atoms
        assert "C8" not in primary_atoms
        assert set(a for a in primary_atoms if a not in BACKBONE_ATOMS) == {
            "N1",
            "C2",
            "C6",
        }
        assert "C1'" in primary_atoms
        # cross-type partner DC -> DG (pyrimidine -> purine): anchors renamed
        # N1 -> N9, C2 -> C4, C6 -> C8
        assert set(a for a in partner_atoms if a not in BACKBONE_ATOMS) == {
            "N9",
            "C4",
            "C8",
        }
        # DNA backbone has no 2'-hydroxyl
        assert "O2'" not in primary_atoms


def test_mutate_double_same_type_keeps_ring(dna_model_list, monkeypatch):
    """A purine->purine base-pair mutation keeps the shared ring atoms."""
    with tempfile.TemporaryDirectory() as tmpdir:
        mut_fname = Path(golden_data, dna_model_list[0].file_name)
        monkeypatch.chdir(path=tmpdir)
        # B2 DG -> DA (purine->purine) ; partner B37 DC -> DT (pyrimidine)
        mut_pdb_fname = mutate(mut_fname, "B", 2, "DA", "B", 37, "DT")

        primary_atoms = []
        with open(mut_pdb_fname) as fh:
            for line in fh:
                if (
                    line.startswith("ATOM")
                    and line[21] == "B"
                    and int(line[22:26]) == 2
                ):
                    assert line[17:20] == " DA"
                    primary_atoms.append(line[12:16].strip().replace("*", "'"))
        # shared purine ring atoms preserved to keep base orientation
        for atom in PURINE_BASE_ATOMS:
            assert atom in primary_atoms
        # guanine-specific atoms O6/N2 dropped (rebuilt as N6 for A)
        assert "O6" not in primary_atoms
        assert "N2" not in primary_atoms


def test_mutate_raises_on_bad_base(dna_model_list, monkeypatch):
    """An unmappable target base raises a KeyError."""
    with tempfile.TemporaryDirectory() as tmpdir:
        mut_fname = Path(golden_data, dna_model_list[0].file_name)
        monkeypatch.chdir(path=tmpdir)
        with pytest.raises(KeyError):
            mutate(mut_fname, "B", 2, "HOH", "B", 37, "DT")


"""Test the base-pair (double) mutation execution"""


def test_is_cross_type():
    """A cross-type mutation changes the base ring type."""
    assert is_cross_type("DG", "DC") is True  # purine -> pyrimidine
    assert is_cross_type("DC", "DG") is True  # pyrimidine -> purine
    assert is_cross_type("DG", "DA") is False  # purine -> purine
    assert is_cross_type("DC", "DT") is False  # pyrimidine -> pyrimidine


# One-pass and two-pass wild-type baselines used by the run tests. They are
# deliberately different so that the tests can verify which baseline is used.
_NATIVE_1STEP = (-100.0, -10.0, -250.0, -5.0, 1000.0)
_NATIVE_2STEP = (-105.0, -11.0, -255.0, -6.0, 1050.0)


def _make_bp_job(ori, target, p_ori, p_target, output_mutants=False):
    """Helper to build a ModelBasePairMutation with fixed native scores."""
    return ModelBasePairMutation(
        model_path=Path("model.pdb"),
        model_id="m",
        chain="B",
        resid=2,
        ori_resname=ori,
        target_resname=target,
        partner_chain="B",
        partner_resid=37,
        partner_ori_resname=p_ori,
        partner_target_resname=p_target,
        native_scores=_NATIVE_1STEP,
        native_scores_2step=_NATIVE_2STEP,
        output_mutants=output_mutants,
    )


def test_bp_mutation_same_type_single_cns_call(mocker, monkeypatch):
    """A same-ring-type base-pair mutation uses a single CNS call."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)

        def fake_mutate(pdb_f, tc, tr, tt, pc, pr, pt):
            p = Path("mut.pdb")
            p.write_text("ATOM\n")
            return p

        mocker.patch(
            "haddock.modules.analysis.dnascan.dnascan.mutate", side_effect=fake_mutate
        )
        mock_calc = mocker.patch(
            "haddock.modules.analysis.dnascan.dnascan.calc_score",
            return_value=(-90.0, -9.0, -200.0, -4.0, 900.0),
        )
        # DG->DA (purine->purine), partner DC->DT (pyrimidine->pyrimidine)
        result = _make_bp_job("DG", "DA", "DC", "DT").run()
        assert result.success
        assert mock_calc.call_count == 1
        assert result.mutant_scores == (-90.0, -9.0, -200.0, -4.0, 900.0)
        # same-type mutants use the one-pass wild-type baseline
        assert result.delta_scores[0] == pytest.approx(_NATIVE_1STEP[0] - (-90.0))


def test_bp_mutation_cross_type_two_cns_calls(mocker, monkeypatch):
    """A cross-type base-pair mutation uses two sequential CNS calls."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)

        applied = []

        def fake_mutate_residues(pdb_f, mutations):
            applied.append(dict(mutations))
            (chain, resid), target = next(iter(mutations.items()))
            p = Path(f"step_{chain}{resid}{target}.pdb")
            p.write_text("ATOM\n")
            return p

        mocker.patch(
            "haddock.modules.analysis.dnascan.dnascan._mutate_residues",
            side_effect=fake_mutate_residues,
        )

        step_scores = iter(
            [
                (-50.0, -5.0, -100.0, -2.0, 500.0),  # step 1 (discarded)
                (-90.0, -9.0, -200.0, -4.0, 900.0),  # step 2 (kept)
            ]
        )

        def fake_calc(pdb_f, run_dir=None, outputpdb=False, **kwargs):
            if outputpdb:
                Path(f"{Path(pdb_f).stem}_hs.pdb").write_text("ATOM\n")
            return next(step_scores)

        mock_calc = mocker.patch(
            "haddock.modules.analysis.dnascan.dnascan.calc_score",
            side_effect=fake_calc,
        )
        # DG->DC (purine->pyrimidine), partner DC->DG (pyrimidine->purine)
        result = _make_bp_job("DG", "DC", "DC", "DG").run()

        assert result.success
        assert mock_calc.call_count == 2
        # the scores/energies of the second (final) CNS call are kept
        assert result.mutant_scores == (-90.0, -9.0, -200.0, -4.0, 900.0)
        # cross-type mutants are compared to the two-pass wild-type baseline
        assert result.delta_scores[0] == pytest.approx(_NATIVE_2STEP[0] - (-90.0))
        # step 1 mutates the purine->pyrimidine nucleotide (primary DG->DC)
        assert applied[0] == {("B", 2): "DC"}
        # step 2 mutates the pyrimidine->purine nucleotide (partner DC->DG)
        assert applied[1] == {("B", 37): "DG"}
        # the first CNS call regularises the intermediate (writes minimized pdb)
        assert mock_calc.call_args_list[0].kwargs["outputpdb"] is True


def test_bp_mutation_cross_type_step_order_partner_purine(mocker, monkeypatch):
    """When the partner is the purine->pyrimidine one, it is mutated first."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)

        applied = []

        def fake_mutate_residues(pdb_f, mutations):
            applied.append(dict(mutations))
            (chain, resid), target = next(iter(mutations.items()))
            p = Path(f"step_{chain}{resid}{target}.pdb")
            p.write_text("ATOM\n")
            return p

        mocker.patch(
            "haddock.modules.analysis.dnascan.dnascan._mutate_residues",
            side_effect=fake_mutate_residues,
        )

        def fake_calc(pdb_f, run_dir=None, outputpdb=False, **kwargs):
            if outputpdb:
                Path(f"{Path(pdb_f).stem}_hs.pdb").write_text("ATOM\n")
            return (-90.0, -9.0, -200.0, -4.0, 900.0)

        mocker.patch(
            "haddock.modules.analysis.dnascan.dnascan.calc_score", side_effect=fake_calc
        )
        # primary DC->DG (pyrimidine->purine), partner DG->DC (purine->pyrimidine)
        result = _make_bp_job("DC", "DG", "DG", "DC").run()

        assert result.success
        # step 1 is the purine->pyrimidine one, which here is the partner (B37)
        assert applied[0] == {("B", 37): "DC"}
        # step 2 is the pyrimidine->purine primary (B2)
        assert applied[1] == {("B", 2): "DG"}


def test_compute_native_baselines_two_passes(
    mocker, dna_model_list, params, monkeypatch
):
    """The wild-type baseline is scored with one and with two CNS passes."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        scanner = InterfaceScanner(model=dna_model_list[0], params=params)

        scores = iter([(-100.0, -10.0, -250.0, -5.0, 1000.0), _NATIVE_2STEP])

        def fake_calc(pdb_f, run_dir=None, outputpdb=False, **kwargs):
            if outputpdb:
                Path(f"{Path(pdb_f).stem}_hs.pdb").write_text("ATOM\n")
            return next(scores)

        mock_calc = mocker.patch(
            "haddock.modules.analysis.dnascan.dnascan.calc_score", side_effect=fake_calc
        )
        native, native_2step = scanner._compute_native_baselines()

        assert mock_calc.call_count == 2
        assert native == (-100.0, -10.0, -250.0, -5.0, 1000.0)
        assert native_2step == _NATIVE_2STEP
        # the first pass writes the energy-minimized WT, the second re-scores it
        assert mock_calc.call_args_list[0].kwargs["outputpdb"] is True
        assert mock_calc.call_args_list[1].kwargs["outputpdb"] is False
        # the intermediate minimized wild-type PDB is cleaned up
        assert not Path(f"{Path(scanner.model_path).stem}_hs.pdb").exists()


def test_write_delta_score_to_pdb(dna_model_list, results_by_model, monkeypatch):
    """Both nucleotides of a base pair receive the (mean) delta score."""
    with tempfile.TemporaryDirectory() as tmpdir:
        mut_fname = Path(golden_data, dna_model_list[0].file_name)
        monkeypatch.chdir(path=tmpdir)
        mut_pdb_fname = mutate(mut_fname, "B", 2, "DA", "B", 37, "DT")
        dna_model_list[0].file_name = mut_pdb_fname
        dna_model_list[0].path = tmpdir
        add_delta_to_bfactor_obj = AddDeltaBFactor(
            dna_model_list[0],
            tmpdir,
            results_by_model["protdna_complex_1"],
        )
        add_delta_to_bfactor_obj.reorder_results()
        expected = np.mean([5.401, 4.401])
        # both members of the base pair (B2 and its partner B37) get the score
        assert np.isclose(add_delta_to_bfactor_obj.model_results["B-2"], expected)
        assert np.isclose(add_delta_to_bfactor_obj.model_results["B-37"], expected)


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
    scores = calc_score(Path(golden_data, "protdna_complex_1.pdb"), run_dir="tmp")
    assert scores == (-106.7376, -29.5808, -316.542, -13.8484, 1494.73)


def test_dnascan_group_scan_by_cluster(dna_input_list, results_by_model, monkeypatch):
    """Group scan data per cluster, keyed by mutation, with partner metadata."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        clt_scan, clt_pops = group_scan_by_cluster(dna_input_list, results_by_model)
        assert "unclustered" in clt_scan.keys()
        assert clt_pops["unclustered"] == 2
        # The two mutations of base pair B2 must be reported separately
        idents = set(clt_scan["unclustered"].keys())
        assert "B-2-DG-DA" in idents
        assert "B-2-DG-DC" in idents
        # partner metadata is carried along
        entry = clt_scan["unclustered"]["B-2-DG-DA"]
        assert entry["partner_chain"] == "B"
        assert entry["partner_resid"] == 37
        assert entry["partner_target_resname"] == "DT"


def test_dnascan_cluster_full_outputs(dna_input_list, results_by_model, monkeypatch):
    """Test ClusterOutputer writes tsv and plot."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        clt_scan, clt_pops = group_scan_by_cluster(dna_input_list, results_by_model)
        ClusterOutputer(
            clt_scan["unclustered"],
            "unclustered",
            clt_pops["unclustered"],
            scan_residue="DNA base pair",
            generate_plot=True,
            offline=False,
        ).run()
        tsv = Path("scan_clt_unclustered.tsv")
        assert tsv.exists()
        assert Path("scan_clt_unclustered.html").exists()
        # partner columns must be present in the cluster tsv
        df = pd.read_csv(tsv, sep="\t", comment="#")
        assert "partner_target_resname" in df.columns


def test_write_scan_out_with_mutation_results(
    successful_mutation_result, successful_mutation_result_2, monkeypatch
):
    """Test write_scan_out with base-pair (double) mutation results."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        results = [successful_mutation_result, successful_mutation_result_2]
        write_scan_out(results, "protdna_complex_1")
        output_file = Path("scan_protdna_complex_1.tsv")
        assert output_file.exists()
        content = output_file.read_text()
        assert "# `dnascan` results for protdna_complex_1" in content
        assert "delta_score" in content
        assert "z_score" in content
        df = pd.read_csv(output_file, sep="\t", comment="#")
        assert len(df) == 2
        assert df["ori_resname"].tolist() == ["DG", "DG"]
        assert set(df["end_resname"].tolist()) == {"DA", "DC"}
        # partner (double mutation) columns are reported
        assert df["partner_res"].tolist() == [37, 37]
        assert set(df["partner_end_resname"].tolist()) == {"DT", "DG"}
