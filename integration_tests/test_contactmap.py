"""Integration-test of the CONTact MAP module."""

import os
import tracemalloc
import pytest
import tempfile
from pathlib import Path

from integration_tests.conftest import generate_synthetic_pdb

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.contactmap import HaddockModule as CMapModule
from haddock.modules.analysis.contactmap import DEFAULT_CONFIG as CONTMAP_CONF
from haddock.modules.analysis.contactmap.contmap import ContactsMap

from integration_tests import GOLDEN_DATA


@pytest.fixture
def contactmap():
    """Return contmap module."""
    with tempfile.TemporaryDirectory() as tmpdir:
        preset_contactmap = CMapModule(
            order=0,
            path=Path(tmpdir),
            initial_params=CONTMAP_CONF,
        )
        yield preset_contactmap


class MockPreviousIO:
    """A mocking class holding the retrieve_models method."""

    def __init__(self, path):
        self.path = path

    def retrieve_models(individualize: bool = False):
        """Provide a set of models."""
        models = [
            PDBFile(
                Path(GOLDEN_DATA, "contactmap_rigidbody_3_cltid_None.pdb"),
                path=GOLDEN_DATA,
            ),
            PDBFile(
                Path(GOLDEN_DATA, "contactmap_rigidbody_5_clt_1.pdb"), path=GOLDEN_DATA
            ),
            PDBFile(
                Path(GOLDEN_DATA, "contactmap_rigidbody_7_clt_1.pdb"), path=GOLDEN_DATA
            ),
        ]
        # set models cluster ids
        models[0].clt_id = None
        models[1].clt_id = 1
        models[1].clt_rank = 2
        models[2].clt_id = 1
        models[2].clt_rank = 2
        return models


def test_contactmap_example(contactmap, mocker):
    """Test the contact map module run."""
    # mock the previous_io behavior
    contactmap.previous_io = MockPreviousIO
    # mock the export_io_models function
    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models",
        return_value=None,
    )
    # Run the module
    contactmap.run()
    # check outputs
    output_bp = contactmap.path
    # clt_id == None
    clustNone_tsv_fpath = (
        f"{output_bp}/Unclustered_contactmap_rigidbody_3_cltid_None_contacts.tsv"  # noqa : E501
    )
    clustNone_html_fpath = (
        f"{output_bp}/Unclustered_contactmap_rigidbody_3_cltid_None_heatmap.html"  # noqa : E501
    )
    assert os.path.exists(clustNone_tsv_fpath)
    assert Path(clustNone_tsv_fpath).stat().st_size != 0
    assert os.path.exists(clustNone_html_fpath)
    assert Path(clustNone_html_fpath).stat().st_size != 0
    Path(clustNone_tsv_fpath).unlink(missing_ok=False)
    Path(clustNone_html_fpath).unlink(missing_ok=False)

    # clt_id == 1
    clust1_tsv_fpath = f"{output_bp}/cluster2_contacts.tsv"
    clust1_html_fpath = f"{output_bp}/cluster2_heatmap.html"
    assert os.path.exists(clust1_tsv_fpath)
    assert Path(clust1_tsv_fpath).stat().st_size != 0
    assert os.path.exists(clust1_html_fpath)
    assert Path(clust1_html_fpath).stat().st_size != 0
    Path(clust1_tsv_fpath).unlink(missing_ok=False)
    Path(clust1_html_fpath).unlink(missing_ok=False)


def test_contactmap_memory_scaling(monkeypatch):
    """Verify ContactsMap workflow doesn't scale memory quadratically."""
    sizes = [100, 200]  # Atom counts to test
    memory_usage = []

    # Set up default params
    params = {
        "ca_ca_dist_threshold": 9.0,
        "shortest_dist_threshold": 7.5,
        "color_ramp": "Greys",
        "single_model_analysis": False,
        "topX": 10,
        "generate_heatmap": False,
        "cluster_heatmap_datatype": "shortest-cont-probability",
        "generate_chordchart": False,
        "chordchart_datatype": "shortest-dist",
        "offline": False,
    }

    for n_atoms in sizes:
        with (
            tempfile.TemporaryDirectory() as tempdir,
            tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as pdb_file,
            tempfile.NamedTemporaryFile(suffix="_output", delete=False) as output_file,
        ):
            monkeypatch.chdir(tempdir)

            # Generate synthetic PDB
            pdb_path = Path(pdb_file.name)
            generate_synthetic_pdb(n_atoms, pdb_path)
            print(pdb_path)

            # Track memory during entire ContactsMap run
            tracemalloc.start()

            contactmap = ContactsMap(
                model=pdb_path,
                output=Path(output_file.name),
                params=params,
            )
            res_contacts, heavy_contacts = contactmap.run()

            _, peak_mem = tracemalloc.get_traced_memory()
            tracemalloc.stop()

            # Verify results were generated
            assert isinstance(res_contacts, list)
            assert isinstance(heavy_contacts, list)
            assert len(res_contacts) > 0

            # Record memory usage in megabytes
            peak_mb = peak_mem / (1024**2)
            memory_usage.append(peak_mb)

            # Cleanup
            pdb_path.unlink(missing_ok=True)
            Path(output_file.name).unlink(missing_ok=True)

    # Check that memory usage does NOT increase with atom count
    # NOTE: If implementation is memory-efficient (e.g., chunked processing),
    # memory should stay relatively constant regardless of input size
    baseline_memory = memory_usage[0]

    for size, mem in zip(sizes[1:], memory_usage[1:]):
        # Memory should not increase significantly from baseline
        # Allow 50% increase to account for overhead and variation
        memory_increase = mem - baseline_memory
        percent_increase = (memory_increase / baseline_memory) * 100

        assert memory_increase < baseline_memory * 0.5, (
            f"Memory increased! Baseline ({sizes[0]} atoms): {baseline_memory:.2f}MB → "
            f"Current ({size} atoms): {mem:.2f}MB "
            f"(+{percent_increase:.1f}%). "
            f"Expected memory to stay constant with efficient implementation."
        )
