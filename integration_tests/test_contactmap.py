"""Integration-test of the CONTact MAP module."""

import glob
import os
import tempfile
from pathlib import Path

import psutil
import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.contactmap import DEFAULT_CONFIG as CONTMAP_CONF
from haddock.modules.analysis.contactmap import HaddockModule as CMapModule
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


@pytest.mark.skipif(
    psutil.virtual_memory().available < 3000000000,
    reason="not enough memory to run this test",
)
def test_contactmap_example(contactmap, monkeypatch, mocker):
    """Test the contact map module run."""
    # mock the previous_io behavior
    contactmap.previous_io = MockPreviousIO
    # mock the export_io_models function
    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models",
        return_value=None,
    )
    monkeypatch.chdir(contactmap.path)
    # Run the module
    contactmap.run()
    # check outputs
    output_bp = contactmap.path
    # clt_id == None
    clustNone_tsv_fpath = (
        f"{output_bp}/Unclustered_contactmap_rigidbody_3_cltid_None_contacts.tsv"
    )
    clustNone_html_fpath = (
        f"{output_bp}/Unclustered_contactmap_rigidbody_3_cltid_None_heatmap.html"
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


def test_contactmap_low_memory(contactmap, monkeypatch, mocker):
    """Test the contact map module fails gracefully with insufficient memory."""
    contactmap.previous_io = MockPreviousIO
    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models",
        return_value=None,
    )

    mocker.patch(
        "haddock.modules.analysis.contactmap.get_available_memory",
        return_value=0.0,
    )
    mocker.patch(
        "haddock.modules.analysis.contactmap.get_necessary_memory",
        return_value=1.0,
    )
    monkeypatch.chdir(contactmap.path)

    # Run the module - should skip execution due to low memory
    contactmap.run()

    # Check that the directory is empty
    ls = list(glob.glob(f"{contactmap.path}/*"))
    assert len(ls) == 0
