"""Integration-test of the CONTact MAP module."""
import os
import pytest
import pytest_mock  # noqa : F401
import tempfile
from pathlib import Path

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.contactmap import HaddockModule as CMapModule
from haddock.modules.analysis.contactmap import DEFAULT_CONFIG as CONTMAP_CONF

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


class MockPreviousIO():
    """A mocking class holding the retrieve_models method."""

    def __init__(self, path):
        self.path = path

    def retrieve_models(individualize: bool = False):
        """Provide a set of models."""
        models = [
            PDBFile(Path(GOLDEN_DATA, "contactmap_rigidbody_3_cltid_None.pdb"),
                    path=GOLDEN_DATA),
            PDBFile(Path(GOLDEN_DATA, "contactmap_rigidbody_5_clt_1.pdb"),
                    path=GOLDEN_DATA),
            PDBFile(Path(GOLDEN_DATA, "contactmap_rigidbody_7_clt_1.pdb"),
                    path=GOLDEN_DATA),
            ]
        # set models cluster ids
        models[0].clt_id = None
        models[1].clt_id = 1
        models[2].clt_id = 1
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
    clustNone_tsv_fpath = f'{output_bp}/Unclustered_contmap_contactmap_rigidbody_3_cltid_None_contacts.tsv'  # noqa : E501
    clustNone_html_fpath = f'{output_bp}/Unclustered_contmap_contactmap_rigidbody_3_cltid_None_heatmap.html'  # noqa : E501
    assert os.path.exists(clustNone_tsv_fpath)
    assert Path(clustNone_tsv_fpath).stat().st_size != 0
    assert os.path.exists(clustNone_html_fpath)
    assert Path(clustNone_html_fpath).stat().st_size != 0
    Path(clustNone_tsv_fpath).unlink(missing_ok=False)
    Path(clustNone_html_fpath).unlink(missing_ok=False)

    # clt_id == 1
    clust1_tsv_fpath = f'{output_bp}/cluster1_contmap_contacts.tsv'
    clust1_html_fpath = f'{output_bp}/cluster1_contmap_heatmap.html'
    assert os.path.exists(clust1_tsv_fpath)
    assert Path(clust1_tsv_fpath).stat().st_size != 0
    assert os.path.exists(clust1_html_fpath)
    assert Path(clust1_html_fpath).stat().st_size != 0
    Path(clust1_tsv_fpath).unlink(missing_ok=False)
    Path(clust1_html_fpath).unlink(missing_ok=False)
