"""Test the rmsdmatrix module."""
import o
import os
import tempfile
from pathlib import Path

import numpy as np
import pytest


from haddock.modules.analysis.ilrmsdmatrix import DEFAULT_CONFIG as ilrmsd_pars
from haddock.modules.analysis.ilrmsdmatrix import (
    HaddockModule as IlrmsdmatrixModule,
)
from haddock.modules.analysis.ilrmsdmatrix.ilrmsd import (
    ContactJob,
    Contact,
)

from . import golden_data
from .test_module_caprieval import protprot_input_list, protprot_onechain_list  # noqa : F401

from haddock.modules.analysis.ilrmsdmatrix import \
    DEFAULT_CONFIG as ILRMSD_DEFAULT_PARAMS
from haddock.modules.analysis.ilrmsdmatrix import \
    HaddockModule as IlrmsdmatrixModule
from haddock.modules.analysis.ilrmsdmatrix.ilrmsd import Contact, ContactJob



@pytest.fixture(name="params")
def fixture_params():
    """???"""
    return {"receptor_chain": "A", "ligand_chains": ["B"]}


@pytest.fixture
def contact_obj(protprot_input_list, params):  # noqa : F811    
    """Return example Contact object."""
    contact_obj = Contact(
        model_list=protprot_input_list,
        output_name="contact",
        path=Path("."),
        core=0,
        contact_distance_cutoff=5.0,
        params=params,
    )

def ilrmsdmatrix():
    """Return ilrmsdmatrix module."""
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        yield IlrmsdmatrixModule(
            order=1,
            path=Path("."),
            initial_params=ILRMSD_DEFAULT_PARAMS,
        )



@pytest.fixture(name="contact_obj")
def fixture_contact_obj(protprot_input_list, params):
    """Return example Contact object."""
    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)
        yield Contact(
            model_list=protprot_input_list,
            output_name="contact",
            path=Path("."),
            core=0,
            contact_distance_cutoff=5.0,
            params=params,
        )


@pytest.fixture(name="contact_job_obj")
def fixture_contact_job_obj(contact_obj, params):
    """Return example ContactJob object."""
    contact_job_obj = ContactJob(
        Path("contact_output"),
        params,
        contact_obj,
    )
    yield contact_job_obj
    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)
        yield ContactJob(
            Path("contact_output"),
            params,
            contact_obj,
        )


def test_contact(contact_obj):
    """Test the contact class."""
    contact_obj.run()
    contact_obj.output()
    assert Path(contact_obj.output_name).exists()
    exp_lig_res = np.array([10, 11, 12, 16, 17, 48, 51, 52, 53, 54, 56, 57])
    # assert that they are all equal
    assert np.array_equal(contact_obj.unique_lig_res, exp_lig_res)
    # assert contact_obj.unique_lig_res == exp_lig_res
    exp_rec_res = np.array(
        [37, 38, 39, 40, 43, 44, 45, 69, 71, 72, 75, 90, 93, 94, 96, 132]
    )
    assert np.array_equal(contact_obj.unique_rec_res, exp_rec_res)


def test_contact_cutoff(contact_obj):
    """Test the contact class with reduced contact cutoff."""
    # Modify cutoff of the object
    contact_obj.contact_distance_cutoff = 3.9
    contact_obj.run()
    contact_obj.output()
    assert Path(contact_obj.output_name).exists()
    exp_lig_res = np.array([11, 12, 16, 17, 48, 51, 52, 56, 57])
    # assert that they are all equal
    assert np.array_equal(contact_obj.unique_lig_res, exp_lig_res)
    # assert contact_obj.unique_lig_res == exp_lig_res
    exp_rec_res = np.array([38, 39, 40, 45, 69, 71, 72, 90, 94, 96, 132])
    assert np.array_equal(contact_obj.unique_rec_res, exp_rec_res)


def test_contact_job(contact_job_obj):
    """Test ContactJob class."""
    contact_job_obj.run()
    exp_output = contact_job_obj.contact_obj.output_name
    assert Path(exp_output).exists()


@pytest.fixture
def ilrmsdmatrix():
    """Return ilrmsdmatrix module."""
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        yield IlrmsdmatrixModule(
            order=1,
            path=Path(tmpdir),
            initial_params=ilrmsd_pars,
        )


def test_ilrmsdmatrix_init(ilrmsdmatrix):
    """Test ilrmsdmatrix module initialization."""
    ilrmsdmatrix.__init__(
        order=42,
        path=Path("0_anything"),
        initial_params=ilrmsd_pars,
        initial_params=ILRMSD_DEFAULT_PARAMS,
    )
    # Once a module is initialized, it should have the following attributes
    assert ilrmsdmatrix.path == Path("0_anything")
    assert ilrmsdmatrix._origignal_config_file == ILRMSD_DEFAULT_PARAMS
    assert type(ilrmsdmatrix.params) == dict
    assert len(ilrmsdmatrix.params) != 0

@pytest.mark.skip(reason="BEING REFACTORED")
def test_ilrmsdmatrix_run(ilrmsdmatrix, mocker):
    """Test ilrmsdmatrix run method."""
    ilrmsdmatrix.previous_io = MockPreviousIO(path=ilrmsdmatrix.path)
    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models",
        return_value=None,
    )
    ilrmsdmatrix.run()
    assert Path(ilrmsdmatrix.path, "ilrmsd.matrix").exists()
    assert Path(ilrmsdmatrix.path, "receptor_contacts.con").exists()
    with open(Path(ilrmsdmatrix.path, "ilrmsd.matrix")) as f:
        assert f.readline() == f"1 2 16.715{os.linesep}"
    with open(Path(ilrmsdmatrix.path, "receptor_contacts.con")) as f:
        lines = f.readlines()
        assert (
            lines[0]
            == f"A 37 38 39 40 43 44 45 69 71 72 75 90 93 94 96 132{os.linesep}"
        )  # noqa : E501
        assert lines[1] == f"B 10 11 12 16 17 48 51 52 53 54 56 57{os.linesep}"


class MockPreviousIO:
    """Mock previous IO class."""

    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = False):
        """Retrieve models."""
        shutil.copy(
            Path(golden_data, "protprot_complex_1.pdb"),
            Path(".", "protprot_complex_1.pdb"),
        )
        shutil.copy(
            Path(golden_data, "protprot_complex_2.pdb"),
            Path(".", "protprot_complex_2.pdb"),
        )
        model_list = [
            PDBFile(file_name="protprot_complex_1.pdb", path="."),
            PDBFile(file_name="protprot_complex_2.pdb", path="."),
        ]
        return model_list
