"""Test the rmsdmatrix module."""
import os
from pathlib import Path

import numpy as np
import pytest
import pytest_mock
import shutil
import tempfile


from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.ilrmsdmatrix import DEFAULT_CONFIG as ilrmsd_pars
from haddock.modules.analysis.ilrmsdmatrix import HaddockModule as IlrmsdmatrixModule
from haddock.modules.analysis.ilrmsdmatrix import get_index_list
from haddock.modules.analysis.ilrmsdmatrix.ilrmsd import (
    ContactJob,
    Contact,
    ilRMSDJob,
    ilRMSD,
    get_pair,
   )

from . import golden_data

from .test_module_caprieval import protprot_input_list, protprot_onechain_list

@pytest.fixture
def params():
    return {"receptor_chain":"A", "ligand_chains":["B"]}

@pytest.fixture
def contact_obj(protprot_input_list, params):
    """Return example alascan module."""
    contact_obj = Contact(
        model_list=protprot_input_list,
        output_name="contact",
        path=Path("."),
        core=0,
        params=params,
    )

    yield contact_obj

@pytest.fixture
def contact_job_obj(contact_obj, params):
    """Return example alascan module."""
    contact_job_obj = ContactJob(
        Path("contact_output"),
        params,
        contact_obj,
        )
    yield contact_job_obj

@pytest.fixture
def ilrmsd_obj(protprot_input_list, params):
    """Return example ilRMSD object."""
    ilrmsd_obj = ilRMSD(
        model_list=protprot_input_list,
        output_name="ilrmsd",
        path=Path("."),
        core=0,
        params=params,
        start_ref=0,
        end_ref=1,
        filter_resdic=None,
        start_mod=1,
        npairs=1,
    )

    yield ilrmsd_obj

@pytest.fixture
def ilrmsdjob_obj(ilrmsd_obj, params):
    """Return example ilRMSDJob object."""
    ilrmsdjob_obj = ilRMSDJob(
        Path("ilrmsd_output"),
        params,
        ilrmsd_obj,
        )
    yield ilrmsdjob_obj

@pytest.fixture
def ilrmsd_obj_onechain(protprot_onechain_list, params):
    """Return example ilRMSD object."""
    ilrmsd_obj_onechain = ilRMSD(
        model_list=protprot_onechain_list + protprot_onechain_list,
        output_name="ilrmsd",
        path=Path("."),
        core=0,
        params=params,
        start_ref=0,
        filter_resdic=None,
        start_mod=1,
        npairs=1,
    )

    yield ilrmsd_obj_onechain

def test_contact(contact_obj):
    """Test the contact class."""
    contact_obj.run()
    contact_obj.output()
    assert Path(contact_obj.output_name).exists()
    exp_lig_res = np.array([10, 11, 12, 16, 17, 48, 51, 52, 53, 54, 56, 57])
    # assert that they are all equal
    assert np.array_equal(contact_obj.unique_lig_res, exp_lig_res)
    #assert contact_obj.unique_lig_res == exp_lig_res
    exp_rec_res = np.array(
        [ 37,  38,  39,  40,  43,  44,  45,  69,  71,  72,  75,  90,  93, 94,  96, 132])
    assert np.array_equal(contact_obj.unique_rec_res, exp_rec_res)


def test_contact_job(contact_job_obj, contact_obj):
    """Test ContactJob class."""
    contact_job_obj.run()
    exp_output = contact_job_obj.contact_obj.output_name
    assert Path(exp_output).exists()


def test_ilrmsd(ilrmsd_obj):
    """Test the ilRMSD calc."""
    assert ilrmsd_obj.atoms != []
    ilrmsd_obj.run()
    exp_data = np.array([[1., 2., 20.93705795]])
    assert np.allclose(ilrmsd_obj.data, exp_data)
    ilrmsd_obj.output()
    assert Path(ilrmsd_obj.output_name).exists()
    # remove stuff
    os.unlink(ilrmsd_obj.output_name)
    

def test_ilrmsd_obj_onechain(ilrmsd_obj_onechain):
    """Test the ilRMSD class."""
    ilrmsd_obj_onechain.run()
    # all zeros
    assert ilrmsd_obj_onechain.data[0][0] == 0.0
    assert ilrmsd_obj_onechain.data[0][1] == 0.0
    assert ilrmsd_obj_onechain.data[0][2] == 0.0


def test_ilrmsdjob(ilrmsdjob_obj):
    """Test the ilRMSDJob class."""
    ilrmsdjob_obj.run()
    exp_output = ilrmsdjob_obj.ilrmsd_obj.output_name
    assert Path(exp_output).exists()
    # remove stuff
    os.unlink(ilrmsdjob_obj.ilrmsd_obj.output_name)


def test_get_pair():
    """Test get_pair function."""
    obs_pair = get_pair(10, 10)
    assert obs_pair == (1, 3)
    with pytest.raises(ValueError):
        get_pair(-10, 11)


def test_get_index_list():
    assert True

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
    ilrmsdmatrix.__init__(
        order=42,
        path=Path("0_anything"),
        initial_params=ilrmsd_pars,
    )

    # Once a module is initialized, it should have the following attributes
    assert ilrmsdmatrix.path == Path("0_anything")
    assert ilrmsdmatrix._origignal_config_file == ilrmsd_pars
    assert type(alascan.params) == dict
    assert len(alascan.params) != 0

def test_ilrmsdmatrix_run(ilrmsdmatrix, mocker, contact_obj):
    ilrmsdmatrix.previous_io = MockPreviousIO()
    mocker.patch("haddock.modules.BaseHaddockModule.export_io_models", return_value = None)
    mocker.patch("haddock.libs.libparallel.Scheduler.run", return_value = None)
    mocker.patch("haddock.libs.libparallel.Scheduler", return_value = None)
    ilrmsdmatrix.run()

class MockPreviousIO:
    def retrieve_models(self, individualize: bool = False):
        shutil.copy(Path("..", "tests", "golden_data", "protprot_complex_1.pdb"), Path(".", "protprot_complex_1.pdb"))
        shutil.copy(Path("..", "tests", "golden_data", "protprot_complex_2.pdb"), Path(".", "protprot_complex_2.pdb"))
        model_list = [
            PDBFile(file_name="protprot_complex_1.pdb", path="."),
            PDBFile(file_name="protprot_complex_2.pdb", path="."),
        ]
        return model_list

    
