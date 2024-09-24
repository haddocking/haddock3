"""Test the rmsdmatrix module."""

import tempfile
from pathlib import Path

import numpy as np
import pytest
import pytest_mock  # noqa : F401

from haddock.modules.analysis.ilrmsdmatrix import DEFAULT_CONFIG as ilrmsd_pars
from haddock.modules.analysis.ilrmsdmatrix import \
    HaddockModule as IlrmsdmatrixModule
from haddock.modules.analysis.ilrmsdmatrix.ilrmsd import Contact, ContactJob

from .test_module_caprieval import (  # noqa : F401
    protprot_input_list,
    protprot_onechain_list,
    )


@pytest.fixture
def params():
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

    yield contact_obj


@pytest.fixture
def contact_job_obj(contact_obj, params):
    """Return example ContactJob object."""
    contact_job_obj = ContactJob(
        Path("contact_output"),
        params,
        contact_obj,
    )
    yield contact_job_obj


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
    with tempfile.TemporaryDirectory() as tmpdir:
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
    )
    # Once a module is initialized, it should have the following attributes
    assert ilrmsdmatrix.path == Path("0_anything")
    assert ilrmsdmatrix._origignal_config_file == ilrmsd_pars
    assert type(ilrmsdmatrix.params) == dict
    assert len(ilrmsdmatrix.params) != 0
