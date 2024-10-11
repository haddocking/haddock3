import shutil
import tempfile
import gzip
from pathlib import Path

import pytest

from haddock.libs.libontology import Format, PDBFile, Persistent
from haddock.modules.sampling.rigidbody import (
    DEFAULT_CONFIG as DEFAULT_RIGIDBODY_CONFIG,
)
from haddock.modules.sampling.rigidbody import HaddockModule as RigidbodyModule

from tests import golden_data as testgolden_data
from integration_tests import GOLDEN_DATA


@pytest.fixture
def rigidbody_module():
    with tempfile.TemporaryDirectory() as tmpdir:
        rigidbody = RigidbodyModule(
            order=0, path=Path(tmpdir), initial_params=DEFAULT_RIGIDBODY_CONFIG
        )
        yield rigidbody


class MockPreviousIO:
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, crossdock: bool = False):
        shutil.copy(
            Path(testgolden_data, "e2aP_1F3G_haddock.pdb"),
            Path(".", "e2aP_1F3G_haddock.pdb"),
        )
        shutil.copy(
            Path(testgolden_data, "e2aP_1F3G_haddock.psf"),
            Path(self.path, "e2aP_1F3G_haddock.psf"),
        )
        shutil.copy(
            Path(testgolden_data, "hpr_ensemble_1_haddock.pdb"),
            Path(self.path, "hpr_ensemble_1_haddock.pdb"),
        )
        shutil.copy(
            Path(testgolden_data, "hpr_ensemble_1_haddock.psf"),
            Path(self.path, "hpr_ensemble_1_haddock.psf"),
        )
        model_list = [
            [
                PDBFile(
                    file_name="e2aP_1F3G_haddock.pdb",
                    path=self.path,
                    topology=[
                        Persistent(
                            file_name="e2aP_1F3G_haddock.psf",
                            path=self.path,
                            file_type=Format.TOPOLOGY,
                        )
                    ],
                ),
                PDBFile(
                    file_name="hpr_ensemble_1_haddock.pdb",
                    path=self.path,
                    topology=[
                        Persistent(
                            file_name="hpr_ensemble_1_haddock.psf",
                            path=self.path,
                            file_type=Format.TOPOLOGY,
                        )
                    ],
                ),
            ]
        ]

        return model_list

    def output(self):
        return None


def test_restraints_rigidbody(rigidbody_module):

    rigidbody_module.previous_io = MockPreviousIO(path=rigidbody_module.path)
    rigidbody_module.params["sampling"] = 1
    rigidbody_module.params["mol_fix_origin_1"] = False
    rigidbody_module.params["mol_fix_origin_2"] = False
    rigidbody_module.params["ambig_fname"] = Path(GOLDEN_DATA, "ambig.tbl")
    rigidbody_module.params["unambig_fname"] = Path(GOLDEN_DATA, "unambig.tbl")
    rigidbody_module.params["hbond_fname"] = Path(GOLDEN_DATA, "hbond.tbl")
    rigidbody_module.params["mode"] = "local"
    rigidbody_module.params["debug"] = True

    rigidbody_module.run()

    assert Path(rigidbody_module.path, "rigidbody_1.pdb").exists()
    assert Path(rigidbody_module.path, "rigidbody_1.out.gz").exists()
    assert Path(rigidbody_module.path, "rigidbody_1.inp").exists()

    assert Path(rigidbody_module.path, "rigidbody_1.pdb").stat().st_size > 0
    assert Path(rigidbody_module.path, "rigidbody_1.out.gz").stat().st_size > 0
    assert Path(rigidbody_module.path, "rigidbody_1.inp").stat().st_size > 0

    file_content = gzip.open(Path(rigidbody_module.path, "rigidbody_1.out.gz"), 'rt').read()
    assert 'NOEPRI: RMS diff. class AMBI' in file_content
    assert 'NOEPRI: RMS diff. class DIST' in file_content
    assert 'NOEPRI: RMS diff. class HBON' in file_content

