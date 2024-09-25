import shutil
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import Format, PDBFile, Persistent
from haddock.modules.sampling.rigidbody import \
    DEFAULT_CONFIG as DEFAULT_RIGIDBODY_CONFIG
from haddock.modules.sampling.rigidbody import HaddockModule as RigidbodyModule
from tests import golden_data


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
            Path(golden_data, "e2aP_1F3G_haddock.pdb"),
            Path(self.path, "e2aP_1F3G_haddock.pdb"),
        )
        shutil.copy(
            Path(golden_data, "e2aP_1F3G_haddock.psf"),
            Path(self.path, "e2aP_1F3G_haddock.psf"),
        )
        shutil.copy(
            Path(golden_data, "hpr_ensemble_1_haddock.pdb"),
            Path(self.path, "hpr_ensemble_1_haddock.pdb"),
        )
        shutil.copy(
            Path(golden_data, "hpr_ensemble_1_haddock.psf"),
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
                            path=".",
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
                            path=".",
                            file_type=Format.TOPOLOGY,
                        )
                    ],
                ),
            ]
        ]

        return model_list

    def output(self):
        return None


def test_rigidbody_local(rigidbody_module):

    sampling = 2
    rigidbody_module.previous_io = MockPreviousIO(path=rigidbody_module.path)
    rigidbody_module.params["sampling"] = sampling
    rigidbody_module.params["cmrest"] = True
    rigidbody_module.params["mol_fix_origin_1"] = True
    rigidbody_module.params["mol_fix_origin_2"] = False
    rigidbody_module.params["mode"] = "local"
    rigidbody_module.params["debug"] = True

    rigidbody_module.run()

    for i in range(1, sampling + 1):
        assert Path(rigidbody_module.path, f"rigidbody_{i}.pdb").exists()
        assert Path(rigidbody_module.path, f"rigidbody_{i}.out.gz").exists()
        assert Path(rigidbody_module.path, f"rigidbody_{i}.inp").exists()
        assert not Path(rigidbody_module.path, f"rigidbody_{i}.seed").exists()

        assert Path(rigidbody_module.path, f"rigidbody_{i}.pdb").stat().st_size > 0
        assert Path(rigidbody_module.path, f"rigidbody_{i}.out.gz").stat().st_size > 0
        assert Path(rigidbody_module.path, f"rigidbody_{i}.inp").stat().st_size > 0


def test_rigidbody_mpi(rigidbody_module):

    sampling = 2
    rigidbody_module.previous_io = MockPreviousIO(path=rigidbody_module.path)
    rigidbody_module.params["sampling"] = sampling
    rigidbody_module.params["ntrials"] = 1
    rigidbody_module.params["cmrest"] = True
    rigidbody_module.params["mol_fix_origin_1"] = True
    rigidbody_module.params["mol_fix_origin_2"] = False
    rigidbody_module.params["mode"] = "mpi"
    rigidbody_module.params["ncores"] = 1
    rigidbody_module.params["debug"] = True

    rigidbody_module.run()

    for i in range(1, sampling + 1):
        assert Path(rigidbody_module.path, f"rigidbody_{i}.pdb").exists()
        assert Path(rigidbody_module.path, f"rigidbody_{i}.out.gz").exists()
        assert Path(rigidbody_module.path, f"rigidbody_{i}.inp").exists()
        assert not Path(rigidbody_module.path, f"rigidbody_{i}.seed").exists()

        assert Path(rigidbody_module.path, f"rigidbody_{i}.pdb").stat().st_size > 0
        assert Path(rigidbody_module.path, f"rigidbody_{i}.out.gz").stat().st_size > 0
        assert Path(rigidbody_module.path, f"rigidbody_{i}.inp").stat().st_size > 0


def test_rigidbody_debug(rigidbody_module):

    sampling = 2
    rigidbody_module.previous_io = MockPreviousIO(path=rigidbody_module.path)
    rigidbody_module.params["sampling"] = sampling
    rigidbody_module.params["cmrest"] = True
    rigidbody_module.params["mol_fix_origin_1"] = True
    rigidbody_module.params["mol_fix_origin_2"] = False

    rigidbody_module.run()

    for i in range(1, sampling + 1):
        assert Path(rigidbody_module.path, f"rigidbody_{i}.pdb").exists()
        assert not Path(rigidbody_module.path, f"rigidbody_{i}.out.gz").exists()
        assert not Path(rigidbody_module.path, f"rigidbody_{i}.inp").exists()
        assert not Path(rigidbody_module.path, f"rigidbody_{i}.seed").exists()

        assert Path(rigidbody_module.path, f"rigidbody_{i}.pdb").stat().st_size > 0
