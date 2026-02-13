import re
import shutil
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import Format, PDBFile, Persistent
from haddock.modules.sampling.rigidbody import (
    DEFAULT_CONFIG as DEFAULT_RIGIDBODY_CONFIG,
)
from haddock.modules.sampling.rigidbody import HaddockModule as RigidbodyModule
from integration_tests import GOLDEN_DATA

@pytest.fixture(name="rigidbody_module")
def fixture_rigidbody_module():
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
            Path(GOLDEN_DATA, "gcn4_A_haddock.pdb"),
            Path(self.path, "gcn4_A_haddock.pdb"),
        )
        shutil.copy(
            Path(GOLDEN_DATA, "gcn4_A_haddock.psf"),
            Path(self.path, "gcn4_A_haddock.psf"),
        )
        shutil.copy(
            Path(GOLDEN_DATA, "gcn4_B_haddock.pdb"),
            Path(self.path, "gcn4_B_haddock.pdb"),
        )
        shutil.copy(
            Path(GOLDEN_DATA, "gcn4_B_haddock.psf"),
            Path(self.path, "gcn4_B_haddock.psf"),
        )
        shutil.copy(
            Path(GOLDEN_DATA, "custom_sym_restraints.tbl"),
            Path(self.path, "custom_sym_restraints.tbl"),
        )
        model_list = [
            [
                PDBFile(
                    file_name="gcn4_A_haddock.pdb",
                    path=self.path,
                    topology=[
                        Persistent(
                            file_name="gcn4_A_haddock.psf",
                            path=".",
                            file_type=Format.TOPOLOGY,
                        )
                    ],
                ),
                PDBFile(
                    file_name="gcn4_B_haddock.pdb",
                    path=self.path,
                    topology=[
                        Persistent(
                            file_name="gcn4_B_haddock.psf",
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
    rigidbody_module.params["sym_on"] = True
    # Copy the symmetry restraints file to the module path before setting the parameter
    rigidbody_module.previous_io.retrieve_models()
    rigidbody_module.params["symtbl_fname"] = Path(rigidbody_module.path, "custom_sym_restraints.tbl")

    rigidbody_module.run()

    for i in range(1, sampling + 1):
        assert Path(rigidbody_module.path, f"rigidbody_{i}.pdb").exists()
        assert Path(rigidbody_module.path, f"rigidbody_{i}.pdb").stat().st_size > 0

        pattern = re.compile(r"REMARK Symmetry energy:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)")
        found = False
        with open(Path(rigidbody_module.path, f"rigidbody_{i}.pdb"), "r") as f:
            for line in f:
                match = pattern.search(line)
                if match:
                   value = float(match.group(1))
                   assert value > 0, f"Symmetry energy not > 0 (found {value})"
                   found = True

        assert found, "No 'REMARK Symmetry energy' line found in file"

