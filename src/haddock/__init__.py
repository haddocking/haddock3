"""HADDOCK3 library."""
from pathlib import Path


haddock3_source_path = Path(__file__).resolve().parent
haddock3_repository_path = haddock3_source_path.parents[1]
toppar_path = Path(haddock3_source_path, "cns", "toppar")


FCC_path = Path(
    Path(__file__).resolve().parents[1],
    'fcc',
    )


# version
version = "3.0.0"
v_major, v_minor, v_patch = current_version.split('.')

contact_us = 'https://github.com/haddocking/haddock3/issues'
