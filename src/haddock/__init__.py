"""HADDOCK3 library."""
from pathlib import Path


haddock3_source_path = Path(__file__).resolve().parent
haddock3_repository_path = haddock3_source_path.parents[1]
toppar_path = Path(haddock3_source_path, "cns", "toppar")


FCC_path = Path(
    Path(__file__).resolve().parents[1],
    'fcc',
    )

v_major = "3"
v_minor = "0"
v_patch = "beta"
v_release = "unreleased"

current_version = f"{v_major}.{v_minor}.{v_patch}-{v_release}"
contact_us = 'https://github.com/haddocking/haddock3/issues'
