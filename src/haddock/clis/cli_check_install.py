"""Check if HADDOCK3 was installed correctly."""

import importlib.resources
import logging
import platform
import sys
from pathlib import Path

import haddock


logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s - %(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


BINARY_DIR = Path(importlib.resources.files(haddock) / "bin")  # type: ignore

BINARIES = {
    "CNS": "cns",
    "FCC": "contact_fcc",
    "Fast-RMSDMatrix": "fast-rmsdmatrix",
}


def check_binary(software: str, binary_name: str) -> None:
    """Check if the binary is installed."""
    expected_binary_location = Path(BINARY_DIR, binary_name)
    if expected_binary_location.exists():
        logging.info(f"✅ {software} binary `{binary_name}` OK")
    else:
        logging.info(
            f"⚠️ {software} binary `{binary_name}` not properly installed "
            "- define it manually in the configuration file"
        )


def main():
    """Check if the installation was successful."""
    system = platform.system().lower()
    machine = platform.machine().lower()
    python_version = f"{sys.version_info.major}.{sys.version_info.minor}"

    logging.info("Checking HADDOCK3 installation...")
    logging.info(f"ℹ️ Python{python_version} running on {machine}-{system}")
    for software, binary_name in BINARIES.items():
        check_binary(software, binary_name)


if __name__ == "__main__":
    main()
