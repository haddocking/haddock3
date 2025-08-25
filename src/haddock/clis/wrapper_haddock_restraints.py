#!/usr/bin/env python3
import sys
from pathlib import Path

from importlib.resources import files


def main():

    binary = Path(files("haddock").joinpath("bin/haddock-restraints"))  # type: ignore
    sys.argv[0] = str(binary)  # Replace argv[0] with actual binary path

    # Execute the binary
    import subprocess

    sys.exit(subprocess.run([str(binary)] + sys.argv[1:]).returncode)


if __name__ == "__main__":
    main()
