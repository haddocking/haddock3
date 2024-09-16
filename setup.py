#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""Setup dot py."""
import os
import platform
import shutil
import subprocess
import sys
import urllib.request
from os.path import dirname, join
from pathlib import Path

from setuptools import find_packages, setup
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install


class CustomBuildInstall(install):
    """Custom Build and Install class"""

    def run(self):

        install.run(self)


class CustomBuild(build_ext):
    """CustomBuild handles the build of the C/C++ dependencies"""

    def run(self):
        print("Building HADDOCK3 C/C++ binary dependencies...")
        bin_dir = Path(self.get_install_dir(), "haddock", "bin")
        bin_dir.mkdir(exist_ok=True, parents=True)

        dep_dir = Path("src", "haddock", "deps")
        assert dep_dir.exists()

        # Build FCC
        cmd = "g++ -O2 -o contact_fcc contact_fcc.cpp"

        _ = subprocess.run(
            cmd,
            shell=True,
            cwd=dep_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )

        src = Path(dep_dir, "contact_fcc")
        dst = Path(bin_dir, "contact_fcc")

        shutil.move(src, dst)

        # Build fast-rmsdmatrix
        cmd = "gcc -Wall -O3 -march=native -std=c99 -o fast-rmsdmatrix fast-rmsdmatrix.c -lm"
        _ = subprocess.run(
            cmd,
            shell=True,
            cwd=dep_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )

        src = Path(dep_dir, "fast-rmsdmatrix")
        dst = Path(bin_dir, "fast-rmsdmatrix")

        shutil.move(src, dst)

        # Run original build_ext
        build_ext.run(self)

    def get_install_dir(self):
        """Get the directory in which HADDOCK was installed"""
        install_cmd = self.get_finalized_command("install")
        return install_cmd.install_lib  # type: ignore


CNS_BINARIES = {
    "x86_64-linux": "https://surfdrive.surf.nl/files/index.php/s/o8X7zZezIM3P0cE",  # linux
    "x86_64-darwin": "",  # macOs
    "arm64-darwin": "",  # Mac M1/M2
    "aarch64-linux": "",  # linux ARM
    "armv7l-linux": "",  # 32-bit ARM linux, like raspberryPi
}


class CustomInstall(CustomBuildInstall):
    """Custom class to handle the download of the CNS binary"""

    def run(self):
        """Run the installation"""

        # Run the standard installation
        CustomBuildInstall.run(self)

        # Get the installation directory
        if self.install_lib is None:
            print("Something went wrong during installation.")
            sys.exit(1)

        # Set where the cns binary needs to be
        bin_dir = Path(self.install_lib, "haddock", "bin")

        # Create the `bin/` directory
        bin_dir.mkdir(exist_ok=True)

        # Download the binary
        cns_exec = Path(bin_dir, "cns")
        if cns_exec.exists():
            cns_exec.unlink()

        arch = self.get_arch()

        if arch not in CNS_BINARIES:
            print(f"Unknown architecture: {arch}")
            print(
                "Please set the CNS binary manually inside your configuration file as `cns_exec`",
            )
            return

        cns_binary_url = CNS_BINARIES[arch]

        status, msg = self.download_file(cns_binary_url, cns_exec)
        if not status:
            print(msg)
            print(
                "Please set the CNS binary manually inside your configuration file as `cns_exec`"
            )
            return

        os.chmod(cns_exec, 0o755)

    @staticmethod
    def download_file(url, dest) -> tuple[bool, str]:
        """Download a file from a URL"""
        try:
            urllib.request.urlretrieve(url, dest)
            return True, "Download successful"
        except Exception as e:  # pylint: disable=broad-except
            return False, f"Download failed: {e}"

    @staticmethod
    def get_arch():
        """Helper function to figure out the architetchure"""
        system = platform.system().lower()
        machine = platform.machine().lower()

        return f"{machine}-{system}"


with open("requirements.txt", "r", encoding="utf-8") as f:
    requirements = f.read().splitlines()


def read_description(*names, **kwargs) -> str:
    """Read description files."""
    path = join(dirname(__file__), *names)
    with open(path, encoding=kwargs.get("encoding", "utf8")) as fh:
        return fh.read()


readme = read_description("README.md")
changelog = read_description("CHANGELOG.md")
long_description = f"{readme}{os.linesep}{changelog}"


setup(
    name="haddock3",
    version="3.0.0",
    description="HADDOCK3",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="Apache License 2.0",
    author="BonvinLab",
    author_email="bonvinlab.support@uu.nl",
    url="https://github.com/haddocking/haddock3",
    packages=find_packages("src"),
    package_dir={"": "src"},
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # TODO: Update the classifiers - http://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Operating System :: POSIX",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
        "Programming Language :: Python :: 3.9",
    ],
    project_urls={
        "webpage": "https://bonvinlab.org/haddock3",
        "Documentation": "https://github.com/haddocking/haddock3#readme",
        "Changelog": "",
        "Issue Tracker": "https://github.com/haddocking/haddock3/issues",
        "Discussion Forum": "https://github.com/haddocking/haddock3/issues",
    },
    keywords=[
        "Structural Biology",
        "Biochemistry",
        "Docking",
        "Protein docking",
        "Proteins",
    ],
    python_requires=">=3.9, <3.10",
    install_requires=[requirements],
    extras_require={},
    setup_requires=[],
    entry_points={
        "console_scripts": [
            "haddock3 = haddock.clis.cli:maincli",
            "haddock3-mpitask = haddock.clis.cli_mpi:maincli",
            "haddock3-bm = haddock.clis.cli_bm:maincli",
            "haddock3-cfg = haddock.clis.cli_cfg:maincli",
            "haddock3-clean = haddock.clis.cli_clean:maincli",
            "haddock3-copy = haddock.clis.cli_cp:maincli",
            "haddock3-dmn = haddock.clis.cli_dmn:maincli",
            "haddock3-pp = haddock.clis.cli_pp:maincli",
            "haddock3-score = haddock.clis.cli_score:maincli",
            "haddock3-unpack = haddock.clis.cli_unpack:maincli",
            "haddock3-analyse = haddock.clis.cli_analyse:maincli",
            "haddock3-traceback = haddock.clis.cli_traceback:maincli",
            "haddock3-re = haddock.clis.cli_re:maincli",
            "haddock3-restraints = haddock.clis.cli_restraints:maincli",
            "haddock3-check = haddock.clis.cli_check_install:main",
        ]
    },
    cmdclass={"build_ext": CustomBuild, "install": CustomInstall},
)
