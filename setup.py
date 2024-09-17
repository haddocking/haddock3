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

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install


cpp_extensions = [
    Extension(
        "haddock.bin.contact_fcc",
        sources=["src/haddock/deps/contact_fcc.cpp"],
        extra_compile_args=["-O2"],
    ),
    Extension(
        "haddock.bin.fast_rmsdmatrix",
        sources=["src/haddock/deps/fast-rmsdmatrix.c"],
        extra_compile_args=["-Wall", "-O3", "-march=native", "-std=c99"],
        extra_link_args=["-lm"],
    ),
]


class CustomBuild(build_ext):
    def run(self):
        print("Building HADDOCK3 C/C++ binary dependencies...")
        self.build_executable(
            "contact_fcc",
            ["g++", "-O2", "-o", "contact_fcc", "src/haddock/deps/contact_fcc.cpp"],
        )
        self.build_executable(
            "fast-rmsdmatrix",
            [
                "gcc",
                "-Wall",
                "-O3",
                "-march=native",
                "-std=c99",
                "-o",
                "fast-rmsdmatrix",
                "src/haddock/deps/fast-rmsdmatrix.c",
                "-lm",
            ],
        )

        # Run the standard build_ext
        build_ext.run(self)

    def build_executable(self, name, cmd):
        try:
            subprocess.check_call(cmd)
            bin_dir = os.path.join("src", "haddock", "bin")
            os.makedirs(bin_dir, exist_ok=True)
            shutil.move(name, os.path.join(bin_dir, name))
            print(f"Successfully built and moved {name}")
        except subprocess.CalledProcessError as e:
            print(f"Error building {name}: {e}")
            raise


CNS_BINARIES = {
    "x86_64-linux": "https://surfdrive.surf.nl/files/index.php/s/BWa5OimzbNliTi6/download",
    "x86_64-darwin": "https://surfdrive.surf.nl/files/index.php/s/3Fzzte0Zx0L8GTY/download",
    "arm64-darwin": "https://surfdrive.surf.nl/files/index.php/s/bYB3xPWf7iwo07X/download",
    "aarch64-linux": "https://surfdrive.surf.nl/files/index.php/s/3rHpxcufHGrntHn/download",
}


class CustomInstall(install):
    """Custom class to handle the download of the CNS binary"""

    def run(self):
        """Run the installation"""

        # Run the standard installation
        install.run(self)

        # Get the installation directory
        if self.install_lib is None:
            print("Something went wrong during installation.")
            sys.exit(1)

        # Set where the cns binary needs to be
        bin_dir = Path(self.install_lib, "haddock", "bin")
        bin_dir = Path("src", "haddock", "bin")

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
    package_data={"haddock": ["bin/*"]},
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
        ]
    },
    cmdclass={"build_ext": CustomBuild, "install": CustomInstall},
    ext_modules=cpp_extensions,
)
