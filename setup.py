#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""Setup dot py."""
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.request
import warnings
from os.path import dirname, join
from pathlib import Path

from setuptools import SetuptoolsDeprecationWarning, find_packages, setup
from setuptools.command.build_ext import build_ext
from setuptools.command.easy_install import EasyInstallDeprecationWarning
from setuptools.command.install import install


# Add warnings filtering to the Setup Deprecation Warnings
warnings.filterwarnings("ignore", category=SetuptoolsDeprecationWarning)
warnings.filterwarnings("ignore", category=EasyInstallDeprecationWarning)

with open("requirements.txt") as f:
    requirements = f.read().splitlines()


def read(*names, **kwargs) -> str:
    """Read description files."""
    path = join(dirname(__file__), *names)
    with open(path, encoding=kwargs.get("encoding", "utf8")) as fh:
        return fh.read()


# activate once added, do not remove
long_description = "{}\n{}".format(
    read("README.md"),
    read("CHANGELOG.md"),
)


CNS_BINARY_URL = ""


class CustomBuildInstall(install):
    def run(self):
        self.run_command("build_ext")
        install.run(self)


class CustomBuild(build_ext):
    def run(self):
        # Ensure bin directory exists
        self.bin_dir = os.path.join(self.get_install_dir(), "haddock", "bin")
        os.makedirs(self.bin_dir, exist_ok=True)

        # Build FCC
        self.build_submodule(
            name="FCC",
            repo_url="https://github.com/haddocking/fcc.git",
            build_cmd="make",
            src_dir="src",
            binary_name="contact_fcc",
        )

        # Build fast-rmsdmatrix
        self.build_submodule(
            name="fast-rmsdmatrix",
            repo_url="https://github.com/mgiulini/fast-rmsdmatrix.git",
            build_cmd="make fast-rmsdmatrix",
            src_dir="src",
            binary_name="fast-rmsdmatrix",
        )

        # Run original build_ext
        build_ext.run(self)

    def get_install_dir(self):
        install_cmd = self.get_finalized_command("install")
        return install_cmd.install_lib

    def build_submodule(self, name, repo_url, build_cmd, src_dir, binary_name):
        print(f"Building {name}...")
        with tempfile.TemporaryDirectory() as temp_dir:

            # clone the repository
            subprocess.run(["git", "clone", repo_url, temp_dir], check=True)

            # Build
            build_dir = Path(temp_dir, src_dir)
            subprocess.run(
                build_cmd,
                shell=True,
                cwd=build_dir,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            print("pass!")

            # Move the binary
            src_bin = Path(build_dir, binary_name)
            dst_bin = Path(self.bin_dir, binary_name)

            shutil.copy2(src_bin, dst_bin)

            print(dst_bin)


class CNSInstall(CustomBuildInstall):
    """Custom class to handle the download of the CNS binary"""

    def run(self):

        CustomBuildInstall.run(self)

        # Run the standard installation
        # install.run(self)

        # Get the installation directory
        if self.install_lib is None:
            # Something went wrong with the installation
            sys.exit(1)

        bin_dir = Path(self.install_lib, "haddock", "bin")
        cns_exec = Path(bin_dir, "cns")

        # Create the `bin/` directory
        bin_dir.mkdir(exist_ok=True)

        # Download the binary
        if cns_exec.exists():
            cns_exec.unlink()

        urllib.request.urlretrieve(CNS_BINARY_URL, cns_exec)

        os.chmod(cns_exec, 0o755)


setup(
    name="haddock3",
    version="3.0.0",
    description="Haddock 3.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="Apache License 2.0",
    author="HADDOCK",
    author_email="A.M.J.J.Bonvin@uu.nl",
    url="https://github.com/haddocking/haddock3",
    packages=find_packages("src"),
    package_dir={"": "src"},
    # py_modules=[splitext(basename(i))[0] for i in glob("src/*.py")],
    include_package_data=True,
    # package_data={"haddock3": ["fcc/src/*", "fast-rmsdmatrix/src/*"]},
    zip_safe=False,
    classifiers=[
        # complete classifier list:
        # http://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Operating System :: POSIX",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
        "Programming Language :: Python :: 3.9",
    ],
    project_urls={
        "webpage": "https://github.com/haddocking/haddock3",
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
    cmdclass={"build_ext": CustomBuild, "install": CNSInstall},
    # cmdclass={'build_ext': optional_build_ext},
    # ext_modules=[
    #    Extension(
    #        splitext(relpath(path, 'src').replace(os.sep, '.'))[0],
    #        sources=[path],
    #        include_dirs=[dirname(path)]
    #    )
    #    for root, _, _ in os.walk('src')
    #    for path in glob(join(root, '*.c'))
    # ],
)
