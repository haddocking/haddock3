#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""Setup dot py."""
import os
import platform
import subprocess
import sys
import urllib.request
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


CNS_BINARIES = {
    "x86_64-linux": "https://surfdrive.surf.nl/files/index.php/s/BWa5OimzbNliTi6/download",
    "x86_64-darwin": "https://surfdrive.surf.nl/files/index.php/s/3Fzzte0Zx0L8GTY/download",
    "arm64-darwin": "https://surfdrive.surf.nl/files/index.php/s/bYB3xPWf7iwo07X/download",
    "aarch64-linux": "https://surfdrive.surf.nl/files/index.php/s/3rHpxcufHGrntHn/download",
}

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
    """Custom build handles the C/C++ dependencies"""

    def run(self):
        """Run the custom build"""
        print("Building HADDOCK3 C/C++ binary dependencies...")
        self.build_executable(
            name="contact_fcc",
            cmd=["g++", "-O2", "-o", "contact_fcc", "src/haddock/deps/contact_fcc.cpp"],
        )
        self.build_executable(
            name="fast-rmsdmatrix",
            cmd=[
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
        print("Downloading the CNS binary...")
        self.download_cns()

        # Run the standard build
        build_ext.run(self)

    def build_executable(self, name, cmd):
        """Helper function to execute the build command"""
        try:
            subprocess.check_call(cmd)
            # Ensure the source bin directory exists
            src_bin_dir = Path("src", "haddock", "bin")
            src_bin_dir.mkdir(exist_ok=True, parents=True)

            # Move the built executable to the source bin directory
            src_bin_exec = Path(src_bin_dir, name)
            if src_bin_exec.exists():
                src_bin_exec.unlink()

            self.move_file(name, src_bin_exec)

            # If build_lib exists, also copy to there
            if hasattr(self, "build_lib"):
                build_bin_dir = Path(self.build_lib, "haddock", "bin")
                build_bin_dir.mkdir(exist_ok=True, parents=True)
                self.copy_file(Path(src_bin_dir, name), Path(build_bin_dir, name))

            print(f"Successfully built and moved {name}")
        except subprocess.CalledProcessError as e:
            print(f"Error building {name}: {e}")
            raise

    def download_cns(self):
        """Helper function to download the CNS binary"""

        arch = self.get_arch()

        if arch not in CNS_BINARIES:
            print(f"Unknown architecture: {arch}")
            sys.exit(1)

        cns_binary_url = CNS_BINARIES[arch]

        src_bin_dir = Path("src", "haddock", "bin")
        src_bin_dir.mkdir(exist_ok=True, parents=True)

        cns_exec = Path(src_bin_dir, "cns")
        status, msg = self.download_file(cns_binary_url, cns_exec)
        if not status:
            print(msg)
            sys.exit(1)

        # Make it executable
        os.chmod(cns_exec, 0o755)

        # If build_lib exists, also copy to there
        if hasattr(self, "build_lib"):
            install_bin_dir = Path(self.build_lib, "haddock", "bin")
            install_bin_dir.mkdir(exist_ok=True, parents=True)
            self.copy_file(cns_exec, Path(install_bin_dir, "cns"))

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


setup(
    cmdclass={
        "build_ext": CustomBuild,
    },
    ext_modules=cpp_extensions,
)
