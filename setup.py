#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""Setup dot py."""
import os
import platform
import subprocess
import sys
import tarfile
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

HADDOCK_RESTRAINTS_VERSION = "0.9.0"
BASE_GITHUB_URL = f"https://github.com/haddocking/haddock-restraints/releases/download/v{HADDOCK_RESTRAINTS_VERSION}"
HADDOCK_RESTRAINTS_BINARIES = {
    "x86_64-linux": f"{BASE_GITHUB_URL}/haddock-restraints-v{HADDOCK_RESTRAINTS_VERSION}-x86_64-unknown-linux-gnu.tar.gz",
    "x86_64-darwin": f"{BASE_GITHUB_URL}/haddock-restraints-v{HADDOCK_RESTRAINTS_VERSION}-x86_64-apple-darwin.tar.gz",
    "arm64-darwin": f"{BASE_GITHUB_URL}/haddock-restraints-v{HADDOCK_RESTRAINTS_VERSION}-aarch64-apple-darwin.tar.gz",
    "aarch64-linux": f"{BASE_GITHUB_URL}/haddock-restraints-v{HADDOCK_RESTRAINTS_VERSION}-aarch64-unknown-linux-gnu.tar.gz",
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
        # print("Building HADDOCK3 C/C++ binary dependencies...")
        # self.build_executable(
        #     name="contact_fcc",
        #     cmd=["g++", "-O2", "-o", "contact_fcc", "src/haddock/deps/contact_fcc.cpp"],
        # )
        # self.build_executable(
        #     name="fast-rmsdmatrix",
        #     cmd=[
        #         "gcc",
        #         "-Wall",
        #         "-O3",
        #         "-march=native",
        #         "-std=c99",
        #         "-o",
        #         "fast-rmsdmatrix",
        #         "src/haddock/deps/fast-rmsdmatrix.c",
        #         "-lm",
        #     ],
        # )
        print("Downloading the binaries...")
        self.download_binaries()

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

    def download_binaries(self):
        """Helper function to download binaries"""
        arch = self.get_arch()
        binary_configs = [
            (HADDOCK_RESTRAINTS_BINARIES, Path("haddock-restraints.tar.gz")),
            (CNS_BINARIES, Path("cns")),
        ]

        # Determine the target directory
        if hasattr(self, "build_lib"):
            target_bin_dir = Path(self.build_lib, "haddock", "bin")
        else:
            target_bin_dir = Path("src", "haddock", "bin")

        target_bin_dir.mkdir(exist_ok=True, parents=True)

        for binary_dict, filename in binary_configs:
            if arch not in binary_dict:
                print(f"Unknown architecture: {arch}")
                sys.exit(1)
            binary_url = binary_dict[arch]

            download_path = Path(target_bin_dir, filename.name)
            status, msg = self.download_file(binary_url, download_path)
            if not status:
                print(msg)
                sys.exit(1)

            if "".join(filename.suffixes) == ".tar.gz":
                try:
                    with tarfile.open(download_path, "r:gz") as tar:
                        tar.extractall(path=target_bin_dir)
                    # Remove the compressed file
                    download_path.unlink()
                    # The executable should be the filename without .tar.gz suffixes
                    executable = Path(
                        target_bin_dir, filename.with_suffix("").with_suffix("").name
                    )
                except tarfile.TarError as e:
                    print(f"Error extracting {download_path}: {e}")
                    sys.exit(1)
            else:
                executable = download_path

            print(f"Downloaded {filename} to {executable}")

            # Make it executable
            os.chmod(executable, 0o755)

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
