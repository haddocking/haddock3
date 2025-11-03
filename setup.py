#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""Setup dot py."""
import os
import platform
import subprocess
import sys
import tarfile
from urllib.error import HTTPError, URLError
import time
import urllib.request
from pathlib import Path
from typing import cast

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from setuptools.command.develop import develop
from setuptools.command.install import install

CNS_BINARIES = {
    "x86_64-linux": "https://surfdrive.surf.nl/public.php/dav/files/BWa5OimzbNliTi6",
    "x86_64-darwin": "https://surfdrive.surf.nl/public.php/dav/files/3Fzzte0Zx0L8GTY",
    "arm64-darwin": "https://surfdrive.surf.nl/public.php/dav/files/bYB3xPWf7iwo07X",
    "aarch64-linux": "https://surfdrive.surf.nl/public.php/dav/files/3rHpxcufHGrntHn",
}

HADDOCK_RESTRAINTS_VERSION = "0.10.0"
BASE_GITHUB_URL = f"https://github.com/haddocking/haddock-restraints/releases/download/v{HADDOCK_RESTRAINTS_VERSION}"
HADDOCK_RESTRAINTS_BINARIES = {
    "x86_64-linux": f"{BASE_GITHUB_URL}/haddock-restraints-v{HADDOCK_RESTRAINTS_VERSION}-x86_64-unknown-linux-musl.tar.gz",
    "x86_64-darwin": f"{BASE_GITHUB_URL}/haddock-restraints-v{HADDOCK_RESTRAINTS_VERSION}-x86_64-apple-darwin.tar.gz",
    "arm64-darwin": f"{BASE_GITHUB_URL}/haddock-restraints-v{HADDOCK_RESTRAINTS_VERSION}-aarch64-apple-darwin.tar.gz",
    "aarch64-linux": f"{BASE_GITHUB_URL}/haddock-restraints-v{HADDOCK_RESTRAINTS_VERSION}-aarch64-unknown-linux-musl.tar.gz",
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
    """Custom build handles the C/C++ dependencies and downloading of binaries"""

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
        if hasattr(self, "build_lib") and self.build_lib and not self.inplace:
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

            if arch != "x86_64-linux" and filename.name == "cns":
                # Force the download of the linux binary, this is needed for GRID executions
                download_path = Path(target_bin_dir, "cns_linux")
                status, msg = self.download_file(
                    binary_dict["x86_64-linux"], download_path
                )
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
    def download_file(url, dest, max_retries=3, delay=2) -> tuple[bool, str]:
        """Download a file from a URL"""
        for attempt in range(max_retries):
            try:
                urllib.request.urlretrieve(url, dest)
                return True, "Download successful"
            except (URLError, HTTPError) as e:  # pylint: disable=broad-except
                if attempt < max_retries - 1:
                    time.sleep(delay)
                    continue
                return (
                    False,
                    f"Download of {dest} failed after {max_retries} attempts: {e} URL: {url}",
                )
            except Exception as e:  # pylint: disable=broad-except
                return False, f"Download of {dest} failed: {e} URL: {url}"
        return (
            False,
            f"Download of {dest} failed after {max_retries} attempts URL: {url}",
        )

    @staticmethod
    def get_arch():
        """Helper function to figure out the architetchure"""
        system = platform.system().lower()
        machine = platform.machine().lower()

        return f"{machine}-{system}"


class DevelopCommand(develop):
    """Custom develop command for editable installations"""

    # NOTE: This is what is performed when using `pip install -e`

    def run(self):
        # First run the standard develop command
        develop.run(self)

        # Then run our custom build steps
        build_cmd = cast(CustomBuild, self.get_finalized_command("build_ext"))
        # specify we are building in place
        build_cmd.inplace = True

        # Build the executables
        print("Building executables for editable installation...")
        build_cmd.build_executable(
            name="contact_fcc",
            cmd=["g++", "-O2", "-o", "contact_fcc", "src/haddock/deps/contact_fcc.cpp"],
        )
        build_cmd.build_executable(
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

        # Download binaries to source directory
        print("Downloading binaries for editable installation...")
        build_cmd.download_binaries()


class PostInstallCommand(install):
    """Post-installation for installation mode."""

    # NOTE: This is only executed for the installation, for editable see `DevelopCommand`

    def run(self):
        install.run(self)
        self.check_cns_binary()

    def check_cns_binary(self):
        """Execute the CNS binary to ensure it works"""
        cns_path = Path(f"{self.install_lib}/haddock/bin/cns")
        proc = subprocess.run(
            [cns_path],
            input="stop\n",
            text=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        if proc.returncode != 0:
            self.cns_warning()
        if proc.returncode == 0:
            self.cns_valid()

    @staticmethod
    def cns_warning():
        """Warn the user that CNS could not be executed"""
        try:
            with open("/dev/tty", "w") as tty:
                tty.write("\n" + "=" * 79 + "\n")
                tty.write("=" * 79 + "\n")
                tty.write("=" * 79 + "\n")
                tty.write(
                    "\n⚠️ WARNING: The pre-compiled CNS binary could not be executed ⚠️\n\n"
                )
                tty.write("This may be due to missing dependencies on your system.\n")
                tty.write(
                    "Please refer to the installation instructions for troubleshooting steps:\n"
                )
                tty.write(
                    "`https://github.com/haddocking/haddock3/blob/main/docs/CNS.md`\n\n"
                )
                tty.write("=" * 79 + "\n")
                tty.write("=" * 79 + "\n")
                tty.write("=" * 79 + "\n")
                tty.flush()
        except Exception as _:
            # Fallback for systems without /dev/tty
            sys.stderr.write(
                "\n⚠️ WARNING: The pre-compiled CNS binary could not be executed ⚠️\n\n"
            )

    @staticmethod
    def cns_valid():
        """Write a message to the user that CNS works"""
        try:
            with open("/dev/tty", "w") as tty:
                tty.write("\n" + "=" * 79 + "\n")
                tty.write("CNS execution passed ✅ \n")
                tty.write("=" * 79 + "\n")
        except Exception as _:
            # Fallback for systems without /dev/tty
            sys.stdout.write("CNS execution passed ✅ \n")


setup(
    cmdclass={
        "build_ext": CustomBuild,
        "install": PostInstallCommand,
        "develop": DevelopCommand,
    },
    ext_modules=cpp_extensions,
)
