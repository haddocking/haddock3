import os
import platform
import subprocess
import sys
import urllib.request
from pathlib import Path

CNS_BINARIES = {
    "x86_64-linux": "https://surfdrive.surf.nl/files/index.php/s/BWa5OimzbNliTi6/download",
    "x86_64-darwin": "https://surfdrive.surf.nl/files/index.php/s/3Fzzte0Zx0L8GTY/download",
    "arm64-darwin": "https://surfdrive.surf.nl/files/index.php/s/bYB3xPWf7iwo07X/download",
    "aarch64-linux": "https://surfdrive.surf.nl/files/index.php/s/3rHpxcufHGrntHn/download",
}


def build_executable(name, cmd):
    try:
        subprocess.check_call(cmd)
        src_bin_dir = Path("src", "haddock", "bin")
        src_bin_dir.mkdir(exist_ok=True, parents=True)
        src_bin_exec = Path(src_bin_dir, name)
        if src_bin_exec.exists():
            src_bin_exec.unlink()
        os.rename(name, src_bin_exec)
        print(f"Successfully built and moved {name}")
    except subprocess.CalledProcessError as e:
        print(f"Error building {name}: {e}")
        raise


def download_cns():
    arch = get_arch()
    if arch not in CNS_BINARIES:
        print(f"Unknown architecture: {arch}")
        sys.exit(1)
    cns_binary_url = CNS_BINARIES[arch]
    src_bin_dir = Path("src", "haddock", "bin")
    src_bin_dir.mkdir(exist_ok=True, parents=True)
    cns_exec = Path(src_bin_dir, "cns")
    status, msg = download_file(cns_binary_url, cns_exec)
    if not status:
        print(msg)
        sys.exit(1)
    os.chmod(cns_exec, 0o755)


def download_file(url, dest):
    try:
        urllib.request.urlretrieve(url, dest)
        return True, "Download successful"
    except Exception as e:
        return False, f"Download failed: {e}"


def get_arch():
    system = platform.system().lower()
    machine = platform.machine().lower()
    return f"{machine}-{system}"


def build():
    print("Building HADDOCK3 C/C++ binary dependencies...")
    build_executable(
        name="contact_fcc",
        cmd=["g++", "-O2", "-o", "contact_fcc", "src/haddock/deps/contact_fcc.cpp"],
    )
    build_executable(
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
    download_cns()


if __name__ == "__main__":
    build()
