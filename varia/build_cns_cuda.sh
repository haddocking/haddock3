#!/usr/bin/env bash
# ==============================================================================
# build_cns_cuda.sh: Automated build script for cns_solve_CUDA in HADDOCK3
# ==============================================================================
# This script compiles the CUDA-accelerated CNS binary (cns_solve_CUDA)
# by detecting target GPU architectures, configuring nvcc and gfortran flags,
# and linking the non-bonded energy kernels.
#
# Reference: https://github.com/chitosan/cns_solve_CUDA
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HADDOCK3_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BIN_DIR="${HADDOCK3_ROOT}/src/haddock/bin"
TARGET_BIN="${BIN_DIR}/cns_solve_CUDA"
BUILD_DIR="${HADDOCK3_ROOT}/build/cns_cuda"

echo "================================================================="
echo " HADDOCK3 cns_solve_CUDA Build Tool"
echo "================================================================="

# 1. Check prerequisite compilers
if ! command -v nvcc &> /dev/null; then
    echo "[-] Error: nvcc (NVIDIA CUDA Compiler) not found in PATH."
    echo "    Please load the CUDA module or install the CUDA Toolkit."
    exit 1
fi

if ! command -v gfortran &> /dev/null; then
    echo "[-] Error: gfortran (GNU Fortran Compiler) not found in PATH."
    echo "    Please install gfortran (e.g., sudo apt install gfortran)."
    exit 1
fi

NVCC_VERSION=$(nvcc --version | grep "release" | awk '{print $5}' | tr -d ',')
GFORTRAN_VERSION=$(gfortran -dumpversion)
echo "[+] Detected nvcc version: ${NVCC_VERSION}"
echo "[+] Detected gfortran version: ${GFORTRAN_VERSION}"

# 2. Detect target GPU architecture
CUDA_ARCH_FLAGS=""
if [[ -n "${CUDA_ARCH:-}" ]]; then
    echo "[+] Using user-specified CUDA_ARCH: ${CUDA_ARCH}"
    CUDA_ARCH_FLAGS="-gencode arch=compute_${CUDA_ARCH},code=sm_${CUDA_ARCH}"
elif command -v nvidia-smi &> /dev/null; then
    COMPUTE_CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n1 | tr -d '.')
    if [[ -n "${COMPUTE_CAP}" ]]; then
        echo "[+] Detected active GPU compute capability: sm_${COMPUTE_CAP}"
        CUDA_ARCH_FLAGS="-gencode arch=compute_${COMPUTE_CAP},code=sm_${COMPUTE_CAP}"
    fi
fi

if [[ -z "${CUDA_ARCH_FLAGS}" ]]; then
    echo "[+] Defaulting to multi-architecture compilation (T4/A10G/A100/H100):"
    echo "    sm_75 (T4), sm_80 (A100), sm_86 (A10G/RTX3090), sm_89 (RTX4090/L4), sm_90 (H100)"
    CUDA_ARCH_FLAGS="\
-gencode arch=compute_75,code=sm_75 \
-gencode arch=compute_80,code=sm_80 \
-gencode arch=compute_86,code=sm_86 \
-gencode arch=compute_89,code=sm_89 \
-gencode arch=compute_90,code=sm_90"
fi

# 3. Setup build directory
mkdir -p "${BUILD_DIR}"
mkdir -p "${BIN_DIR}"
cd "${BUILD_DIR}"

# 4. Clone or update cns_solve_CUDA repository if not provided locally
if [[ ! -d "cns_solve_CUDA" ]]; then
    echo "[+] Cloning cns_solve_CUDA source repository..."
    git clone --depth 1 https://github.com/chitosan/cns_solve_CUDA.git
fi

cd cns_solve_CUDA

# 5. Compile CUDA acceleration source
echo "[+] Compiling CUDA nonbonded kernels with nvcc..."
NVCC_FLAGS="-O3 --use_fast_math -Xcompiler -fPIC ${CUDA_ARCH_FLAGS}"

# Check for CUDA source file in repository
if [[ -f "cns_cuda_nb.cu" ]]; then
    nvcc ${NVCC_FLAGS} -c cns_cuda_nb.cu -o cns_cuda_nb.o
elif [[ -f "cuda_nb.cu" ]]; then
    nvcc ${NVCC_FLAGS} -c cuda_nb.cu -o cns_cuda_nb.o
else
    # Find any .cu files
    CU_FILES=$(find . -maxdepth 2 -name "*.cu")
    if [[ -n "${CU_FILES}" ]]; then
        for cu_f in ${CU_FILES}; do
            obj_name="$(basename "${cu_f}" .cu).o"
            echo "    Compiling ${cu_f} -> ${obj_name}"
            nvcc ${NVCC_FLAGS} -c "${cu_f}" -o "${obj_name}"
        done
    else
        echo "[-] Note: No .cu files found in cns_solve_CUDA clone root."
    fi
fi

# 6. Build CNS binary if Makefile exists
if [[ -f "Makefile" ]]; then
    echo "[+] Running Makefile to link cns_solve_CUDA executable..."
    make -j"$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)"
fi

# 7. Locate generated binary and link to haddock bin
COMPILED_BIN=$(find . -maxdepth 3 -type f -name "cns_solve_CUDA" -o -name "cns_solve_cuda.bin" | head -n1)
if [[ -n "${COMPILED_BIN}" && -f "${COMPILED_BIN}" ]]; then
    cp "${COMPILED_BIN}" "${TARGET_BIN}"
    chmod +x "${TARGET_BIN}"
    echo "[+] Successfully installed cns_solve_CUDA to: ${TARGET_BIN}"
    echo "    Exporting environment: export CNS_CUDA_EXEC=${TARGET_BIN}"
else
    echo "[!] Warning: Built components ready, but primary binary not assembled automatically."
    echo "    Please verify Makefile configuration in ${BUILD_DIR}/cns_solve_CUDA."
fi

echo "================================================================="
echo " Build process finished."
echo "================================================================="
