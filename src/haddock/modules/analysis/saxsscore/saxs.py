from pathlib import Path
import subprocess
import math
from os import linesep

from haddock.libs.libontology import PDBFile


def run_crysol(
    atsas_path: Path, input_f: Path, saxs_data: Path, lm: float, ns: float
) -> None:
    """Use CRYSOL to calculate fit of a HADDOCK model to SAXS data."""
    
    crysol_exec = f"{atsas_path}/bin/crysol"
    cmd = f"{crysol_exec} {input_f} {saxs_data} -lm {lm} -ns {ns}"
    subprocess.call(
        cmd,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
        env={"ATSAS": atsas_path},
    )

    # TODO: Remove unecessary files?
    # suffixes_to_be_removed = [".log", ".alm", ".int", ".abs"]
    # for suffix in suffixes_to_be_removed:
    #     Path(input_f).with_suffix(suffix).unlink()


def read_chi2(fit_file: Path) -> float:
    """Get Chi^2 value from CRYSOL .fit file."""

    with open(fit_file, "r") as ff:
        fit_header = ff.readline()
        chi2_str = fit_header.partition("Chi^2:")[2].partition(r"\s")[0]
        try:
            chi2 = float(chi2_str)
        except ValueError:
            chi2 = float("nan")

    return chi2


def calculate_haddocksaxs_score(
    score: float, chi2: float, w_haddock: float, w_saxs: float) -> float:
    """Calculate HADDOCKsaxs based on previous HADDOCK score and Chi^2."""

    chi = math.sqrt(chi2)
    return (score * w_haddock) + (chi * w_saxs)


def output_saxsscore(output_models: list[PDBFile], output_fname: Path) -> None:
    """Write Chi^2 values and HADDOCKsaxs scores to output file."""

    text_generator = (
        f"{pdb.file_name}\t{pdb.ori_name}\t{pdb.md5}\t{pdb.chi2}\t{pdb.score}"
        for pdb in output_models
    )

    with open(output_fname, "w") as fh:
        fh.write(
            "\t".join(("structure", "original_name", "md5", "Chi^2", "score"))
            + linesep
        )
        fh.write(linesep.join(text_generator))


def test_saxsmodule_run():
    pass
