"""Set of functionalities to run prodigy."""


import subprocess
from abc import ABC, abstractmethod

import numpy as np

from haddock.core.exceptions import ModuleError
from haddock.core.typing import FilePath, Optional, ParamDict, Union
from haddock.libs.libontology import PDBFile


# Define conversion constants
R_GAZ_CONST = 8.314
CALS_TO_JOULES = 4.184
KELVIN_TO_CELCIUS = 273.15


class CheckInstall:
    """Verify that installation of prodigy is ok."""

    def __init__(self):
        """Verify that installation of prodigy is ok.

        Raises
        ------
        ModuleError
            Raised when prodigy subprocess failed.
        """
        try:
            self.run_helps()
        except FileNotFoundError:
            raise ModuleError(
                "Issue detected with the installation of Prodigy. "
                "Consider installing it with: "
                "pip install prodigy-prot prodigy-lig"
                )

    @staticmethod
    def run_helps() -> None:
        """Run prodigy CLI with help."""
        for ext in ("",):#"-lig"):
            prodigy_help = subprocess.Popen(
                [f"prodigy{ext}", "--help"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                )
            _help_string = prodigy_help.stdout.read().decode('utf-8')
        return None


class ProdigyWorker(ABC):
    def __init__(self, model: FilePath, params: ParamDict):
        self.model = model
        self.topKd = params["to_pkd"]
        self.temperature = params["temperature"]
        self.acc_cutoff = params["accessibility_cutoff"]
        self.dist_cutoff = self.set_distance_cutoff(params["distance_cutoff"])
        self.score = None
        self.error = None

    def _run(self) -> None:
        """Main computation function."""
        prodigy_stdout = self.evaluate_complex()
        deltaG = self.extract_deltaG(prodigy_stdout)
        self.score = self.pkd_converter(deltaG)

    def run(self) -> None:
        """Wrapper of the run function."""
        try:
            self._run()
        except ModuleError as error:
            self.error = error
    
    def evaluate_complex(self) -> str:
        """Run prodigy and return stdout.

        Returns
        -------
        prodigy_stdout : str
            Standard output of prodigy (in quiet mode).
        """
        # Build command line
        cmd_ = self.gen_cmd_line()
        # Subprocess run
        prodigy_stdout = self.subprocess_run(cmd_)
        return prodigy_stdout
    
    def pkd_converter(self, deltaG: float) -> float:
        if self.topKd:
            return self.deltaG_to_pKd(
                deltaG,
                kelvins=self.temperature + KELVIN_TO_CELCIUS,
                )
        else:
            return deltaG
    
    @staticmethod
    def deltaG_to_pKd(deltaG: float, kelvins: float = 298.3) -> float:
        """Convert a deltaG (in Kcal/mol) to pKd.

        Parameters
        ----------
        deltaG : float
            DeltaG value (in Kcal/mol)
        kelvins : float, optional
            Temperature at which to perform the conversion.
            By default 298.3 (25 degrees Celcius)

        Returns
        -------
        pKd : float
            The corresponding pKd value at given temperature.
        """
        # Convert Kcals to joules
        deltaG_joules = deltaG * CALS_TO_JOULES * 1000
        # Use formula
        Kd = np.exp(-deltaG_joules / (R_GAZ_CONST * kelvins))
        # Convert to pKd
        _pKd = np.log10(Kd)
        # Reduce number of decimals
        pKd = round(_pKd, 2)
        return -pKd

    @staticmethod
    def subprocess_run(cmd_: list[Union[str, FilePath]]) -> Optional[str]:
        print(" ".join([str(v) for v in cmd_]))
        # Run it
        run = subprocess.Popen(
            cmd_,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            )
        # Check that run went smooth
        stderr = run.stderr.read().decode("utf-8")
        if stderr != "":
            raise ModuleError(stderr)
        # Gather standard output
        stdout = run.stdout.read().decode("utf-8")
        return stdout

    @abstractmethod
    def set_distance_cutoff(_dist_cutoff: Optional[float]) -> float:
        pass

    @abstractmethod
    def extract_deltaG(stdout: str) -> str:
        pass

    @abstractmethod
    def gen_cmd_line(self) -> list[Union[str, FilePath]]:
        pass


class ProdigyProt(ProdigyWorker):
    def __init__(self, model: FilePath, params: ParamDict):
        super().__init__(model, params)

    def gen_cmd_line(self) -> list[Union[str, FilePath]]:
        # Build command line
        cmd_ = [
            "prodigy",
            self.model,
            "--distance-cutoff", str(self.dist_cutoff),
            "--acc-threshold", str(self.acc_cutoff),
            "--temperature", str(self.temperature),
            "--quiet",
            ]
        return cmd_

    @staticmethod
    def set_distance_cutoff(_dist_cutoff: Optional[float]) -> float:
        """Set distance cutoff to use for analysis of contacts.

        Parameters
        ----------
        _dist_cutoff : Optional[float]
            Distance set by the user. Can be None if undefined.

        Returns
        -------
        dist_cutoff : float
            The distance cutoff to use for this run.
            Set to default software value (5.5) if undefined.
        """
        # User defined value
        if not np.isnan(_dist_cutoff):
            dist_cutoff = _dist_cutoff
        # Not defined case
        else:
            dist_cutoff = 5.5
        return dist_cutoff

    @staticmethod
    def extract_deltaG(prodigy_predictions: str) -> float:
        """Extract score from prodigy-prot

        This function is ment to be used when the flag --quiet is used.

        Parameters
        ----------
        prodigy_predictions : str
            Standard output of prodigy (prodigy-protein)

        Returns
        -------
        deltaG: float
            The predicted deltaG
        """
        deltaG_str = prodigy_predictions.strip().split()[-1]
        deltaG = float(deltaG_str)
        return round(deltaG, 3)


class ProdigyLig(ProdigyWorker):
    def __init__(self, model: FilePath, params: ParamDict):
        super().__init__(model, params)

    def gen_cmd_line(self) -> list[Union[str, FilePath]]:
        # Build command line
        cmd_ = [
            "prodigy_lig",
            "--chains", ",".join(chains), f"{ligchain}:{ligname}",
            "--distance_cutoff", self.dist_cutoff,
            "--input_file", self.model,
            ]
        return cmd_

    @staticmethod
    def set_distance_cutoff(_dist_cutoff: Optional[float]) -> float:
        """Set distance cutoff to use for analysis of contacts.

        Parameters
        ----------
        _dist_cutoff : Optional[float]
            Distance set by the user. Can be None if undefined.

        Returns
        -------
        dist_cutoff : float
            The distance cutoff to use for this run.
            Set to default software value (10.5) if undefined.
        """
        # User defined value
        if not np.isnan(_dist_cutoff):
            dist_cutoff = _dist_cutoff
        # Not defined case
        else:
            dist_cutoff = 10.5
        return dist_cutoff

    @staticmethod
    def extract_deltaG(prodigy_predictions: str) -> float:
        """Extract delatG from prodigy-lig predictions.

        Parameters
        ----------
        prodigy_lig_predictions : str
            Standard output of prodigy-lig.

        Returns
        -------
        deltaG : float
            The predicted deltaG
        """
        # Point line and value
        prediction_lines = prodigy_predictions.split('\n')
        second_line = prediction_lines[1].strip()
        predicted_DG_string = second_line.split()[-1]
        # Cast to float
        _deltaG = float(predicted_DG_string)
        # Reduce number of decimals
        deltaG = round(_deltaG, 2)
        return deltaG


class ModelScore:
    """Simple class for holding score for a model."""

    def __init__(self, model_index: int):
        self.index = model_index
        self.score = None
        self.error = None


class AnyProdigyJob:
    """Managing the computation of prodigy scores within haddock3."""

    def __init__(self, model: PDBFile, params: ParamDict, index: int = 1):
        worker = self.get_worker(params["scoring_mode"])
        self.worker = worker(model.rel_path, params)
        self.score = ModelScore(index)

    def _run(self) -> ModelScore:
        """Main computation function."""
        try:
            self.worker.run()
        except Exception as e:
            print(e)
        else:
            self.score.score = self.worker.score
            self.score.error = self.worker.error
        return self.score

    def run(self) -> ModelScore:
        """Wrapper to the _run function."""
        return self._run()

    @staticmethod
    def get_worker(mode: str) -> Union[ProdigyProt, ProdigyLig]:
        """Find appropriated class that will handle the computation.

        Parameters
        ----------
        mode : str
            Prodigy scoring mode, either protein or ligand

        Returns
        -------
        Union[ProdigyProt, ProdigyLig]
            The prodigy class that will handle the computation
        """
        if mode == "protein":
            return ProdigyProt
        else:
            return ProdigyLig
