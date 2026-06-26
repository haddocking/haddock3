"""Set of functionalities to run prodigy."""

import logging
from abc import ABC, abstractmethod

import numpy as np

from haddock import log
from haddock.core.exceptions import ModuleError
from haddock.core.typing import FilePath, Optional, ParamDict
from haddock.libs.libontology import PDBFile


# Define conversion constants
R_GAS_CONST = 8.31446261815324 	# (J x mol−1 x K−1)
CALS_TO_JOULES = 4.184
KELVIN_TO_CELCIUS = 273.15


class CheckInstall:
    """Verify that installation of prodigy is ok."""

    def __init__(self):
        """Verify that installation of prodigy is ok.

        Raises
        ------
        ModuleError
            Raised when prodigy libraries could not be loaded.
        """
        try:
            self.import_prodigy_libs()
        except ModuleNotFoundError:
            raise ModuleError(
                "Issue detected with the installation of Prodigy. "
                "Consider installing it with: "
                "pip install prodigy-prot prodigy-lig"
                )

    @staticmethod
    def import_prodigy_libs() -> None:
        """Import prodigy prot and lig libraries."""
        import prodigy_prot  # noqa : F401
        import prodigy_lig  # noqa : F401
        return None


class ProdigyWorker(ABC):
    """Prodigy Worker class."""

    def __init__(self, model: FilePath, params: ParamDict) -> None:
        # Use by both prodigy -prot and -lig
        self.model = model
        self.topKd = params["to_pkd"]
        self.temperature = params["temperature"]
        self.dist_cutoff = params["distance_cutoff"]
        # Output values
        self.score: Optional[float] = None
        self.error: Optional[Exception] = None

    def _run(self) -> None:
        """Evaluate complex and compute score."""
        # Patching the logging to not print everything on screen
        logging.getLogger("Prodigy").setLevel(logging.CRITICAL)
        # Running the predictions
        self.score = self.pkd_converter(self.evaluate_complex())

    def run(self) -> None:
        """Execute the _run method."""
        try:
            self._run()
        except ModuleError as error:
            self.error = error

    def pkd_converter(self, deltaG: float) -> float:
        """Decide if deltaG must be converted to pKd.

        Parameters
        ----------
        deltaG : float
            Input DeltaG value

        Returns
        -------
        score : float
            The converted DeltaG to pKd if self.topKd is true,
            else the input DeltaG.
        """
        if self.topKd:
            score = self.deltaG_to_pKd(
                deltaG,
                kelvins=self.temperature + KELVIN_TO_CELCIUS,
                )
        else:
            score = deltaG
        return round(score, 3)
    
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
        Kd = np.exp(-deltaG_joules / (R_GAS_CONST * kelvins))
        # Convert to pKd
        _pKd = np.log10(Kd)
        # Reduce number of decimals
        pKd = round(_pKd, 2)
        return -pKd

    @abstractmethod
    def evaluate_complex(self) -> float:
        """Logic to evaluate a complex using prodigy."""
        pass


class ModelScore:
    """Simple class for holding score for a model."""

    def __init__(self, model_index: int) -> None:
        """Initiate models scores."""
        self.index = model_index
        self.score = None
        self.error = None


class ProdigyBaseJob(ABC):
    """Managing the computation of prodigy scores within haddock3."""

    def __init__(
            self,
            model: PDBFile,
            params: ParamDict,
            index: int = 1,
            ) -> None:
        """Initiate a worker."""
        self.score = ModelScore(index)
        worker = self.get_worker()
        self.worker = worker(model.rel_path, params)

    def _run(self) -> ModelScore:
        """Run the worker and retrieve output values."""
        try:
            self.worker.run()
        except Exception as e:
            log.error(e)
            self.worker.error = str(e)
            self.score.error = str(e)
        else:
            self.score.score = self.worker.score
            self.score.error = self.worker.error
        return self.score

    def run(self) -> ModelScore:
        """Execute the _run method."""
        return self._run()

    @staticmethod
    @abstractmethod
    def get_worker() -> object:
        """Return the appropriate worker."""
        pass
