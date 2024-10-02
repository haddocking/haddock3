"""Set of functionalities to run prodigy."""

from abc import ABC, abstractmethod

import numpy as np

from haddock import log
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
        except ModuleNotFoundError:
            raise ModuleError(
                "Issue detected with the installation of Prodigy. "
                "Consider installing it with: "
                "pip install prodigy-prot prodigy-lig"
                )

    @staticmethod
    def run_helps() -> None:
        """Run prodigy CLI with help."""
        import prodigy_prot.predict_IC  # noqa : F401
        #import prodigy_lig.prodigy_lig  # noqa : F401
        return None


class ProdigyWorker(ABC):
    """Superclass of Prodigy Worker."""

    def __init__(self, model: FilePath, params: ParamDict) -> None:
        # Use by both prodigy -prot and -lig
        self.model = model
        self.topKd = params["to_pkd"]
        self.temperature = params["temperature"]
        self.dist_cutoff = self.set_distance_cutoff(params["distance_cutoff"])
        self.chains = params["chains"]
        # Only used by prodigy-prot
        self.acc_cutoff = params["accessibility_cutoff"]
        # Only used by prodigy-lig
        self.lig_resname = params["ligand_resname"]
        self.lig_chain = params["ligand_chain"]
        # Output values
        self.score: Optional[float] = None
        self.error: Optional[Exception] = None

    def _run(self) -> None:
        """Evaluate complex and compute score."""
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
        Kd = np.exp(-deltaG_joules / (R_GAZ_CONST * kelvins))
        # Convert to pKd
        _pKd = np.log10(Kd)
        # Reduce number of decimals
        pKd = round(_pKd, 2)
        return -pKd

    @abstractmethod
    def evaluate_complex(self) -> float:
        """Logic to evaluate a complex using prodigy."""
        pass

    @staticmethod
    @abstractmethod
    def set_distance_cutoff(_dist_cutoff: Optional[float]) -> float:
        """Logic to set the distance cutoff."""
        pass


class ProdigyProtein(ProdigyWorker):
    """Class managing the computation of protein-protein with prodigy."""

    def __init__(self, model: FilePath, params: ParamDict) -> None:
        """Instanciate the class with superclass."""
        super().__init__(model, params)

    def evaluate_complex(self) -> float:
        """Evaluate a complex with prodigy-prot.

        Returns
        -------
        deltaG : float
            The computed DeltaG of the input complex.
        """
        from prodigy_prot.predict_IC import parse_structure, Prodigy
        structure, _n_chains, _n_res = parse_structure(self.model)
        prodigy = Prodigy(structure, self.chains, self.temperature)
        prodigy.predict(
            distance_cutoff=self.dist_cutoff,
            acc_threshold=self.acc_cutoff,
            )
        deltaG = prodigy.ba_val
        return deltaG

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


class ProdigyLigand(ProdigyWorker):
    """Class managing the computation of protein-ligand with prodigy."""

    def __init__(self, model: FilePath, params: ParamDict) -> None:
        """Instanciate the class with superclass."""
        super().__init__(model, params)

    def evaluate_complex(self) -> float:
        """Evaluate a complex with prodigy-lig.

        Returns
        -------
        deltaG : float
            The computed DeltaG of the input complex.
        """
        from prodigy_lig.prodigy_lig import (
            basename,
            extract_electrostatics,
            FastMMCIFParser,
            PDBParser,
            ProdigyLig,
            splitext,
            )
        fname, s_ext = splitext(basename(self.model))
        if s_ext in (".pdb", ".ent", ):
            parser = PDBParser(QUIET=1)
        elif s_ext == ".cif":
            parser = FastMMCIFParser(QUIET=1)

        with open(self.model) as in_file:
            # try to set electrostatics from input file if not provided by user
            electrostatics = False
            try:
                electrostatics = extract_electrostatics(in_file)
            except Exception:
                pass
            prodigy_lig = ProdigyLig(
                parser.get_structure(fname, in_file),
                chains=[
                    ":".join(self.chains),
                    f"{self.lig_chain}:{self.lig_resname}",
                    ],
                electrostatics=electrostatics,
                cutoff=self.dist_cutoff,
                )
        prodigy_lig.predict()
        deltaG = prodigy_lig.dg_elec if prodigy_lig.dg_elec else prodigy_lig.dg
        return deltaG

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


class ModelScore:
    """Simple class for holding score for a model."""

    def __init__(self, model_index: int) -> None:
        """Initiate models scores."""
        self.index = model_index
        self.score = None
        self.error = None


class AnyProdigyJob:
    """Managing the computation of prodigy scores within haddock3."""

    def __init__(
            self,
            model: PDBFile,
            params: ParamDict,
            index: int = 1,
            ) -> None:
        """Initiate a worker."""
        worker = self.get_worker(params["scoring_mode"])
        self.worker = worker(model.rel_path, params)
        self.score = ModelScore(index)

    def _run(self) -> ModelScore:
        """Run the worker and retrieve output values."""
        try:
            self.worker.run()
        except Exception as e:
            log.error(e)
        else:
            self.score.score = self.worker.score
            self.score.error = self.worker.error
        return self.score

    def run(self) -> ModelScore:
        """Execute the _run method."""
        return self._run()

    @staticmethod
    def get_worker(mode: str) -> Union[ProdigyProtein, ProdigyLigand]:
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
            return ProdigyProtein
        else:
            return ProdigyLigand
