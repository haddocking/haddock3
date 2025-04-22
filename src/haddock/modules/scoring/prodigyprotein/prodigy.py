"""Set of functionalities to run prodigy-protein.

DevNotes:
The magic happens in haddock/libs/libprodigy.py
"""

from haddock.core.typing import FilePath, Optional, ParamDict
from haddock.libs.libontology import PDBFile
from haddock.libs.libprodigy import ProdigyBaseJob, ProdigyWorker


class ProdigyProtein(ProdigyWorker):
    """Class managing the computation of protein-protein with prodigy."""

    def __init__(self, model: FilePath, params: ParamDict) -> None:
        """Instantiate the class with superclass."""
        super().__init__(model, params)
        # Hold specific parameters
        self.acc_cutoff = params["accessibility_cutoff"]
        self.chains = self.convert_chains(params["chains"])

    @staticmethod
    def convert_chains(chains: list[str]) -> Optional[list[str]]:
        """Read and convert the `chains` parameter.

        Parameters
        ----------
        chains : Optional[list[str]]
            Content of the chains parameter as provided by the user.

        Returns
        -------
        Optional[list[str]]
            None or list of strings containing the chains.
        """
        if len(chains) == 0:
            return None
        else:
            # Removes potential spaces if user is inputing something wrong
            # e.g.: "A, B" needs to be "A,B"
            return [c.replace(" ", "") for c in chains]

    def evaluate_complex(self) -> float:
        """Evaluate a complex with prodigy-prot.

        Returns
        -------
        deltaG : float
            The computed DeltaG of the input complex.
        """
        from prodigy_prot.modules.prodigy import Prodigy
        from prodigy_prot.modules.parsers import parse_structure
        structure, _n_chains, _n_res = parse_structure(self.model)
        prodigy = Prodigy(structure, self.chains, self.temperature)
        prodigy.predict(
            distance_cutoff=self.dist_cutoff,
            acc_threshold=self.acc_cutoff,
            )
        deltaG = prodigy.ba_val
        return deltaG


class ProdigyJob(ProdigyBaseJob):
    """Managing the computation of prodigy scores within haddock3."""

    def __init__(
            self,
            model: PDBFile,
            params: ParamDict,
            index: int = 1,
            ) -> None:
        """Instantiate the class with superclass."""
        super().__init__(model, params, index)
    
    @staticmethod
    def get_worker() -> object:
        """Return the appropriate worker."""
        return ProdigyProtein
