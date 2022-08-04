"""HADDOCK library custom errors."""


class HaddockError(Exception):
    """Error in HADDOCK3."""

    pass


class ConfigurationError(HaddockError):
    """Error in the configuration file."""

    pass


class ModuleError(HaddockError):
    """Error in a HADDOCK3 module."""

    pass


class StepError(HaddockError):
    """Error in a HADDOCK3 workflow step."""

    pass


class JobRunningError(HaddockError):
    """General job running error."""

    pass


class CNSRunningError(HaddockError):
    """CNS run error."""

    pass


class HaddockModuleError(HaddockError):
    """General error in a HADDOCK3 module."""

    pass


class SetupError(HaddockError):
    """Set up error."""

    pass


class HaddockTermination(HaddockError):
    """Terminates HADDOCK."""

    pass
