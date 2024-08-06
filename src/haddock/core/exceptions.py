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


class DependencyError(ModuleError):
    """Error throw when required dependency not satisfied."""

    def __init__(
            self,
            msg: str = "",
            module: str = "",
            dependency: str = "",
            ):
        self.message = msg
        self.module = module
        self.dependency = dependency

    def __str__(self) -> str:
        additions: str = ""
        if self.module:
            additions += f"Module `{self.module}` -> "
        if self.dependency:
            additions += f"Required dependency `{self.dependency}`"
        return f"{self.message} {additions}"
