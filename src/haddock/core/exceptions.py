"""HADDOCK library custom errors"""


class HaddockError(Exception):
    pass


class ConfigurationError(HaddockError):
    pass


class ModuleError(HaddockError):
    pass


class StepError(HaddockError):
    pass


class JobRunningError(HaddockError):
    pass


class CNSRunningError(HaddockError):
    pass


class HaddockModuleError(HaddockError):
    pass


class SetupError(HaddockError):
    pass
