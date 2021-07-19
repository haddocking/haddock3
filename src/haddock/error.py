"""HADDOCK library custom errors"""


class HaddockError(Exception):
    pass


class ConfigurationError(HaddockError):
    pass


class RecipeError(HaddockError):
    pass


class CNSRunningError(HaddockError):
    pass


class HaddockModuleError(HaddockError):
    pass
