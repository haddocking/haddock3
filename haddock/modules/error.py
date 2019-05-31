"""Library errors"""


class HaddockError(Exception):
    """Main library error class"""
    pass


class CNSError(HaddockError):
    """CNS Engine errors"""
    pass
