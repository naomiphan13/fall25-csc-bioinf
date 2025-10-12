class InvalidFileError(Exception):
    """
    Indicates that the file is not suitable for the requested action,
    either because the file does not contain the required data or
    because the file is malformed.
    """

    pass
