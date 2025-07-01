from importlib.metadata import version, PackageNotFoundError

# Whole-package versioning
try:
    __version__ = version('autoDCR')
except PackageNotFoundError:
    __version__ = "unknown"
