from importlib.metadata import version

from .ranks import get_rankings

__all__ = [
    "get_rankings"
]
__version__ = version("pyucell")
