from importlib.metadata import version

from .ranks import get_rankings
from .scoring import compute_ucell_scores

__all__ = [
    "get_rankings",
    "compute_ucell_scores"
]
__version__ = version("pyucell")
