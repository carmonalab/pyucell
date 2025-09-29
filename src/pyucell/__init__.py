from importlib.metadata import version

from .knn import smooth_knn_scores
from .ranks import get_rankings
from .scoring import compute_ucell_scores

__all__ = ["get_rankings", "compute_ucell_scores", "smooth_knn_scores"]
__version__ = version("pyucell")
