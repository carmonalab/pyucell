import numpy as np
import scanpy as sc
from scipy import sparse


def smooth_knn_scores(
    adata, obs_columns, k=10, use_rep="X_pca", decay: float = 0.1, up_only: bool = False, graph_key=None, suffix="_kNN"
):
    """
    Smooth per-cell signature scores by weight-averaging over k nearest neighbors.

    Parameters
    ----------
    adata : AnnData
        AnnData with .obs columns containing signature scores.
    obs_columns : list of str
        Names of columns in adata.obs with the scores to smooth.
    k : int
        Number of neighbors (ignored if graph_key provided).
    use_rep : str
        Representation to compute neighbors (e.g., 'X_pca', 'X').
    decay : float
        Decay parameter (0<decay<1) for weighting nearest neighbors:
        weight = (1 - decay)^i for i-th neighbor.
    up_only : bool
        If True, ensures smoothed scores are at least as high as the original.
    graph_key : str or None
        If provided, the name of a precomputed graph in adata.obsp to use
        (e.g. 'connectivities' or 'my_graph'). If None, compute neighbors with scanpy.
    suffix : str
        Optional suffix for smoothed scores. By default '_kNN' is appended.

    Returns
    -------
    None. Adds smoothed scores to adata.obs.
    """
    if not (0 < decay < 1):
        raise ValueError("decay must be between 0 and 1")

    n_cells = adata.n_obs

    # Decide which graph to use
    if graph_key is None:
        # compute neighbors if not already done
        sc.pp.neighbors(adata, n_neighbors=k, use_rep=use_rep)
        graph_key = "connectivities"  # scanpy default
    if graph_key not in adata.obsp:
        raise ValueError(f"Graph '{graph_key}' not found in adata.obsp")

    # Get adjacency
    C = adata.obsp[graph_key]
    if not sparse.issparse(C):
        C = sparse.csr_matrix(C)

    for col in obs_columns:
        x = adata.obs[col].to_numpy(dtype=np.float32)
        smoothed = np.zeros(n_cells, dtype=np.float32)
        for i in range(n_cells):
            neighbors = C[i].indices
            conn = C[i].data
            # order neighbors by descending connectivity
            order = np.argsort(conn)[::-1]
            neighbors_ordered = neighbors[order]
            # include self first
            all_idx = np.insert(neighbors_ordered, 0, i)
            w = (1 - decay) ** np.arange(len(all_idx))
            w = w / w.sum()
            smoothed[i] = (x[all_idx] * w).sum()
        if up_only:
            smoothed = np.maximum(smoothed, x)
        adata.obs[f"{col}{suffix}"] = smoothed
