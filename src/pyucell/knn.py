import numpy as np
import scanpy as sc
from scipy import sparse

from pyucell._torch_utils import import_torch, resolve_device


def smooth_knn_scores(
    adata,
    obs_columns,
    k=10,
    use_rep="X_pca",
    decay: float = 0.1,
    up_only: bool = False,
    graph_key=None,
    suffix="_kNN",
    device: str | None = None,
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
    device : str | None, optional
        ``None`` (default) → CPU sparse matmul.
        ``"cpu"`` / ``"cuda"`` / ``"mps"`` / ``"auto"`` → torch sparse matmul
        on that device. Requires the ``pyucell[gpu]`` extra.

    Returns
    -------
    None. Adds smoothed scores to adata.obs.
    """
    if not (0 < decay < 1):
        raise ValueError("decay must be between 0 and 1")

    n_cells = adata.n_obs

    if graph_key is None:
        sc.pp.neighbors(adata, n_neighbors=k, use_rep=use_rep)
        graph_key = "connectivities"
    if graph_key not in adata.obsp:
        raise ValueError(f"Graph '{graph_key}' not found in adata.obsp")

    C = adata.obsp[graph_key]
    if not sparse.issparse(C):
        C = sparse.csr_matrix(C)
    C = C.tocsr()

    W = _build_weight_matrix(C, decay=decay, n_cells=n_cells)

    # Stack score columns into (n_cells, n_cols) for a single matmul
    score_arr = np.stack(
        [adata.obs[col].to_numpy(dtype=np.float32) for col in obs_columns],
        axis=1,
    )

    if device is not None:
        torch = import_torch()
        dev = resolve_device(device)
        W_coo = W.tocoo()
        idx_t = torch.tensor(np.vstack([W_coo.row, W_coo.col]), dtype=torch.long, device=dev)
        val_t = torch.tensor(W_coo.data, dtype=torch.float32, device=dev)
        W_t = torch.sparse_coo_tensor(idx_t, val_t, size=W.shape).coalesce()
        X_t = torch.tensor(score_arr, dtype=torch.float32, device=dev)
        smoothed = torch.sparse.mm(W_t, X_t).cpu().numpy()
    else:
        smoothed = (W @ score_arr).astype(np.float32, copy=False)

    if up_only:
        smoothed = np.maximum(smoothed, score_arr)

    for j, col in enumerate(obs_columns):
        adata.obs[f"{col}{suffix}"] = smoothed[:, j]


def _build_weight_matrix(C: sparse.csr_matrix, decay: float, n_cells: int) -> sparse.csr_matrix:
    """Construct the row-normalised neighbour-weighting matrix W (n_cells x n_cells).

    Each row i has the original cell as the first entry (weight 1) and its
    neighbours from ``C[i, :]`` ordered by descending connectivity with
    weights ``(1 - decay)^r`` for r = 1..k. Rows are then normalised to sum to 1.
    """
    indptr = C.indptr
    indices = C.indices
    data = C.data
    n_nb_per_row = np.diff(indptr).astype(np.int64)
    n_total = int(n_nb_per_row.sum())

    if n_total == 0:
        # No edges anywhere: smoothing is the identity.
        return sparse.identity(n_cells, format="csr", dtype=np.float32)

    all_rows = np.repeat(np.arange(n_cells, dtype=np.int64), n_nb_per_row)
    # Sort by row asc, then by connectivity desc, jointly.
    key = np.lexsort((-data, all_rows))
    sorted_cols = indices[key]
    sorted_rows = all_rows[key]

    # Within-row 0-based position after sorting.
    row_starts = np.concatenate(([0], np.cumsum(n_nb_per_row)[:-1]))
    within_row_pos = np.arange(n_total, dtype=np.int64) - row_starts[sorted_rows]

    # Neighbour weight at within-row position p → (1-decay)^(p+1) (self gets p=0).
    nb_weights = (1.0 - decay) ** (within_row_pos + 1).astype(np.float64)

    self_rows = np.arange(n_cells, dtype=np.int64)
    self_cols = self_rows
    self_weights = np.ones(n_cells, dtype=np.float64)

    rows_all = np.concatenate([self_rows, sorted_rows])
    cols_all = np.concatenate([self_cols, sorted_cols])
    data_all = np.concatenate([self_weights, nb_weights]).astype(np.float32)

    W = sparse.csr_matrix((data_all, (rows_all, cols_all)), shape=(n_cells, n_cells))

    # Row-normalise.
    row_sums = np.asarray(W.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1.0
    inv = sparse.diags(1.0 / row_sums)
    return (inv @ W).astype(np.float32)
