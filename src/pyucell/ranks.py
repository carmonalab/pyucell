import numpy as np
from anndata import AnnData
from scipy import sparse
from scipy.stats import rankdata


def get_rankings(
    data,
    layer: str = None,
    max_rank: int = 1500,
    ties_method: str = "average",
) -> sparse.csr_matrix:
    """
    Compute per-cell ranks of genes for an AnnData object.

    Parameters
    ----------
    data : AnnData | np.ndarray | sparse matrix
        Either an AnnData object (cells x genes) or directly a 2D matrix.
    layer : str, optional
        Only used if input is AnnData. Which layer to use (None = adata.X).
    max_rank : int, optional
        Cap ranks at this value (ranks > max_rank are dropped for sparsity).
    ties_method : str, optional
        Passed to scipy.stats.rankdata.

    Returns
    -------
    ranks : csr_matrix of shape (genes, cells)
        Sparse matrix of ranks.
    """
    # Load matrix
    if isinstance(data, AnnData):
        X = data.layers[layer] if layer else data.X
    else:
        X = data

    n_cells, n_genes = X.shape

    # Store COO components per cell in lists of arrays
    data_parts = []
    row_parts = []
    col_parts = []

    for j in range(n_cells):
        col = X[j, :]
        if sparse.issparse(col):
            col = col.toarray().ravel()
        else:
            col = np.asarray(col, dtype=float)

        # missing values
        np.nan_to_num(col, copy=False)

        # Only rank non-zero elements
        nz_idx = np.nonzero(col)[0]
        if len(nz_idx) == 0:
            continue

        nz_vals = col[nz_idx]
        ranks = rankdata(-nz_vals, method=ties_method).astype(np.int32)

        keep_mask = ranks <= max_rank
        kept_idx = nz_idx[keep_mask]
        kept_ranks = ranks[keep_mask]

        if len(kept_idx) > max_rank:
            kept_idx = kept_idx[:max_rank]
            kept_ranks = kept_ranks[:max_rank]

        n = len(kept_idx)
        if n == 0:
            continue

        # Convert to small NumPy arrays per cell
        data_parts.append(kept_ranks)
        row_parts.append(kept_idx)
        col_parts.append(np.full(n, j, dtype=np.int32))

    # All zeros
    if not data_parts:
        return sparse.csr_matrix((n_genes, n_cells), dtype=np.int32)

    # Concatenate arrays only once at the end
    data_arr = np.concatenate(data_parts).astype(np.int32)
    rows_arr = np.concatenate(row_parts).astype(np.int32)
    cols_arr = np.concatenate(col_parts).astype(np.int32)

    ranks_mat = sparse.csr_matrix((data_arr, (rows_arr, cols_arr)), shape=(n_genes, n_cells))
    return ranks_mat
