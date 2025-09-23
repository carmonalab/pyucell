from warnings import warn

import numpy as np
from anndata import AnnData
from scipy import sparse
from scipy.stats import rankdata


def _parse_sig(sig):
    pos, neg = [], []
    for g in sig:
        gs = str(g).strip()
        if gs.endswith("+"):
            pos.append(gs[:-1])
        elif gs.endswith("-"):
            neg.append(gs[:-1])
        else:
            pos.append(gs)
    return pos, neg


def _calculate_U(ranks, idx, max_rank: int = 1500):
    sum_pos = np.array(ranks[idx, :].sum(axis=0)).ravel()
    npos = len(idx)
    s_min = npos * (npos + 1) / 2.0
    s_max = npos * max_rank
    U = sum_pos - s_min
    Umax = s_max - s_min
    score = 1.0 - (U / Umax)
    return score


def get_rankings(
    adata: AnnData,
    layer: str = None,
    max_rank: int = 1500,
    ties_method: str = "average",
) -> sparse.csr_matrix:
    """
    Compute per-cell ranks of genes for an AnnData object.

    Parameters
    ----------
    - adata: AnnData with cells x genes.
    - layer: which layer to use (None = adata.X)
    - max_rank: cap ranks at this value
    Returns sparse matrix (genes x cells) of ranks.
    """
    X = adata.layers[layer] if layer else adata.X
    n_cells, n_genes = X.shape

    # Convert to csc for fast column slicing
    is_sparse = sparse.issparse(X)
    Xarr = X.toarray() if is_sparse else np.asarray(X)

    # Allocate vectors, at most max_rank entries per cell
    n_cells, n_genes = X.shape
    nnz_per_cell = max_rank 
    nnz_total = n_cells * nnz_per_cell

    data = np.empty(nnz_total, dtype=np.int32)
    rows = np.empty(nnz_total, dtype=np.int32)
    cols = np.empty(nnz_total, dtype=np.int32)

    #Calculate ranks, while keeping the matrix sparse
    ptr = 0
    for j in range(n_cells):
        col = Xarr[j, :].astype(float)
        col[np.isnan(col)] = -np.inf
        ranks = rankdata(-col, method=ties_method)
        mask = ranks <= max_rank  #mask out ranks to impose sparsity
        idx = np.nonzero(mask)[0]
        rks = ranks[idx].astype(np.int32)
        n = len(idx)

        data[ptr:ptr+n] = rks
        rows[ptr:ptr+n] = idx
        cols[ptr:ptr+n] = j
        ptr += n

    # slice arrays to actual size
    data = data[:ptr]
    rows = rows[:ptr]
    cols = cols[:ptr]

    ranks_mat = sparse.coo_matrix((data, (rows,cols)), shape=(n_genes,n_cells)).tocsr()
    
    return ranks_mat
