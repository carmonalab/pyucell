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
    chunk_cells: int = 200,
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
    if is_sparse:
        Xcsc = X.tocsc()
    else:
        Xarr = np.asarray(X)

    data, rows, cols = [], [], []

    for cstart in range(0, n_cells, chunk_cells):
        cend = min(n_cells, cstart + chunk_cells)
        if is_sparse:
            chunk = Xcsc[cstart:cend, :].T.toarray()  # genes x chunk
        else:
            chunk = Xarr[cstart:cend, :].T  # genes x chunk
        # now chunk is genes x cells_in_chunk
        for j in range(chunk.shape[1]):
            col = chunk[:, j].astype(float)
            col[np.isnan(col)] = -np.inf
            ranks = rankdata(-col, method=ties_method)
            ranks_clipped = np.minimum(ranks, max_rank).astype(np.int32)
            idx = np.nonzero(ranks_clipped > 0)[0]
            data.extend(ranks_clipped[idx])
            rows.extend(idx)
            cols.extend([cstart + j] * len(idx))

    ranks_mat = sparse.coo_matrix(
        (np.array(data, dtype=np.int32), (np.array(rows, dtype=np.int32), np.array(cols, dtype=np.int32))),
        shape=(n_genes, n_cells),
    ).tocsr()
    return ranks_mat


def get_signatures_ucell(
    adata: AnnData,
    signatures: dict,
    ranks: sparse.csr_matrix = None,
    layer: str = None,
    max_rank: int = 1500,
    w_neg: float = 1.0,
    prefix: str = "_UCell",
) -> AnnData:
    """
    Compute UCell scores for AnnData.

    Parameters
    ----------
    - adata: AnnData with cells x genes
    - signatures: dict {name: [genes]} (genes may end with + or -)
    - ranks: precomputed ranks from store_rankings_ucell (optional)
    - layer: use adata.layers[layer] instead of .X
    Adds new columns to adata.obs named '<sig><prefix>'.
    """
    if ranks is None:
        ranks = get_rankings(adata, layer=layer, max_rank=max_rank)

    gene_names = list(adata.var_names)
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    n_cells = adata.n_obs
    ranks_csr = ranks.tocsr()

    for sname, sig_genes in signatures.items():
        pos_genes, neg_genes = _parse_sig(sig_genes)
        pos_idx = [gene_to_idx[g] for g in pos_genes if g in gene_to_idx]
        neg_idx = [gene_to_idx[g] for g in neg_genes if g in gene_to_idx]

        if not pos_idx and not neg_idx:
            warn(f"Signature '{sname}' has no matching genes.", stacklevel=2)

        # positive part
        if pos_idx:
            pos_score = _calculate_U(ranks=ranks_csr, idx=pos_idx, max_rank=max_rank)
        else:
            pos_score = np.full(n_cells, np.nan)

        # negative part
        if neg_idx:
            neg_score = _calculate_U(ranks=ranks_csr, idx=neg_idx, max_rank=max_rank)
            uscore = pos_score - w_neg * neg_score
        else:
            uscore = pos_score

        adata.obs[f"{sname}{prefix}"] = uscore

    return adata
