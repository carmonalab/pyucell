import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
from scipy.stats import rankdata

from pyucell._torch_utils import import_torch, resolve_device, to_torch_dense


def get_rankings(
    data,
    layer: str = None,
    max_rank: int = 1500,
    ties_method: str = "average",
    device: str | None = None,
) -> AnnData:
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
        Passed to scipy.stats.rankdata on the CPU path. The torch backend
        (``device != None``) only supports ``"min"`` and ``"ordinal"``.
    device : str | None, optional
        ``None`` (default) → CPU path using ``scipy.stats.rankdata``.
        ``"cpu"`` / ``"cuda"`` / ``"mps"`` / ``"auto"`` → PyTorch backend on
        that device. Requires the ``pyucell[gpu]`` extra.

    Returns
    -------
    ranks_adata : AnnData
        An AnnData object of shape (n_cells, n_genes) containing the sparse matrix of ranks.
    """
    # Load matrix and preserve index names if available
    if isinstance(data, AnnData):
        X = data.layers[layer] if layer else data.X
        obs_names = data.obs_names
        var_names = data.var_names
    else:
        X = data
        obs_names = None
        var_names = None

    if device is not None:  # pragma: no cover
        dev = resolve_device(device)
        ranks_dense = _rankings_torch(X, max_rank=max_rank, ties_method=ties_method, device=dev)
        ranks_np = ranks_dense.cpu().numpy()
        ranks_mat = sparse.csr_matrix(ranks_np, dtype=np.int32)
    else:
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

        # All zeros fallback: Build empty sparse matrix safely
        if not data_parts:
            ranks_mat = sparse.csr_matrix((n_genes, n_cells), dtype=np.int32)
        else:
            # Concatenate arrays only once at the end
            data_arr = np.concatenate(data_parts).astype(np.int32)
            rows_arr = np.concatenate(row_parts).astype(np.int32)
            cols_arr = np.concatenate(col_parts).astype(np.int32)
            ranks_mat = sparse.csr_matrix((data_arr, (rows_arr, cols_arr)), shape=(n_genes, n_cells))

    # Wrap into standard AnnData (cells x genes)
    ranks_adata = AnnData(
        X=ranks_mat.T,
        obs=pd.DataFrame(index=obs_names) if obs_names is not None else None,
        var=pd.DataFrame(index=var_names) if var_names is not None else None
    )

    return ranks_adata


def _rankings_torch(X, max_rank: int, ties_method: str, device):  # pragma: no cover
    """Compute the (n_genes, n_cells) rank matrix as a dense torch tensor.

    Returns a dense int32 tensor on ``device``. Zero-valued (and capped)
    entries are represented as 0, matching the CSR convention used downstream.
    """
    torch = import_torch()

    if ties_method not in ("min", "ordinal"):
        raise ValueError(
            f"ties_method={ties_method!r} is not supported with device={device}. "
            "Use ties_method='min' or 'ordinal' for the torch backend "
            "(scipy's 'average' is not yet implemented on GPU)."
        )

    Xd = to_torch_dense(X, device=device, dtype=torch.float32)  # (n_cells, n_genes)
    n_cells, n_genes = Xd.shape

    neg_inf = torch.tensor(float("-inf"), dtype=Xd.dtype, device=device)
    masked = torch.where(Xd > 0, Xd, neg_inf)

    # Sort descending. argsort gives ordinal ranks directly; we add a 'min'
    # adjustment in a second pass for tied groups.
    order = torch.argsort(masked, dim=1, descending=True, stable=True)  # (n_cells, n_genes)

    arange_g = torch.arange(1, n_genes + 1, dtype=torch.int32, device=device).expand(n_cells, -1)

    if ties_method == "ordinal":
        ordinal_ranks = arange_g
    else:  # 'min'
        sorted_vals, _ = torch.sort(masked, dim=1, descending=True, stable=True)
        is_new = torch.ones_like(sorted_vals, dtype=torch.bool)
        is_new[:, 1:] = sorted_vals[:, 1:] != sorted_vals[:, :-1]
        positions = torch.arange(1, n_genes + 1, dtype=torch.int32, device=device).expand(n_cells, -1)
        start_pos = positions * is_new.to(torch.int32)
        ordinal_ranks, _ = torch.cummax(start_pos, dim=1)

    ranks = torch.empty((n_cells, n_genes), dtype=torch.int32, device=device)
    ranks.scatter_(1, order, ordinal_ranks)

    # Drop zeros (sparse marker) and anything beyond max_rank
    zero_mask = Xd == 0
    ranks = torch.where(zero_mask | (ranks > max_rank), torch.zeros_like(ranks), ranks)

    # Return as (n_genes, n_cells) to match the rest of the pipeline.
    return ranks.T.contiguous()