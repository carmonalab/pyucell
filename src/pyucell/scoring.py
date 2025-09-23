from warnings import warn
from anndata import AnnData
from scipy import sparse
from scipy.stats import rankdata
import numpy as np
from typing import Dict, List
from pyucell.ranks import get_rankings

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

def _prepare_sig_indices(signatures: dict, genes: np.ndarray):
    """
    Parse all signatures once and map to gene indices.

    Returns:
        sig_indices: dict of signature_name -> list of gene indices
    """
    gene_idx = {g: i for i, g in enumerate(genes)}
    sig_indices = {}

    for sig_name, sig_genes in signatures.items():
        pos_genes, _ = _parse_sig(sig_genes)
        idx = [gene_idx[g] for g in pos_genes if g in gene_idx]
        sig_indices[sig_name] = idx

    return sig_indices

def _calculate_U(ranks, idx, max_rank: int = 1500):

    if sparse.issparse(ranks):
        ranks_dense = ranks[idx, :].toarray()
    else:
        ranks_dense = np.asarray(ranks[idx, :])

    # Replace zeros (sparse) with max_rank
    ranks_dense[ranks_dense == 0] = max_rank
 
    # Sum ranks per cell
    rank_sum = np.array(ranks_dense.sum(axis=0)).ravel()
    lgt = len(idx)

    s_min = lgt * (lgt + 1) / 2.0
    s_max = lgt * max_rank
    U = rank_sum - s_min
    Umax = s_max - s_min
    score = 1.0 - (U / Umax)
    return score


def _score_chunk(ranks: sparse.csr_matrix,
    sig_indices: dict,
    max_rank: int = 1500):

    n_genes, n_cells = ranks.shape
    n_signatures = len(sig_indices)
    scores = np.zeros((n_cells, n_signatures), dtype=np.float32)

    for j, (sig_name, idx) in enumerate(sig_indices.items()):
        if len(idx) == 0:
            continue
        scores[:, j] = _calculate_U(ranks, idx, max_rank=max_rank)

    return scores


def compute_ucell_scores(
    adata: AnnData,
    signatures: Dict[str, List[str]],
    layer: str = None,
    max_rank: int = 1500,
    ties_method: str = "average",
    chunk_size: int = 500,
) -> AnnData:
    """
    Compute UCell scores for an AnnData object using cell-wise chunking.
    Stores results in adata.obs.
    """
    genes = adata.var_names.to_numpy()
    n_cells = adata.n_obs
    n_signatures = len(signatures)
    scores_all = np.zeros((n_cells, n_signatures), dtype=np.float32)

    # Precompute signature indices once
    sig_indices = _prepare_sig_indices(signatures, genes)

    # iterate over cell chunks
    for start in range(0, n_cells, chunk_size):
        end = min(start + chunk_size, n_cells)
        print(f"Processing cells {start}-{end}...")
        if layer:
            chunk_X = adata.layers[layer][start:end, :]
        else:
            chunk_X = adata.X[start:end, :]

        # compute ranks
        ranks_chunk = get_rankings(chunk_X, max_rank=max_rank, ties_method=ties_method)
        # compute UCell scores
        scores_chunk = _score_chunk(ranks_chunk, sig_indices, max_rank=max_rank)
        scores_all[start:end, :] = scores_chunk

    # store in adata.obs
    for j, sig_name in enumerate(signatures.keys()):
        adata.obs[sig_name] = scores_all[:, j]

    return adata
