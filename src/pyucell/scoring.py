import numpy as np
from anndata import AnnData
from joblib import Parallel, delayed
from scipy import sparse

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


def _prepare_sig_indices(signatures: dict, genes: np.ndarray, missing_genes: str = "impute"):
    """
    Map signature genes to indices in the dataset.

    Parameters
    ----------
    signatures : dict
        Dictionary of signature_name -> list of genes.
    genes : np.ndarray
        Array of gene names in adata.var_names.
    missing_genes : str
        "impute": missing genes get a placeholder -1 (to be treated as max_rank)
        "skip": missing genes are simply removed

    Returns
    -------
    sig_indices : dict
        signature_name -> list of indices (or -1 for missing genes if impute)
    """
    gene_idx = {g: i for i, g in enumerate(genes)}
    sig_indices = {}

    gene_idx = {g: i for i, g in enumerate(genes)}
    sig_indices = {}

    for sig_name, sig_genes in signatures.items():
        pos_genes, neg_genes = _parse_sig(sig_genes)

        if missing_genes == "impute":
            pos_idx = [gene_idx.get(g, -1) for g in pos_genes]
            neg_idx = [gene_idx.get(g, -1) for g in neg_genes]
        elif missing_genes == "skip":
            pos_idx = [gene_idx[g] for g in pos_genes if g in gene_idx]
            neg_idx = [gene_idx[g] for g in neg_genes if g in gene_idx]
        else:
            raise ValueError("missing_genes must be 'impute' or 'skip'")

        sig_indices[sig_name] = {"pos": pos_idx, "neg": neg_idx}

    return sig_indices


def _calculate_U(ranks, idx, max_rank: int = 1500):
    idx = np.array(idx)
    lgt = len(idx)
    n_cells = ranks.shape[1]

    # Split indices (missing genes get index -1)
    missing_idx = idx[idx == -1]
    present_idx = idx[idx != -1]

    # Start with sum from missing genes, if any
    rank_sum = np.full(n_cells, len(missing_idx) * max_rank, dtype=np.float32)

    if len(present_idx) > 0:
        present_ranks = ranks[present_idx, :]
        # Always convert to dense safely
        present_ranks = present_ranks.toarray() if sparse.issparse(present_ranks) else np.asarray(present_ranks)
        # Ensure 2D shape even if single row
        if present_ranks.ndim == 1:
            present_ranks = present_ranks[np.newaxis, :]
        present_ranks = present_ranks.astype(np.float32)
        # rank==0 is equivalent to max_rank (for sparsity)
        present_ranks[present_ranks == 0] = max_rank
        rank_sum += present_ranks.sum(axis=0)

    s_min = lgt * (lgt + 1) / 2.0
    s_max = lgt * max_rank
    score = 1.0 - (rank_sum - s_min) / (s_max - s_min)
    return score


def _score_chunk(ranks: sparse.csr_matrix, sig_indices: dict, w_neg: float = 1.0, max_rank: int = 1500):
    n_genes, n_cells = ranks.shape
    n_signatures = len(sig_indices)
    scores = np.zeros((n_cells, n_signatures), dtype=np.float32)

    for j, (_sig_name, idx_dict) in enumerate(sig_indices.items()):
        pos_idx = idx_dict["pos"]
        neg_idx = idx_dict["neg"]

        pos_score = _calculate_U(ranks, pos_idx, max_rank=max_rank) if len(pos_idx) > 0 else 0.0
        neg_score = _calculate_U(ranks, neg_idx, max_rank=max_rank) if len(neg_idx) > 0 else 0.0

        score = pos_score - w_neg * neg_score
        score[score < 0] = 0.0  # clip negatives
        scores[:, j] = score

    return scores


def compute_ucell_scores(
    adata: AnnData,
    signatures: dict[str, list[str]],
    layer: str = None,
    max_rank: int = 1500,
    ties_method: str = "average",
    missing_genes: str = "impute",
    chunk_size: int = 500,
    w_neg: float = 1.0,
    suffix: str = "_UCell",
    n_jobs: int = -1,
):
    """
    Compute UCell scores for an AnnData object.

    Parameters
    ----------
    adata : AnnData
        An AnnData object (cells x genes)
    signatures:  Dict[str, List[str]]
        A dictionary of signatures, where the names of the entries are the signature names
    layer : str, optional
        Which layer to use (None = adata.X).
    max_rank : int, optional
        Cap ranks at this value (ranks > max_rank are dropped for sparsity).
    ties_method : str, optional
        Passed to scipy.stats.rankdata.
    missing_genes : str
        "impute": missing genes get a placeholder -1 (to be treated as max_rank)
        "skip": missing genes are simply removed
    chunk_size : int, optional
        The size of the blocks of cells to be processed at once. Avoids having large
        dense matrices in memory
    w_neg : float
        Weight on negative gene sets, when using signatures with positive and negative genes
    suffix : str, optional
        Suffix to append to column names in adata.obs.
    n_jobs : int, optional
        Number of parallel jobs

    Returns
    -------
    Adds signature scores in adata.obs

    """
    genes = adata.var_names.to_numpy()
    n_cells = adata.n_obs
    n_signatures = len(signatures)
    scores_all = np.zeros((n_cells, n_signatures), dtype=np.float32)

    # Precompute signature indices once
    sig_indices = _prepare_sig_indices(signatures, genes, missing_genes=missing_genes)

    # Split indices into chunks
    starts = list(range(0, n_cells, chunk_size))
    chunks = [(s, min(s + chunk_size, n_cells)) for s in starts]

    # Iterate over cell chunks
    def process_chunk(start, end):
        if layer:
            chunk_X = adata.layers[layer][start:end, :]
        else:
            chunk_X = adata.X[start:end, :]
        # compute ranks
        ranks_chunk = get_rankings(chunk_X, max_rank=max_rank, ties_method=ties_method)
        # get UCell scores for chunk
        scores_chunk = _score_chunk(ranks_chunk, sig_indices, w_neg=w_neg, max_rank=max_rank)
        return (start, end, scores_chunk)

    # Run chunks in serial or parallel
    if n_jobs == 1:
        results = [process_chunk(start, end) for start, end in chunks]
    else:
        results = Parallel(n_jobs=n_jobs, backend="loky")(delayed(process_chunk)(start, end) for start, end in chunks)

    # Merge results back
    scores_all = np.zeros((n_cells, n_signatures), dtype=np.float32)
    for start, end, scores_chunk in results:
        scores_all[start:end, :] = scores_chunk

    # Store scores in adata.obs with suffix
    for j, sig_name in enumerate(signatures.keys()):
        adata.obs[f"{sig_name}{suffix}"] = scores_all[:, j]
