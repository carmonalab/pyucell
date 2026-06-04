import numpy as np
import pandas as pd
from anndata import AnnData
from joblib import Parallel, delayed
from scipy import sparse

from pyucell._torch_utils import import_torch, resolve_device
from pyucell.ranks import _rankings_torch, get_rankings


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
    if s_max == s_min:
        return np.zeros(n_cells, dtype=np.float32)
    score = 1.0 - (rank_sum - s_min) / (s_max - s_min)
    return score


def _score_chunk(ranks: sparse.csr_matrix, sig_indices: dict, w_neg: float = 1.0, max_rank: int = 1500):
    n_genes, n_cells = ranks.shape
    n_signatures = len(sig_indices)
    scores = np.zeros((n_cells, n_signatures), dtype=np.float32)

    for j, (_sig_name, idx_dict) in enumerate(sig_indices.items()):
        pos_idx = idx_dict["pos"]
        neg_idx = idx_dict["neg"]

        pos_score = _calculate_U(ranks, pos_idx, max_rank=max_rank) if len(pos_idx) > 0 else np.zeros(n_cells, dtype=np.float32)
        neg_score = _calculate_U(ranks, neg_idx, max_rank=max_rank) if len(neg_idx) > 0 else np.zeros(n_cells, dtype=np.float32)

        score = pos_score - w_neg * neg_score
        score[score < 0] = 0.0  # clip negatives
        scores[:, j] = score

    return scores


def _calculate_U_torch(ranks_dense, idx_tensor, n_missing: int, max_rank: int):  # pragma: no cover
    """Torch equivalent of ``_calculate_U`` operating on a dense (n_genes, n_cells) tensor.

    ``idx_tensor`` already excludes the ``-1`` placeholders; ``n_missing`` is
    their count.
    """
    torch = import_torch()
    lgt = idx_tensor.numel() + n_missing
    n_cells = ranks_dense.shape[1]
    device = ranks_dense.device

    rank_sum = torch.full((n_cells,), float(n_missing * max_rank), dtype=torch.float32, device=device)
    if idx_tensor.numel() > 0:
        present = ranks_dense.index_select(0, idx_tensor).to(torch.float32)
        max_rank_t = torch.tensor(float(max_rank), dtype=torch.float32, device=device)
        present = torch.where(present == 0, max_rank_t, present)
        rank_sum = rank_sum + present.sum(dim=0)

    s_min = lgt * (lgt + 1) / 2.0
    s_max = lgt * max_rank
    return 1.0 - (rank_sum - s_min) / (s_max - s_min)


def _score_chunk_torch(ranks_dense, sig_index_tensors, w_neg: float, max_rank: int):  # pragma: no cover
    torch = import_torch()
    n_cells = ranks_dense.shape[1]
    n_signatures = len(sig_index_tensors)
    device = ranks_dense.device
    scores = torch.zeros((n_cells, n_signatures), dtype=torch.float32, device=device)

    for j, (_sig_name, idx_dict) in enumerate(sig_index_tensors.items()):
        pos_idx, pos_missing = idx_dict["pos"]
        neg_idx, neg_missing = idx_dict["neg"]

        pos_score = (
            _calculate_U_torch(ranks_dense, pos_idx, pos_missing, max_rank)
            if (pos_idx.numel() + pos_missing) > 0
            else None
        )
        neg_score = (
            _calculate_U_torch(ranks_dense, neg_idx, neg_missing, max_rank)
            if (neg_idx.numel() + neg_missing) > 0
            else None
        )

        if pos_score is None and neg_score is None:
            continue
        if pos_score is None:
            score = -w_neg * neg_score
        elif neg_score is None:
            score = pos_score
        else:
            score = pos_score - w_neg * neg_score

        scores[:, j] = torch.clamp(score, min=0.0)

    return scores


def _build_sig_index_tensors(sig_indices, device):  # pragma: no cover
    """Move signature indices onto ``device``; split into (present_idx_tensor, n_missing)."""
    torch = import_torch()
    out = {}
    for sig_name, idx_dict in sig_indices.items():
        out[sig_name] = {}
        for key in ("pos", "neg"):
            idx_arr = np.asarray(idx_dict[key], dtype=np.int64)
            n_missing = int((idx_arr == -1).sum())
            present = idx_arr[idx_arr != -1]
            out[sig_name][key] = (
                torch.as_tensor(present, dtype=torch.long, device=device),
                n_missing,
            )
    return out


def compute_ucell_scores(
    adata: AnnData,
    signatures: dict[str, list[str]],
    layer: str = None,
    max_rank: int = 1500,
    ties_method: str = "average",
    missing_genes: str = "impute",
    chunk_size: int | None = None,
    w_neg: float = 1.0,
    suffix: str = "_UCell",
    n_jobs: int = -1,
    device: str | None = "cpu",
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
        Passed to scipy.stats.rankdata on the CPU path. The torch backend
        (``device != None``) only supports ``"min"`` and ``"ordinal"``.
    missing_genes : str
        "impute": missing genes get a placeholder -1 (to be treated as max_rank)
        "skip": missing genes are simply removed
    chunk_size : int, optional
        Size of the cell blocks processed at once. Defaults to 500 on CPU and
        5000 when ``device`` is set (GPUs benefit from larger batches).
    w_neg : float
        Weight on negative gene sets, when using signatures with positive and negative genes
    suffix : str, optional
        Suffix to append to column names in adata.obs.
    n_jobs : int, optional
        Number of parallel jobs (ignored when ``device`` is set; GPU chunks
        run serially to avoid multi-process CUDA init).
    device : str | None, optional
        ``"cpu"`` or ``None`` (default): CPU path with joblib parallelism (avoids PyTorch).
        ``"cuda"`` / ``"mps"`` / ``"auto"``: PyTorch backend for hardware acceleration.
        Requires the ``pyucell[gpu]`` extra.

    Returns
    -------
    Adds signature scores in adata.obs

    """
    genes = adata.var_names.to_numpy()
    n_cells = adata.n_obs
    n_signatures = len(signatures)
    scores_all = np.zeros((n_cells, n_signatures), dtype=np.float32)

    use_torch = device is not None and device not in ("cpu", "CPU")

    sig_indices = _prepare_sig_indices(signatures, genes, missing_genes=missing_genes)

    if chunk_size is None:
        chunk_size = 5000 if use_torch else 500

    starts = list(range(0, n_cells, chunk_size))
    chunks = [(s, min(s + chunk_size, n_cells)) for s in starts]

    if use_torch:  # pragma: no cover
        dev = resolve_device(device)
        sig_index_tensors = _build_sig_index_tensors(sig_indices, dev)

        for start, end in chunks:
            chunk_X = adata.layers[layer][start:end, :] if layer else adata.X[start:end, :]
            ranks_dense = _rankings_torch(chunk_X, max_rank=max_rank, ties_method=ties_method, device=dev)
            scores_chunk = _score_chunk_torch(ranks_dense, sig_index_tensors, w_neg=w_neg, max_rank=max_rank)
            scores_all[start:end, :] = scores_chunk.cpu().numpy()
    else:

        def process_chunk(start, end):
            chunk_X = adata.layers[layer][start:end, :] if layer else adata.X[start:end, :]
            ranks_chunk = get_rankings(chunk_X, max_rank=max_rank, ties_method=ties_method)
            scores_chunk = _score_chunk(ranks_chunk, sig_indices, w_neg=w_neg, max_rank=max_rank)
            return (start, end, scores_chunk)

        if n_jobs == 1:
            results = [process_chunk(start, end) for start, end in chunks]
        else:
            results = Parallel(n_jobs=n_jobs, backend="loky")(
                delayed(process_chunk)(start, end) for start, end in chunks
            )

        for start, end, scores_chunk in results:
            scores_all[start:end, :] = scores_chunk

    cols = [f"{sig}{suffix}" for sig in signatures.keys()]
    scores_df = pd.DataFrame(scores_all, index=adata.obs_names, columns=cols)

    # Drop columns from adata.obs that are about to be overwritten
    adata.obs = adata.obs.drop(columns=scores_df.columns, errors="ignore")
    # Concatenate into adata.obs
    adata.obs = pd.concat([adata.obs, scores_df], axis=1)
